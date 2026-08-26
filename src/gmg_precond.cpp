/*! \file gmg_precond.cpp
 *  \brief Parallel geometric multigrid preconditioner for sparse matrices
 *         arising from structured-grid finite-difference discretizations
 *         (companion to rbsor_precond.cpp).
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 */

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <omp.h>
#include "gmg_precond.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


/* Portable CAS-loop atomic add for double (std::atomic<double> has no
 * fetch_add() until C++20). Duplicated from the identical helper in
 * scharge.cpp rather than shared, to keep this file dependency-free --
 * same rationale as color_graph() being duplicated from
 * RBSOR_Precond rather than shared with it.
 */
namespace {
inline void atomic_add_double( double &target, double val )
{
    static_assert( sizeof(double) == sizeof(uint64_t), "unexpected double size" );
    uint64_t *addr = reinterpret_cast<uint64_t *>( &target );
    uint64_t old_bits = __atomic_load_n( addr, __ATOMIC_RELAXED );
    for( ;; ) {
        double old_val, new_val;
        std::memcpy( &old_val, &old_bits, sizeof(double) );
        new_val = old_val + val;
        uint64_t new_bits;
        std::memcpy( &new_bits, &new_val, sizeof(double) );
        if( __atomic_compare_exchange_n( addr, &old_bits, new_bits, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED ) )
            break;
    }
}
}


/* *****************************************************************************
 * Construction / bookkeeping
 */

GMG_Precond::GMG_Precond( uint32_t nx, uint32_t ny, uint32_t nz, uint32_t dim,
                           uint32_t nlevels, uint32_t npre, uint32_t npost,
                           uint32_t ncoarse_sweeps, double coarse_rtol, double w,
                           double coarse_w, uint32_t ncycles, uint32_t nthreads )
    : _nx0(nx), _ny0(ny), _nz0(nz), _dim(dim), _nlevels_req(nlevels),
      _npre(npre), _npost(npost), _ncoarse_sweeps(ncoarse_sweeps), _coarse_rtol(coarse_rtol),
      _w(w), _coarse_w(coarse_w == 0.0 ? w : coarse_w), _ncycles(ncycles), _nthreads(nthreads),
      _last_coarse_sweeps(0), _prepared(false)
{
    if( dim < 1 || dim > 3 )
        throw( ErrorDim( ERROR_LOCATION, "dim must be 1, 2 or 3" ) );
    if( nx < 1 || (dim >= 2 && ny < 1) || (dim >= 3 && nz < 1) )
        throw( ErrorDim( ERROR_LOCATION, "invalid grid size" ) );
    if( npre < 1 || npost < 1 )
        throw( ErrorDim( ERROR_LOCATION, "npre/npost must be at least 1" ) );
    if( ncoarse_sweeps < 1 )
        throw( ErrorDim( ERROR_LOCATION, "ncoarse_sweeps must be at least 1" ) );
    if( coarse_rtol <= 0.0 || coarse_rtol >= 1.0 )
        throw( ErrorDim( ERROR_LOCATION, "coarse_rtol must be in (0,1)" ) );
    if( w <= 0.0 || w >= 2.0 )
        throw( ErrorDim( ERROR_LOCATION, "relaxation factor w must be in (0,2)" ) );
    if( coarse_w < 0.0 || coarse_w >= 2.0 )
        throw( ErrorDim( ERROR_LOCATION, "relaxation factor coarse_w must be in [0,2), 0 meaning \"same as w\"" ) );
    if( ncycles < 1 )
        throw( ErrorDim( ERROR_LOCATION, "ncycles must be at least 1" ) );
}


GMG_Precond::~GMG_Precond()
{
}


GMG_Precond *GMG_Precond::copy( void ) const
{
    // Same philosophy as RBSOR_Precond::copy(): the hierarchy is cheap
    // to rebuild from the matrix at first use, so we hand back a
    // fresh, unprepared preconditioner with the same parameters rather
    // than deep-copy the (possibly large) level data.
    GMG_Precond *pc = new GMG_Precond( _nx0, _ny0, _nz0, _dim, _nlevels_req,
                                        _npre, _npost, _ncoarse_sweeps, _coarse_rtol, _w,
                                        _coarse_w, _ncycles, _nthreads );
    if( !_node_map.empty() )
        pc->set_node_map( _node_map );
    return( pc );
}


void GMG_Precond::clear( void )
{
    _prepared = false;
    _level.clear();
    _x.clear();
    _b.clear();
    _r.clear();
    _lvl0_diag_ref.clear();
    _lvl1_lin_val.clear();
    _delta_src.clear();
    _delta_dst.clear();
    _delta_coef.clear();
    _delta_group.clear();
}


std::string GMG_Precond::typestring( void ) const
{
    return( std::string( "GMG" ) );
}


/* *****************************************************************************
 * Hierarchy construction
 */

/* Refreshes _level[0]'s VALUES only, in place, reusing the ptr/col
 * layout already established by build_level0() in prepare(). Safe
 * because A's sparsity PATTERN is invariant for the lifetime of this
 * preconditioner (fixed by geometry/discretization), so only VALUES
 * need updating from one construct() call to the next -- e.g. every
 * Newton round and every major cycle, as the nonlinear plasma model's
 * Jacobian diagonal evolves with the current potential. Mirrors
 * build_level0()'s value-assignment loop exactly (same order), just
 * without touching ptr/col or re-pushing a level.
 */
void GMG_Precond::refresh_level0_values( const CRowMatrix &A )
{
    Level &lev = _level[0];
    int nfull = (int)_node_map.size();

    for( int m = 0; m < nfull; m++ ) {
        int32_t r = _node_map[m];
        int p = lev.ptr[m];
        if( r < 0 ) {
            lev.val[p] = 1.0;
        } else {
            int a0 = A.ptr(r), a1 = A.ptr(r+1);
            for( int a = a0; a < a1; a++, p++ )
                lev.val[p] = A.val(a);
        }
    }
}


/* Embeds matrix A (A.rows() == number of *free* nodes, "dof") into a
 * level-0 operator sized nx0*ny0*nz0 (the full mesh), using
 * _node_map/_row_to_node built by prepare(). Every eliminated
 * (fixed/Dirichlet) mesh slot becomes a trivial decoupled identity
 * row (diagonal 1, no off-diagonal, always solves to x=0); every free
 * mesh slot gets A's row for that node, with column indices remapped
 * from matrix-row space back to mesh-index space via _row_to_node so
 * that all of level 0's neighbours are addressed the same way as
 * every coarser level (by mesh index).
 */
void GMG_Precond::build_level0( const CRowMatrix &A )
{
    int nfull = (int)_node_map.size();

    Level lev;
    lev.nx = _nx0;
    lev.ny = _ny0;
    lev.nz = _nz0;
    lev.n  = nfull;

    lev.ptr.assign( nfull+1, 0 );
    for( int m = 0; m < nfull; m++ ) {
        int32_t r = _node_map[m];
        int cnt = (r < 0 ? 1 : A.ptr(r+1)-A.ptr(r));
        lev.ptr[m+1] = lev.ptr[m] + cnt;
    }

    lev.col.resize( lev.ptr.back() );
    lev.val.resize( lev.ptr.back() );

    for( int m = 0; m < nfull; m++ ) {
        int32_t r = _node_map[m];
        int p = lev.ptr[m];
        if( r < 0 ) {
            lev.col[p] = m;
            lev.val[p] = 1.0;
        } else {
            int a0 = A.ptr(r), a1 = A.ptr(r+1);
            for( int a = a0; a < a1; a++, p++ ) {
                lev.col[p] = _row_to_node[ A.col(a) ];
                lev.val[p] = A.val(a);
            }
        }
    }

    _level.push_back( std::move(lev) );
}


bool GMG_Precond::next_level_size( const Level &fine, uint32_t &cnx, uint32_t &cny, uint32_t &cnz ) const
{
    // Same coarsening convention as EpotMGSolver::prepare_mg_geom():
    // size[d] = (size[d]+1)/2, only exact (i.e. keeps the coarse grid
    // nodes coincident with every other fine grid node) when size[d]
    // is odd. Stop the hierarchy rather than throwing if that's not
    // the case, or if the grid would become too small to be useful.
    cnx = fine.nx; cny = fine.ny; cnz = fine.nz;

    if( fine.nx % 2 == 0 )
        return( false );
    cnx = (fine.nx+1)/2;

    if( _dim >= 2 ) {
        if( fine.ny % 2 == 0 )
            return( false );
        cny = (fine.ny+1)/2;
    }
    if( _dim >= 3 ) {
        if( fine.nz % 2 == 0 )
            return( false );
        cnz = (fine.nz+1)/2;
    }

    // Stop once the coarsest grid would have essentially no interior
    // left -- below this the SOR "coarse solve" converges in a
    // handful of sweeps anyway, so there's nothing to gain from
    // coarsening further, and the (n+1)/2 formula stops making sense
    // once a dimension collapses to 1.
    if( cnx < 3 || (_dim >= 2 && cny < 3) || (_dim >= 3 && cnz < 3) )
        return( false );

    return( true );
}


/* fine_node_map, when non-NULL, is _node_map (the level-0 mesh-node ->
 * matrix-row map, see set_node_map()): entries whose fine index maps
 * to an eliminated (fixed/Dirichlet, node_map<0) node are dropped
 * entirely from both P and R (no renormalization -- same convention
 * already used for out-of-range/boundary neighbours). This matters
 * for two reasons:
 *   - Without it, prolong_add() would write nonzero corrections into
 *     fixed mesh slots that the smoother otherwise correctly locks at
 *     x=0 forever, corrupting every free neighbour's stencil on the
 *     next smoothing sweep.
 *   - Without it, the Galerkin product A_coarse = R*A*P would pick up
 *     spurious coarse-coarse coupling "through" a fixed node: for a
 *     fixed row i (a pure identity row in A), AP[i][:] = P[i][:]
 *     rather than 0, since A[i][i]=1 is the only nonzero. Dropping
 *     fixed targets from P makes P's row i empty, so AP[i][:] = 0 as
 *     it physically should be (no coupling through solid/eliminated
 *     material).
 * Only the level-0 -> level-1 transfer needs this: coarser levels are
 * pure Galerkin unknowns with no elimination of their own.
 */
void GMG_Precond::build_transfer_operators( const Level &fine, Level &coarse,
                                             const std::vector<int32_t> *fine_node_map ) const
{
    int fnx = fine.nx,   fny = (_dim>=2? fine.ny  :1),   fnz = (_dim>=3? fine.nz  :1);
    int cnx = coarse.nx, cny = (_dim>=2? coarse.ny :1),  cnz = (_dim>=3? coarse.nz :1);

    // c = 1/2^dim, the standard variational full-weighting factor
    // relating restriction to the transpose of prolongation.
    double c = 1.0;
    for( uint32_t d = 0; d < _dim; d++ )
        c *= 0.5;

    // For every coarse node, gather the (fine index, weight) pairs of
    // the multilinear interpolation stencil (identical stencil to the
    // one implied by EpotMGSolver::prolong_3d(): weight 1 at the
    // coincident node, 1/2 at each face-adjacent neighbour, 1/4 at
    // each edge-adjacent neighbour, 1/8 at each corner neighbour, and
    // so on generalized to 1D/2D/3D). Out-of-range fine neighbours are
    // simply omitted (no renormalization), matching
    // EpotMGSolver::restrict_3d()'s boundary-face handling.
    std::vector< std::vector< std::pair<int,double> > > rows( coarse.n );

    for( int K = 0; K < cnz; K++ ) {
        for( int J = 0; J < cny; J++ ) {
            for( int I = 0; I < cnx; I++ ) {

                int crow = K*cnx*cny + J*cnx + I;
                std::vector< std::pair<int,double> > &entries = rows[crow];

                int dkmin = (_dim>=3? -1:0), dkmax = (_dim>=3? 1:0);
                int djmin = (_dim>=2? -1:0), djmax = (_dim>=2? 1:0);

                for( int dk = dkmin; dk <= dkmax; dk++ ) {
                    int fk = 2*K+dk;
                    if( fk < 0 || fk >= fnz ) continue;
                    double wk = (dk==0? 1.0 : 0.5);

                    for( int dj = djmin; dj <= djmax; dj++ ) {
                        int fj = 2*J+dj;
                        if( fj < 0 || fj >= fny ) continue;
                        double wj = (dj==0? 1.0 : 0.5);

                        for( int di = -1; di <= 1; di++ ) {
                            int fi = 2*I+di;
                            if( fi < 0 || fi >= fnx ) continue;
                            double wi = (di==0? 1.0 : 0.5);

                            double weight = wi*wj*wk;
                            int frow = fk*fnx*fny + fj*fnx + fi;
                            if( fine_node_map && (*fine_node_map)[frow] < 0 )
                                continue; // fine target is fixed/eliminated -- never interpolate into it
                            entries.push_back( std::make_pair(frow, weight) );
                        }
                    }
                }
            }
        }
    }

    // Assemble R = c * (the above), coarse.n rows x fine.n cols.
    coarse.r_ptr.assign( coarse.n+1, 0 );
    for( int i = 0; i < coarse.n; i++ )
        coarse.r_ptr[i+1] = coarse.r_ptr[i] + (int)rows[i].size();
    coarse.r_col.resize( coarse.r_ptr.back() );
    coarse.r_val.resize( coarse.r_ptr.back() );
    for( int i = 0; i < coarse.n; i++ ) {
        int p = coarse.r_ptr[i];
        for( size_t k = 0; k < rows[i].size(); k++, p++ ) {
            coarse.r_col[p] = rows[i][k].first;
            coarse.r_val[p] = c*rows[i][k].second;
        }
    }

    // Assemble P as the (unscaled) transpose, fine.n rows x coarse.n cols.
    std::vector<int> col_count( fine.n, 0 );
    for( int i = 0; i < coarse.n; i++ )
        for( size_t k = 0; k < rows[i].size(); k++ )
            col_count[ rows[i][k].first ]++;

    coarse.p_ptr.assign( fine.n+1, 0 );
    for( int i = 0; i < fine.n; i++ )
        coarse.p_ptr[i+1] = coarse.p_ptr[i] + col_count[i];
    coarse.p_col.resize( coarse.p_ptr.back() );
    coarse.p_val.resize( coarse.p_ptr.back() );

    std::vector<int> cursor( fine.n );
    for( int i = 0; i < fine.n; i++ )
        cursor[i] = coarse.p_ptr[i];

    for( int I = 0; I < coarse.n; I++ ) {
        for( size_t k = 0; k < rows[I].size(); k++ ) {
            int f = rows[I][k].first;
            double weight = rows[I][k].second;
            int p = cursor[f]++;
            coarse.p_col[p] = I;
            coarse.p_val[p] = weight;
        }
    }
}


/* Sparse-sparse product C = A*B (A: m x k CSR, B: k x n CSR), used to
 * compute the Galerkin coarse-grid operator. Parallelized over the
 * rows of A/C -- each output row is independent, so this is
 * embarrassingly parallel. Uses a thread-local hash accumulator; the
 * matrices involved here are the (small) per-level operators built
 * once in construct(), so simplicity is favoured over the last bit of
 * performance (a dense marker-array Gustavson accumulator would be
 * the next optimization if this ever shows up in a profile).
 */
static void spgemm( int m,
                     const std::vector<int> &ptrA, const std::vector<int> &colA, const std::vector<double> &valA,
                     const std::vector<int> &ptrB, const std::vector<int> &colB, const std::vector<double> &valB,
                     std::vector<int> &ptrC, std::vector<int> &colC, std::vector<double> &valC )
{
    std::vector< std::vector< std::pair<int,double> > > rowsC( m );

    #pragma omp parallel
    {
        std::unordered_map<int,double> acc;
        #pragma omp for schedule(dynamic,64)
        for( int i = 0; i < m; i++ ) {
            acc.clear();
            for( int pa = ptrA[i]; pa < ptrA[i+1]; pa++ ) {
                int k = colA[pa];
                double va = valA[pa];
                if( va == 0.0 )
                    continue;
                for( int pb = ptrB[k]; pb < ptrB[k+1]; pb++ )
                    acc[ colB[pb] ] += va*valB[pb];
            }

            std::vector< std::pair<int,double> > &row = rowsC[i];
            row.reserve( acc.size() );
            for( std::unordered_map<int,double>::const_iterator it = acc.begin(); it != acc.end(); ++it )
                row.push_back( std::make_pair(it->first, it->second) );
            std::sort( row.begin(), row.end() );
        }
    }

    ptrC.assign( m+1, 0 );
    for( int i = 0; i < m; i++ )
        ptrC[i+1] = ptrC[i] + (int)rowsC[i].size();
    colC.resize( ptrC.back() );
    valC.resize( ptrC.back() );
    for( int i = 0; i < m; i++ ) {
        int p = ptrC[i];
        for( size_t k = 0; k < rowsC[i].size(); k++, p++ ) {
            colC[p] = rowsC[i][k].first;
            valC[p] = rowsC[i][k].second;
        }
    }
}


void GMG_Precond::galerkin_coarsen( const Level &fine, Level &coarse ) const
{
    // AP = A_fine * P   (fine.n x coarse.n)
    std::vector<int> ap_ptr, ap_col;
    std::vector<double> ap_val;
    spgemm( fine.n, fine.ptr, fine.col, fine.val,
            coarse.p_ptr, coarse.p_col, coarse.p_val,
            ap_ptr, ap_col, ap_val );

    // A_coarse = R * AP   (coarse.n x coarse.n)
    spgemm( coarse.n, coarse.r_ptr, coarse.r_col, coarse.r_val,
            ap_ptr, ap_col, ap_val,
            coarse.ptr, coarse.col, coarse.val );
}


/* Greedy 2-coloring of the symmetrized sparsity graph, identical in
 * spirit to RBSOR_Precond::color_graph() (see that file for the full
 * rationale); duplicated here rather than shared so that GMG_Precond
 * has no compile-time dependency on RBSOR_Precond.
 */
void GMG_Precond::color_graph( Level &lev, const std::vector<int32_t> *node_map ) const
{
    /* node_map, when non-NULL, is _node_map: rows whose mesh slot is
     * eliminated (node_map < 0) are left OUT of the colour lists entirely,
     * so the smoother never visits them.
     *
     * They are inert by construction, not merely unimportant. build_level0()
     * gives an eliminated slot a single entry, the identity (m,m)=1, and
     * maps every live row's columns through _row_to_node, which can only
     * ever produce free nodes -- so nothing anywhere points AT an eliminated
     * row. solve() fills their rhs with 0. Relaxing one therefore computes
     *     x[m] <- (1-w)*x[m] + w*0
     * from x[m] = 0, i.e. it does nothing at all, once per sweep, four
     * sweeps per V-cycle, twice per BiCGSTAB iteration.
     *
     * They were previously visited because the existing skip in relax() is
     * on inv_diag == 0, and an identity row's diagonal is 1, not 0.
     *
     * Only level 0 has eliminated slots -- coarser levels are pure Galerkin
     * unknowns -- so node_map is passed there and NULL everywhere else. A
     * coarse row that ends up entirely decoupled still gets zero diagonal
     * and is still caught by relax()'s existing inv_diag test.
     */
    int n = lev.n;
    std::vector<uint8_t> color( n, 0xFF );

    for( int i = 0; i < n; i++ ) {
        if( node_map && (*node_map)[i] < 0 )
            continue;
        bool red_taken = false, black_taken = false;
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            int j = lev.col[p];
            if( j == i || j < 0 || j >= n )
                continue;
            if( color[j] == 0 ) red_taken = true;
            else if( color[j] == 1 ) black_taken = true;
        }
        if( !red_taken ) color[i] = 0;
        else if( !black_taken ) color[i] = 1;
        else color[i] = 2;
    }

    for( int i = 0; i < n; i++ ) {
        if( node_map && (*node_map)[i] < 0 ) continue;
        if( color[i] == 2 ) continue;
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            int j = lev.col[p];
            if( j == i || j < 0 || j >= n )
                continue;
            if( color[j] == color[i] )
                color[ std::max(i,j) ] = 2;
        }
    }

    lev.red.clear(); lev.black.clear(); lev.leftover.clear();
    for( int i = 0; i < n; i++ ) {
        if( node_map && (*node_map)[i] < 0 ) continue;
        if( color[i] == 0 ) lev.red.push_back( (uint32_t)i );
        else if( color[i] == 1 ) lev.black.push_back( (uint32_t)i );
        else lev.leftover.push_back( (uint32_t)i );
    }

    if( !lev.leftover.empty() ) {
        ibsimu.message( 2 ) << "GMG_Precond: " << lev.leftover.size()
                            << " of " << n
                            << " rows could not be 2-coloured at this level, "
                            << "processing them serially\n";
        ibsimu.flush();
    }
}


/* Stores the *raw* reciprocal diagonal 1/aii, not w/aii -- the
 * relaxation factor is applied explicitly by relax()/relax_serial()
 * instead, so that the same precomputed inv_diag can be reused with
 * different w for pre-/post-smoothing (w=1, plain Gauss-Seidel -- a
 * good high-frequency error *smoother*) versus the coarsest-level
 * solve (w>1 over-relaxation can give substantially faster asymptotic
 * *convergence*, which is what actually matters there; see
 * coarse_w in the constructor).
 */
void GMG_Precond::find_diagonals( Level &lev ) const
{
    int n = lev.n;
    lev.diag_idx.assign( n, -1 );
    lev.inv_diag.assign( n, 0.0 );
    for( int i = 0; i < n; i++ ) {
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            if( lev.col[p] == i ) {
                lev.diag_idx[i] = p;
                if( lev.val[p] != 0.0 )
                    lev.inv_diag[i] = 1.0 / lev.val[p];
                break;
            }
        }
    }
}


/* Precompute the scatter tables used by construct() to apply the
 * level0->level1 Galerkin update cheaply. Only level 0's diagonal
 * changes between construct() calls (see the class-level comment on
 * the corresponding members in gmg_precond.hpp), so writing
 * A0 = A0_ref + diag(delta) (A0_ref being whatever level 0 looked like
 * when this was last called from prepare(), delta the change since
 * then) gives, via the Galerkin definition A1 = R0*A0*P0:
 *
 *   A1 = R0*A0_ref*P0 + R0*diag(delta)*P0
 *      = A1_ref                + R0*diag(delta)*P0
 *
 * A1_ref is exactly the level-1 operator galerkin_coarsen() already
 * built in prepare() (cached below as _lvl1_lin_val). For the second
 * term, using R0 = c*P0^T (c = 1/2^dim, see build_transfer_operators()):
 *
 *   R0*diag(delta)*P0 = c * P0^T * diag(delta) * P0
 *
 * so its (I,J) entry is c * sum_m delta[m]*P0[m,I]*P0[m,J] -- a sum
 * over fine rows m, each contributing a small (<= (2^dim)^2) dense
 * outer product of its own P-row with itself. This function walks
 * every fine row once and, for every (I,J) pair its P-row touches,
 * looks up that entry's position in level 1's (fixed-pattern) CSR
 * storage, recording the (fine row, level1.val index, coefficient)
 * triple. construct() then just scatter-accumulates
 * delta[m]*coefficient into level1.val at those indices every call --
 * O(fine.n * average-stencil^2) cheap array writes instead of a
 * general hash-map-based sparse-sparse product over the largest
 * matrix in the hierarchy, every single Newton iteration.
 *
 * The lookup below assumes every (I,J) pair produced this way already
 * exists in level 1's pattern; that holds because A1_ref's own
 * pattern is the union of contributions from every nonzero of A0_ref
 * (off-diagonal *and* diagonal), and A0_ref's diagonal is a fully
 * generic (essentially never exactly zero) physical value at every
 * free node, so the diagonal-only pattern is already a subset of
 * A1_ref's pattern. If some pathological entry is ever missing anyway
 * (e.g. a genuinely zero reference diagonal at some node), the
 * contribution is just dropped rather than throwing: this is a
 * preconditioner, so the very worst case is a slightly less accurate
 * M^-1, never an incorrect BiCGSTAB/Newton result.
 */
void GMG_Precond::build_delta_update_tables( void )
{
    _lvl0_diag_ref.clear();
    _lvl1_lin_val.clear();
    _delta_src.clear();
    _delta_dst.clear();
    _delta_coef.clear();
    _delta_group.clear();

    if( _level.size() < 2 )
        return; // no coarser level exists to update this way

    const Level &lev0 = _level[0];
    const Level &lev1 = _level[1];

    _lvl0_diag_ref.assign( lev0.n, 0.0 );
    for( int m = 0; m < lev0.n; m++ ) {
        int di = lev0.diag_idx[m];
        if( di >= 0 )
            _lvl0_diag_ref[m] = lev0.val[di];
    }

    _lvl1_lin_val = lev1.val;

    double c = 1.0;
    for( uint32_t d = 0; d < _dim; d++ )
        c *= 0.5;

    int nfine = (int)lev1.p_ptr.size()-1;

    int nthr = omp_get_max_threads();
    std::vector< std::vector<int32_t> > src_buf( nthr );
    std::vector< std::vector<int32_t> > dst_buf( nthr );
    std::vector< std::vector<double> >  coef_buf( nthr );

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<int32_t> &src  = src_buf[tid];
        std::vector<int32_t> &dst  = dst_buf[tid];
        std::vector<double>  &coef = coef_buf[tid];

        #pragma omp for schedule(dynamic,256)
        for( int m = 0; m < nfine; m++ ) {
            int p0 = lev1.p_ptr[m], p1 = lev1.p_ptr[m+1];
            for( int a = p0; a < p1; a++ ) {
                int Ia = lev1.p_col[a];
                double wa = lev1.p_val[a];
                for( int b = p0; b < p1; b++ ) {
                    int Ib = lev1.p_col[b];
                    double wb = lev1.p_val[b];

                    int found = -1;
                    for( int p = lev1.ptr[Ia]; p < lev1.ptr[Ia+1]; p++ ) {
                        if( lev1.col[p] == Ib ) {
                            found = p;
                            break;
                        }
                    }
                    if( found < 0 )
                        continue; // see function comment -- defensive, should not happen

                    src.push_back( m );
                    dst.push_back( found );
                    coef.push_back( c*wa*wb );
                }
            }
        }
    }

    size_t total = 0;
    for( int t = 0; t < nthr; t++ )
        total += src_buf[t].size();
    _delta_src.reserve( total );
    _delta_dst.reserve( total );
    _delta_coef.reserve( total );
    for( int t = 0; t < nthr; t++ ) {
        _delta_src.insert( _delta_src.end(), src_buf[t].begin(), src_buf[t].end() );
        _delta_dst.insert( _delta_dst.end(), dst_buf[t].begin(), dst_buf[t].end() );
        _delta_coef.insert( _delta_coef.end(), coef_buf[t].begin(), coef_buf[t].end() );
    }

    /* GROUP THE TRIPLES BY DESTINATION.
     *
     * The apply loop used to scatter these with atomic_add_double(), which
     * is correct but NOT reproducible: floating-point addition is not
     * associative, so the value landing in lev1.val depended on the order
     * the threads happened to arrive in. That made the preconditioner --
     * and therefore epot -- differ in its last bits between two runs of the
     * same binary on the same input.
     *
     * That is not a harmless last-bit difference. It is the first link in a
     * chain: the trajectory integrator's step controller drives the local
     * error estimate to sit AT its tolerance by construction, so an
     * arbitrarily small change in the field flips accept/reject decisions.
     * Measured on the ECR cut case, two runs of one binary differed by
     * hundreds of ODE steps out of two million in cycle 0 alone, and the
     * iteration then amplified that over cycles.
     *
     * Sorting by destination lets each destination be summed by exactly one
     * thread, in a fixed index order. Deterministic, and faster -- the CAS
     * loop is gone.
     *
     * The sort key is the FULL TRIPLE (dst, src, coef), not dst alone.
     * An earlier version used stable_sort on dst only, reasoning that the
     * build order was deterministic because the per-thread buffers are
     * concatenated in thread order. That reasoning was wrong: the build
     * loop above is schedule(dynamic,256), so WHICH rows land in which
     * thread's buffer varies from run to run, and a stable sort faithfully
     * preserves that varying order within each destination group. The
     * summation order therefore still varied, and so did the result -- the
     * bug survived the fix that was supposed to remove it.
     *
     * Ordering on the whole triple is canonical regardless of how the
     * triples were generated, so it is robust to the schedule clause rather
     * than dependent on it. Any triples identical in all three fields are
     * interchangeable, so ties among them cannot affect the sum.
     */
    const size_t ntot = _delta_dst.size();
    std::vector<size_t> perm( ntot );
    std::iota( perm.begin(), perm.end(), (size_t)0 );
    std::sort( perm.begin(), perm.end(),
               [&]( size_t a, size_t b ) {
                   if( _delta_dst[a] != _delta_dst[b] )
                       return( _delta_dst[a] < _delta_dst[b] );
                   if( _delta_src[a] != _delta_src[b] )
                       return( _delta_src[a] < _delta_src[b] );
                   return( _delta_coef[a] < _delta_coef[b] );
               } );

    std::vector<int32_t> s2( ntot ), d2( ntot );
    std::vector<double>  c2( ntot );
    for( size_t i = 0; i < ntot; i++ ) {
        s2[i] = _delta_src[ perm[i] ];
        d2[i] = _delta_dst[ perm[i] ];
        c2[i] = _delta_coef[ perm[i] ];
    }
    _delta_src.swap( s2 );
    _delta_dst.swap( d2 );
    _delta_coef.swap( c2 );

    _delta_group.clear();
    for( size_t i = 0; i < ntot; ) {
        size_t j = i;
        while( j < ntot && _delta_dst[j] == _delta_dst[i] )
            j++;
        _delta_group.push_back( i );
        i = j;
    }
    _delta_group.push_back( ntot );
}


void GMG_Precond::prepare( const CRowMatrix &A )
{
    if( A.rows() != A.columns() )
        throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );

    int nfull = (int)( _nx0 * (_dim>=2?_ny0:1) * (_dim>=3?_nz0:1) );

    if( _node_map.empty() ) {
        // No elimination map supplied -- fall back to assuming every
        // mesh node has its own row (A.rows() == nx*ny*nz).
        if( A.rows() != nfull )
            throw( ErrorDim( ERROR_LOCATION, "matrix size does not match nx*ny*nz given to "
                              "GMG_Precond, and no node map was supplied via set_node_map(). "
                              "If your matrix assembler eliminates fixed/Dirichlet nodes "
                              "(as EpotMatrixSolver does, giving A.rows() == number of free "
                              "nodes rather than nx*ny*nz), build the mesh-node -> matrix-row "
                              "map from that assembler and pass it to set_node_map() before "
                              "calling construct()." ) );
        _node_map.assign( nfull, 0 );
        for( int m = 0; m < nfull; m++ )
            _node_map[m] = m;
    } else if( (int)_node_map.size() != nfull ) {
        throw( ErrorDim( ERROR_LOCATION, "node map size does not match nx*ny*nz given to GMG_Precond" ) );
    }

    _row_to_node.assign( A.rows(), -1 );
    for( int m = 0; m < nfull; m++ ) {
        int32_t r = _node_map[m];
        if( r < 0 )
            continue;
        if( r >= A.rows() )
            throw( ErrorDim( ERROR_LOCATION, "node map references a matrix row that does not exist" ) );
        if( _row_to_node[r] != -1 )
            throw( ErrorDim( ERROR_LOCATION, "node map assigns two mesh nodes to the same matrix row" ) );
        _row_to_node[r] = m;
    }
    for( int r = 0; r < A.rows(); r++ )
        if( _row_to_node[r] == -1 )
            throw( ErrorDim( ERROR_LOCATION, "node map does not cover every matrix row -- check that "
                              "it was built from the same matrix/geometry pair passed to construct()" ) );

    _level.clear();

    build_level0( A );
    color_graph( _level[0], &_node_map );
    find_diagonals( _level[0] );

    uint32_t maxlevels = (_nlevels_req == 0 ? 0xFFFFFFFF : _nlevels_req);

    while( _level.size() < maxlevels ) {

        const Level &fine = _level.back();
        uint32_t cnx, cny, cnz;
        if( !next_level_size( fine, cnx, cny, cnz ) )
            break;

        Level coarse;
        coarse.nx = cnx;
        coarse.ny = (_dim>=2? cny : 1);
        coarse.nz = (_dim>=3? cnz : 1);
        coarse.n  = (int)coarse.nx*coarse.ny*coarse.nz;

        const std::vector<int32_t> *fine_map = (_level.size() == 1 ? &_node_map : NULL);
        build_transfer_operators( fine, coarse, fine_map );
        galerkin_coarsen( fine, coarse );
        color_graph( coarse, NULL );
        find_diagonals( coarse );

        _level.push_back( std::move(coarse) );
    }

    build_delta_update_tables();

    // Scratch buffers, one pair per level, sized once and reused by
    // every call to solve().
    _x.assign( _level.size(), std::vector<double>() );
    _b.assign( _level.size(), std::vector<double>() );
    _r.assign( _level.size(), std::vector<double>() );
    for( size_t l = 0; l < _level.size(); l++ ) {
        _x[l].assign( _level[l].n, 0.0 );
        _b[l].assign( _level[l].n, 0.0 );
        _r[l].assign( _level[l].n, 0.0 );
    }

    // Report how much of level 0 the smoother actually sweeps. Level 0 is
    // full-mesh sized by construction (the geometric transfer operators
    // need a structured grid), so on a geometry with many eliminated nodes
    // -- Dirichlet electrodes, or a Neumann mask -- the gap between this
    // and the level size is the work the colour lists now skip.
    {
        size_t swept = _level[0].red.size() + _level[0].black.size()
                     + _level[0].leftover.size();
        ibsimu.message( 1 ) << "GMG_Precond: smoother sweeps " << swept
                            << " of " << _level[0].n << " level-0 rows ("
                            << (100.0*swept/_level[0].n) << " %), "
                            << (_level[0].n - swept)
                            << " eliminated rows skipped\n";
    }

    ibsimu.message( 1 ) << "GMG_Precond: built " << _level.size() << " level(s), sizes:";
    for( size_t l = 0; l < _level.size(); l++ )
        ibsimu.message( 1 ) << " " << _level[l].n;
    ibsimu.message( 1 ) << "\n";
    if( nfull > 0 )
        ibsimu.message( 1 ) << "GMG_Precond: " << (nfull-A.rows()) << " of " << nfull
                            << " mesh nodes (" << (100.0*(nfull-A.rows())/nfull)
                            << "%) are fixed/eliminated at level 0\n";
    ibsimu.flush();

    _prepared = true;
}


/* prepare() (called once, lazily, from here on first use) builds
 * everything that depends only on A's sparsity PATTERN: level sizing,
 * the geometric transfer operators P/R, and the red/black/leftover
 * coloring at every level -- all fixed for the lifetime of this
 * preconditioner, since the discretization's nonzero pattern never
 * changes (only the nonlinear plasma model's Jacobian VALUES do, via
 * its diagonal contribution). prepare() also performs an initial
 * values pass (level-0 embedding + Galerkin coarsening) as a side
 * effect of bootstrapping the hierarchy, which is why construct()
 * does not need to repeat that work on the very first call.
 *
 * Every call after the first must still refresh VALUE-dependent state
 * throughout the hierarchy: level 0's embedded values (from A, which
 * generally differs from the matrix prepare() last saw -- e.g. a new
 * Newton round or new major cycle), level 0's diagonal, and every
 * coarser level's Galerkin-coarsened operator + diagonal (rebuilt from
 * the now-current fine values). Skipping this after the first call was
 * a real bug: it silently reused a hierarchy frozen at whatever
 * Jacobian prepare() happened to see once, for the rest of the
 * preconditioner's lifetime, degrading (though never breaking
 * correctness of) BiCGSTAB's convergence rate more and more as the
 * actual matrix drifted from that frozen snapshot.
 */
void GMG_Precond::construct( const CRowMatrix &A )
{
    if( !_prepared ) {
        prepare( A );
        return;
    }

    refresh_level0_values( A );
    find_diagonals( _level[0] );

    if( _level.size() >= 2 ) {
        // Cheap incremental update of level 1's Galerkin operator --
        // see build_delta_update_tables() for the derivation. This
        // replaces what used to be a full general sparse-sparse
        // Galerkin product here (the single most expensive one in the
        // hierarchy, since level 0 is by far the largest level) with
        // an O(nnz) reset + scatter-accumulate.
        Level &lev1 = _level[1];
        lev1.val = _lvl1_lin_val;

        const Level &lev0 = _level[0];
        size_t ntriples = _delta_src.size();

        // One destination per iteration, summed in fixed index order -- no
        // atomics, and bit-reproducible. See the grouping in
        // build_delta_triples().
        (void)ntriples;
        const size_t ngroup = _delta_group.size() - 1;
        #pragma omp parallel for schedule(static)
        for( size_t g = 0; g < ngroup; g++ ) {
            double acc = 0.0;
            for( size_t k = _delta_group[g]; k < _delta_group[g+1]; k++ ) {
                int m = _delta_src[k];
                int di = lev0.diag_idx[m];
                double curdiag = (di >= 0 ? lev0.val[di] : 0.0);
                double delta = curdiag - _lvl0_diag_ref[m];
                if( delta == 0.0 )
                    continue;
                acc += delta*_delta_coef[k];
            }
            if( acc != 0.0 )
                lev1.val[ _delta_dst[ _delta_group[g] ] ] += acc;
        }

        find_diagonals( lev1 );
    }

    for( size_t l = 2; l < _level.size(); l++ ) {
        galerkin_coarsen( _level[l-1], _level[l] );
        find_diagonals( _level[l] );
    }
}


/* *****************************************************************************
 * Smoother (same red-black SOR as RBSOR_Precond, per level)
 */

void GMG_Precond::relax( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                          const std::vector<uint32_t> &rows, double w ) const
{
    const int *ptr = lev.ptr.data();
    const int *col = lev.col.data();
    const double *val = lev.val.data();
    const double *inv_diag = lev.inv_diag.data();
    const double w2 = 1.0-w;

    #pragma omp for schedule(static)
    for( size_t k = 0; k < rows.size(); k++ ) {

        uint32_t i = rows[k];
        double inv = inv_diag[i];
        if( inv == 0.0 )
            continue;

        double sum = b[i];
        int p0 = ptr[i], p1 = ptr[i+1];
        for( int p = p0; p < p1; p++ ) {
            int j = col[p];
            if( j != (int)i )
                sum -= val[p]*x[j];
        }
        x[i] = w*sum*inv + w2*x[i];
    }
}


/* Leftover (un-colourable) rows, updated serially in index order --
 * identical rationale to RBSOR_Precond::sweep_serial(). Called from
 * inside an "omp single" section by the caller.
 */
void GMG_Precond::relax_serial( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                                 const std::vector<uint32_t> &rows, double w ) const
{
    const double w2 = 1.0-w;

    for( size_t k = 0; k < rows.size(); k++ ) {
        uint32_t i = rows[k];
        int di = lev.diag_idx[i];
        if( di < 0 )
            continue;
        double aii = lev.val[di];
        if( aii == 0.0 )
            continue;

        double sum = b[i];
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            int j = lev.col[p];
            if( j == (int)i )
                continue;
            sum -= lev.val[p]*x[j];
        }
        double xnew = sum/aii;
        x[i] = w*xnew + w2*x[i];
    }
}


/* Squared L2 norm of the residual r = b - A*x on level lev, computed
 * in parallel (residual() already is), then reduced to a single
 * scalar visible to every thread via "single ... copyprivate" -- the
 * standard OpenMP idiom for "one thread computes a value, all threads
 * need to see it", used here rather than a reduction clause because
 * this is called from deep inside the recursive vcycle() call chain,
 * not directly inside the enclosing "omp parallel" block, so there is
 * no single shared reduction variable in scope to reduce into.
 */
double GMG_Precond::residual_norm2( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                                     std::vector<double> &r ) const
{
    residual( lev, x, b, r );

    double s = 0.0;
    #pragma omp single copyprivate(s)
    {
        double acc = 0.0;
        for( int i = 0; i < lev.n; i++ )
            acc += r[i]*r[i];
        s = acc;
    }
    return( s );
}


/* Coarsest-level solve: relax until the residual has dropped by a
 * factor of at least _coarse_rtol relative to where it started, or
 * _ncoarse_sweeps sweeps have been done, whichever comes first. This
 * replaces always doing a fixed number of sweeps -- on an easy right
 * hand side (e.g. later BiCGSTAB iterations, once the correction is
 * already small) a handful of sweeps is often enough, while on a hard
 * one _ncoarse_sweeps acts as a safety cap rather than a promise of
 * convergence. _last_coarse_sweeps records how many sweeps were
 * actually used on the most recent call, for tuning/diagnostics (see
 * last_coarse_sweeps()).
 */
void GMG_Precond::smooth_to_convergence( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                                          uint32_t max_sweeps, double rtol, double w ) const
{
    std::vector<double> &r = _r.back(); // scratch buffer sized for the coarsest level

    double rn2_0 = residual_norm2( lev, x, b, r );
    uint32_t s = 0;

    if( rn2_0 > 0.0 ) {
        for( ; s < max_sweeps; s++ ) {
            relax( lev, x, b, lev.red, w );
            relax( lev, x, b, lev.black, w );
            #pragma omp single
            {
                if( !lev.leftover.empty() )
                    relax_serial( lev, x, b, lev.leftover, w );
            }

            double rn2 = residual_norm2( lev, x, b, r );
            if( rn2 <= rtol*rtol*rn2_0 ) {
                s++;
                break;
            }
        }
    }

    #pragma omp single
    _last_coarse_sweeps = s;
}


void GMG_Precond::smooth( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                           uint32_t nsweeps, double w ) const
{
    for( uint32_t s = 0; s < nsweeps; s++ ) {
        relax( lev, x, b, lev.red, w );
        relax( lev, x, b, lev.black, w );
        #pragma omp single
        {
            if( !lev.leftover.empty() )
                relax_serial( lev, x, b, lev.leftover, w );
        }
    }
}


void GMG_Precond::residual( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                             std::vector<double> &r ) const
{
    const int *ptr = lev.ptr.data();
    const int *col = lev.col.data();
    const double *val = lev.val.data();

    #pragma omp for schedule(static)
    for( int i = 0; i < lev.n; i++ ) {
        double sum = b[i];
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum -= val[p]*x[col[p]];
        r[i] = sum;
    }
}


void GMG_Precond::restrict_vec( const Level &coarse, const std::vector<double> &rfine,
                                 std::vector<double> &rcoarse ) const
{
    const int *ptr = coarse.r_ptr.data();
    const int *col = coarse.r_col.data();
    const double *val = coarse.r_val.data();

    #pragma omp for schedule(static)
    for( int i = 0; i < coarse.n; i++ ) {
        double sum = 0.0;
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum += val[p]*rfine[col[p]];
        rcoarse[i] = sum;
    }
}


void GMG_Precond::prolong_add( const Level &coarse, const std::vector<double> &ecoarse,
                                std::vector<double> &xfine ) const
{
    const int *ptr = coarse.p_ptr.data();
    const int *col = coarse.p_col.data();
    const double *val = coarse.p_val.data();
    int nfine = (int)coarse.p_ptr.size()-1;

    #pragma omp for schedule(static)
    for( int i = 0; i < nfine; i++ ) {
        double sum = 0.0;
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum += val[p]*ecoarse[col[p]];
        xfine[i] += sum;
    }
}


/* *****************************************************************************
 * V-cycle
 */

void GMG_Precond::vcycle( uint32_t level, std::vector<double> &x, const std::vector<double> &b ) const
{
    const Level &lev = _level[level];

    if( level+1 == _level.size() ) {
        // Coarsest level: relax until converged (relative to its own
        // starting residual) or the sweep cap is hit. Uses its own
        // relaxation factor (_coarse_w) rather than the smoothing w,
        // since here fast convergence is wanted rather than good
        // high-frequency damping -- see the constructor docs.
        smooth_to_convergence( lev, x, b, _ncoarse_sweeps, _coarse_rtol, _coarse_w );
        return;
    }

    // Pre-smoothing
    smooth( lev, x, b, _npre, _w );

    // Residual on this level
    std::vector<double> &r = _r[level];
    residual( lev, x, b, r );

    // Restrict residual to become the coarse right-hand side
    const Level &clev = _level[level+1];
    std::vector<double> &bc = _b[level+1];
    restrict_vec( clev, r, bc );

    // Zero initial guess for the coarse correction
    std::vector<double> &xc = _x[level+1];
    #pragma omp single
    std::fill( xc.begin(), xc.end(), 0.0 );

    // Recurse to solve (approximately) for the correction
    vcycle( level+1, xc, bc );

    // Prolong correction back up and apply it
    prolong_add( clev, xc, x );

    // Post-smoothing
    smooth( lev, x, b, _npost, _w );
}


void GMG_Precond::solve( Vector &x, const Vector &b ) const
{
    if( !_prepared )
        throw( Error( ERROR_LOCATION, "GMG_Precond::solve() called before construct()" ) );

    int dof = (int)_row_to_node.size();
    if( (int)b.size() != dof )
        throw( ErrorDim( ERROR_LOCATION, "right hand side size does not match the matrix passed to construct()" ) );
    if( x.size() != b.size() )
        x.resize( b.size() );

    // Embed the dof-sized b into the full-mesh-sized level-0 working
    // vector: eliminated mesh slots get rhs=0 (their identity row
    // then keeps them at x=0 throughout, i.e. they are inert). Fixed
    // linear operator, always starting from zero -- same rationale as
    // RBSOR_Precond::solve().
    //
    // These four passes (two full-mesh fills, then the dof-sized
    // scatter in and gather out at the end) used to run single-threaded
    // while everything inside the V-cycle itself was parallel. That is
    // not a negligible tail: solve() is called TWICE per BiCGSTAB
    // iteration, so on a large mesh these serial memory-bound sweeps
    // showed up as a fixed serial cost on every preconditioner
    // application. They are all trivially data-parallel -- the fills
    // write disjoint elements, and _row_to_node is a permutation-like
    // injective map (each dof owns exactly one distinct mesh node), so
    // the scatter/gather have no write conflicts either.
    std::vector<double> &x0 = _x[0];
    std::vector<double> &b0 = _b[0];

    uint32_t nthreads = (_nthreads ? _nthreads : (uint32_t)omp_get_max_threads());

    const size_t n0 = x0.size();
    double *x0p = x0.data();
    double *b0p = b0.data();
    const int32_t *r2n = _row_to_node.data();
    const double *bp = b.get_data();

    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for schedule(static) nowait
        for( size_t n = 0; n < n0; n++ )
            x0p[n] = 0.0;
        #pragma omp for schedule(static)
        for( size_t n = 0; n < b0.size(); n++ )
            b0p[n] = 0.0;

        #pragma omp for schedule(static)
        for( int d = 0; d < dof; d++ )
            b0p[ r2n[d] ] = bp[d];

        for( uint32_t c = 0; c < _ncycles; c++ )
            vcycle( 0, x0, b0 );
    }

    // Extract the dof-sized solution back out of the embedded vector.
    double *xp = x.get_data();
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for( int d = 0; d < dof; d++ )
        xp[d] = x0p[ r2n[d] ];
}