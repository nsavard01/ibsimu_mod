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
#include <cmath>
#include <cstdint>
#include "gmg_precond.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


/* A loop over fewer rows than this runs on the calling thread, because
 * starting an OpenMP team costs more than such a loop gains. Measured
 * with 32 threads on 41k and 198k node meshes: the solve time is flat
 * for 256...2048 and rises from 8192 up.
 */
static const int GMG_PAR_MIN_ROWS = 1024;

/* The coarsest level is solved by dense LU up to this many rows, and by
 * relaxation beyond it. It is factorised on every construct(), which
 * costs 2/3*n^3 operations.
 */
static const int GMG_DIRECT_MAX_ROWS = 256;


/* Calls body(i) for i in [0,n). No iteration may depend on another.
 */
template <class Body>
static inline void for_rows( int n, int nthreads, const Body &body )
{
    if( nthreads > 1 && n >= GMG_PAR_MIN_ROWS ) {
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for( int i = 0; i < n; i++ )
            body( i );
    } else {
        for( int i = 0; i < n; i++ )
            body( i );
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


/* Returns an unprepared preconditioner with the same parameters, like
 * RBSOR_Precond::copy(): the hierarchy is rebuilt at first use.
 */
GMG_Precond *GMG_Precond::copy( void ) const
{
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
    _lvl0_delta.clear();
    _lvl1_lin_val.clear();
    _delta_ptr.clear();
    _delta_src.clear();
    _delta_coef.clear();
    _lu.clear();
    _lu_piv.clear();
}


std::string GMG_Precond::typestring( void ) const
{
    return( std::string( "GMG" ) );
}


int GMG_Precond::nthreads( void ) const
{
    return( _nthreads ? (int)_nthreads : (int)ibsimu.get_thread_count() );
}


/* *****************************************************************************
 * Hierarchy construction
 */

/* Level 0 is matrix A embedded in the full mesh: a free node gets its
 * row of A with the columns translated to mesh indices, an eliminated
 * node gets an identity row. Nothing refers to an eliminated node, so
 * those rows stay at x = 0 and are left out of the smoother, see
 * color_graph().
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


/* Copies the values of A into level 0, whose layout build_level0() made.
 */
void GMG_Precond::refresh_level0_values( const CRowMatrix &A )
{
    Level &lev = _level[0];

    for_rows( lev.n, nthreads(), [&]( int m ) {
        int32_t r = _node_map[m];
        if( r < 0 )
            return; // identity row, set by build_level0()
        int p = lev.ptr[m];
        for( int a = A.ptr(r); a < A.ptr(r+1); a++, p++ )
            lev.val[p] = A.val(a);
    } );
}


/* Size of the next coarser grid, (size+1)/2 in every active direction
 * as in EpotMGSolver::prepare_mg_geom(). Returns false where the
 * hierarchy has to stop: an even size cannot be halved with coincident
 * nodes, and a direction with fewer than 3 nodes has no interior left.
 */
bool GMG_Precond::next_level_size( const Level &fine, uint32_t &cnx, uint32_t &cny, uint32_t &cnz ) const
{
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

    if( cnx < 3 || (_dim >= 2 && cny < 3) || (_dim >= 3 && cnz < 3) )
        return( false );

    return( true );
}


/* Prolongation P (multilinear interpolation: weight 1 at the coincident
 * fine node, 1/2 per direction of offset) and restriction
 * R = (1/2^dim)*P^T from the coarse level to its parent \a fine.
 * Neighbours outside the grid are left out without renormalisation.
 *
 * fine_node_map is _node_map for the transfer from level 0, NULL
 * otherwise. Eliminated fine nodes are left out of P and R: they must
 * receive no correction, and with an empty P row the Galerkin product
 * couples nothing through them.
 */
void GMG_Precond::build_transfer_operators( const Level &fine, Level &coarse,
                                             const std::vector<int32_t> *fine_node_map ) const
{
    int fnx = fine.nx,   fny = (_dim>=2? fine.ny  :1),   fnz = (_dim>=3? fine.nz  :1);
    int cnx = coarse.nx, cny = (_dim>=2? coarse.ny :1),  cnz = (_dim>=3? coarse.nz :1);

    double c = 1.0;
    for( uint32_t d = 0; d < _dim; d++ )
        c *= 0.5;

    // (fine index, weight) pairs of every coarse node
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

                            int frow = fk*fnx*fny + fj*fnx + fi;
                            if( fine_node_map && (*fine_node_map)[frow] < 0 )
                                continue;
                            entries.push_back( std::make_pair(frow, wi*wj*wk) );
                        }
                    }
                }
            }
        }
    }

    // R = c * (the above), coarse.n rows x fine.n columns
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

    // P = the transpose, fine.n rows x coarse.n columns
    std::vector<int> col_count( fine.n, 0 );
    for( int i = 0; i < coarse.n; i++ )
        for( size_t k = 0; k < rows[i].size(); k++ )
            col_count[ rows[i][k].first ]++;

    coarse.p_ptr.assign( fine.n+1, 0 );
    for( int i = 0; i < fine.n; i++ )
        coarse.p_ptr[i+1] = coarse.p_ptr[i] + col_count[i];
    coarse.p_col.resize( coarse.p_ptr.back() );
    coarse.p_val.resize( coarse.p_ptr.back() );

    std::vector<int> cursor( coarse.p_ptr.begin(), coarse.p_ptr.end()-1 );
    for( int I = 0; I < coarse.n; I++ ) {
        for( size_t k = 0; k < rows[I].size(); k++ ) {
            int p = cursor[ rows[I][k].first ]++;
            coarse.p_col[p] = I;
            coarse.p_val[p] = rows[I][k].second;
        }
    }
}


/* Sparsity pattern of C = A*B in CSR form, columns ascending within a
 * row. It depends on the patterns of A and B only, never on their
 * values, so prepare() computes it once.
 */
static void spgemm_pattern( int m, int ncol,
                            const std::vector<int> &ptrA, const std::vector<int> &colA,
                            const std::vector<int> &ptrB, const std::vector<int> &colB,
                            std::vector<int> &ptrC, std::vector<int> &colC )
{
    std::vector<int> mark( ncol, -1 );
    ptrC.assign( m+1, 0 );
    colC.clear();

    for( int i = 0; i < m; i++ ) {
        for( int pa = ptrA[i]; pa < ptrA[i+1]; pa++ ) {
            int k = colA[pa];
            for( int pb = ptrB[k]; pb < ptrB[k+1]; pb++ ) {
                int c = colB[pb];
                if( mark[c] != i ) {
                    mark[c] = i;
                    colC.push_back( c );
                }
            }
        }
        std::sort( colC.begin()+ptrC[i], colC.end() );
        ptrC[i+1] = (int)colC.size();
    }
}


/* Values of C = A*B on the pattern made by spgemm_pattern(). Each row
 * is summed by one thread in a fixed order, so the result does not
 * depend on threading.
 */
static void spgemm_values( int m, int ncol, int nthreads,
                           const std::vector<int> &ptrA, const std::vector<int> &colA, const std::vector<double> &valA,
                           const std::vector<int> &ptrB, const std::vector<int> &colB, const std::vector<double> &valB,
                           const std::vector<int> &ptrC, const std::vector<int> &colC, std::vector<double> &valC )
{
    if( m < GMG_PAR_MIN_ROWS )
        nthreads = 1;

    #pragma omp parallel num_threads(nthreads)
    {
        std::vector<double> acc( ncol, 0.0 );

        #pragma omp for schedule(static)
        for( int i = 0; i < m; i++ ) {
            for( int pa = ptrA[i]; pa < ptrA[i+1]; pa++ ) {
                int k = colA[pa];
                double va = valA[pa];
                for( int pb = ptrB[k]; pb < ptrB[k+1]; pb++ )
                    acc[ colB[pb] ] += va*valB[pb];
            }
            for( int p = ptrC[i]; p < ptrC[i+1]; p++ ) {
                valC[p] = acc[ colC[p] ];
                acc[ colC[p] ] = 0.0;
            }
        }
    }
}


/* Pattern of the Galerkin operator coarse = R*(fine*P) and of the
 * intermediate product fine*P.
 */
void GMG_Precond::galerkin_pattern( const Level &fine, Level &coarse ) const
{
    spgemm_pattern( fine.n, coarse.n, fine.ptr, fine.col,
                    coarse.p_ptr, coarse.p_col, coarse.ap_ptr, coarse.ap_col );
    spgemm_pattern( coarse.n, coarse.n, coarse.r_ptr, coarse.r_col,
                    coarse.ap_ptr, coarse.ap_col, coarse.ptr, coarse.col );
    coarse.ap_val.assign( coarse.ap_col.size(), 0.0 );
    coarse.val.assign( coarse.col.size(), 0.0 );
}


/* Values of the Galerkin operator on the pattern made by
 * galerkin_pattern().
 */
void GMG_Precond::galerkin_values( const Level &fine, Level &coarse ) const
{
    const int nt = nthreads();
    spgemm_values( fine.n, coarse.n, nt, fine.ptr, fine.col, fine.val,
                   coarse.p_ptr, coarse.p_col, coarse.p_val,
                   coarse.ap_ptr, coarse.ap_col, coarse.ap_val );
    spgemm_values( coarse.n, coarse.n, nt, coarse.r_ptr, coarse.r_col, coarse.r_val,
                   coarse.ap_ptr, coarse.ap_col, coarse.ap_val,
                   coarse.ptr, coarse.col, coarse.val );
}


/* Greedy colouring of the rows: every row gets the lowest colour that
 * no row coupled to it has, where two rows are coupled if either one
 * has an entry in the other's column. Rows of one colour can therefore
 * be relaxed in any order, or in parallel, with the same result. A
 * 5-point stencil gives the usual red-black pair, a 9-point stencil
 * (every Galerkin level in 2D) four colours.
 *
 * node_map is _node_map for level 0, NULL otherwise. Eliminated rows
 * get no colour, so the smoother never visits them.
 */
void GMG_Precond::color_graph( Level &lev, const std::vector<int32_t> *node_map ) const
{
    const int n = lev.n;
    const int max_colors = 64;
    std::vector<int8_t> color( n, -1 );
    std::vector<uint64_t> taken( n, 0 ); // colours of the earlier rows that have an entry in column i

    lev.color.clear();
    for( int i = 0; i < n; i++ ) {
        if( node_map && (*node_map)[i] < 0 )
            continue;

        uint64_t t = taken[i];
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            int j = lev.col[p];
            if( j != i && color[j] >= 0 )
                t |= uint64_t(1) << color[j];
        }

        int c = 0;
        while( c < max_colors && (t & (uint64_t(1) << c)) )
            c++;
        if( c == max_colors )
            throw( Error( ERROR_LOCATION, "matrix row coupled to too many others to colour" ) );

        color[i] = c;
        if( c == (int)lev.color.size() )
            lev.color.push_back( std::vector<uint32_t>() );
        lev.color[c].push_back( i );

        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            int j = lev.col[p];
            if( j != i )
                taken[j] |= uint64_t(1) << c;
        }
    }
}


void GMG_Precond::find_diagonals( Level &lev ) const
{
    lev.diag_idx.assign( lev.n, -1 );
    lev.inv_diag.assign( lev.n, 0.0 );
    for( int i = 0; i < lev.n; i++ ) {
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ ) {
            if( lev.col[p] == i ) {
                lev.diag_idx[i] = p;
                break;
            }
        }
    }
}


/* inv_diag = 1/aii, or 0 for a row with no usable diagonal, which the
 * smoother then skips. A coarse node whose fine nodes are all
 * eliminated has such a row. The relaxation factor is applied by the
 * smoother, so that the same inv_diag serves any w.
 */
void GMG_Precond::update_inv_diag( Level &lev ) const
{
    for_rows( lev.n, nthreads(), [&]( int i ) {
        int di = lev.diag_idx[i];
        lev.inv_diag[i] = (di >= 0 && lev.val[di] != 0.0 ? 1.0/lev.val[di] : 0.0);
    } );
}


/* Tables for the cheap update of level 1 in construct(). Only the
 * diagonal of level 0 changes after prepare(), so with delta = (level 0
 * diagonal now) - (level 0 diagonal at prepare())
 *
 *   A1 = R*(A0_ref + diag(delta))*P = A1_ref + c * P^T*diag(delta)*P ,
 *
 * c = 1/2^dim. Entry (I,J) of the last term is
 * c * sum_m delta[m]*P[m,I]*P[m,J], so each fine row m adds
 * c*P[m,I]*P[m,J]*delta[m] to every level 1 entry (I,J) it
 * interpolates from. The (m, coefficient) pairs are stored here per
 * level 1 entry, in CSR form and in ascending m.
 */
void GMG_Precond::build_delta_tables( void )
{
    _lvl0_diag_ref.clear();
    _lvl0_delta.clear();
    _lvl1_lin_val.clear();
    _delta_ptr.clear();
    _delta_src.clear();
    _delta_coef.clear();

    if( _level.size() < 2 )
        return;

    const Level &lev0 = _level[0];
    const Level &lev1 = _level[1];

    _lvl0_diag_ref.assign( lev0.n, 0.0 );
    for( int m = 0; m < lev0.n; m++ ) {
        int di = lev0.diag_idx[m];
        if( di >= 0 )
            _lvl0_diag_ref[m] = lev0.val[di];
    }
    _lvl0_delta.assign( lev0.n, 0.0 );
    _lvl1_lin_val = lev1.val;

    double c = 1.0;
    for( uint32_t d = 0; d < _dim; d++ )
        c *= 0.5;

    // (level 1 entry, fine row, coefficient), in ascending fine row.
    std::vector<int32_t> dst, src;
    std::vector<double> coef;
    for( int m = 0; m < lev0.n; m++ ) {
        int p0 = lev1.p_ptr[m], p1 = lev1.p_ptr[m+1];
        for( int a = p0; a < p1; a++ ) {
            int Ia = lev1.p_col[a];
            for( int b = p0; b < p1; b++ ) {
                int Ib = lev1.p_col[b];
                for( int p = lev1.ptr[Ia]; p < lev1.ptr[Ia+1]; p++ ) {
                    if( lev1.col[p] == Ib ) {
                        dst.push_back( p );
                        src.push_back( m );
                        coef.push_back( c*lev1.p_val[a]*lev1.p_val[b] );
                        break;
                    }
                }
            }
        }
    }

    // Counting sort by level 1 entry; keeps the fine rows ascending.
    _delta_ptr.assign( lev1.val.size()+1, 0 );
    for( size_t k = 0; k < dst.size(); k++ )
        _delta_ptr[ dst[k]+1 ]++;
    for( size_t p = 0; p < lev1.val.size(); p++ )
        _delta_ptr[p+1] += _delta_ptr[p];

    _delta_src.resize( dst.size() );
    _delta_coef.resize( dst.size() );
    std::vector<size_t> cursor( _delta_ptr.begin(), _delta_ptr.end()-1 );
    for( size_t k = 0; k < dst.size(); k++ ) {
        size_t q = cursor[ dst[k] ]++;
        _delta_src[q]  = src[k];
        _delta_coef[q] = coef[k];
    }
}


/* Dense LU factorisation with partial pivoting of the coarsest
 * operator, for coarse_solve(). Leaves _lu empty if the level is too
 * large or the operator is singular; coarse_solve() then relaxes.
 */
void GMG_Precond::factor_coarsest( void )
{
    const Level &lev = _level.back();
    const int n = lev.n;

    _lu.clear();
    if( n > GMG_DIRECT_MAX_ROWS )
        return;

    _lu.assign( (size_t)n*n, 0.0 );
    _lu_piv.assign( n, 0 );
    for( int i = 0; i < n; i++ ) {
        if( lev.inv_diag[i] == 0.0 ) {
            // A row the smoother skips, see update_inv_diag().
            // coarse_solve() keeps its unknown at zero.
            _lu[i*n+i] = 1.0;
            continue;
        }
        for( int p = lev.ptr[i]; p < lev.ptr[i+1]; p++ )
            _lu[ i*n+lev.col[p] ] = lev.val[p];
    }

    for( int k = 0; k < n; k++ ) {
        int piv = k;
        for( int i = k+1; i < n; i++ ) {
            if( fabs( _lu[i*n+k] ) > fabs( _lu[piv*n+k] ) )
                piv = i;
        }
        if( _lu[piv*n+k] == 0.0 ) {
            _lu.clear();
            return;
        }
        _lu_piv[k] = piv;
        if( piv != k )
            std::swap_ranges( &_lu[k*n], &_lu[k*n]+n, &_lu[piv*n] );

        const double inv = 1.0/_lu[k*n+k];
        for( int i = k+1; i < n; i++ ) {
            const double f = _lu[i*n+k]*inv;
            _lu[i*n+k] = f;
            if( f != 0.0 ) {
                for( int j = k+1; j < n; j++ )
                    _lu[i*n+j] -= f*_lu[k*n+j];
            }
        }
    }
}


void GMG_Precond::prepare( const CRowMatrix &A )
{
    if( A.rows() != A.columns() )
        throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );

    int nfull = (int)( _nx0 * (_dim>=2?_ny0:1) * (_dim>=3?_nz0:1) );

    if( _node_map.empty() ) {
        // No map given: every mesh node has its own row.
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
    update_inv_diag( _level[0] );

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
        galerkin_pattern( fine, coarse );
        galerkin_values( fine, coarse );
        color_graph( coarse, NULL );
        find_diagonals( coarse );
        update_inv_diag( coarse );

        _level.push_back( std::move(coarse) );
    }

    build_delta_tables();
    if( _level.size() >= 2 ) {
        // Level 1 is updated through the delta tables from here on, so
        // its intermediate product is not needed again.
        Level &lev1 = _level[1];
        std::vector<int>().swap( lev1.ap_ptr );
        std::vector<int>().swap( lev1.ap_col );
        std::vector<double>().swap( lev1.ap_val );
    }
    factor_coarsest();

    _x.assign( _level.size(), std::vector<double>() );
    _b.assign( _level.size(), std::vector<double>() );
    _r.assign( _level.size(), std::vector<double>() );
    for( size_t l = 0; l < _level.size(); l++ ) {
        _x[l].assign( _level[l].n, 0.0 );
        _b[l].assign( _level[l].n, 0.0 );
        _r[l].assign( _level[l].n, 0.0 );
    }

    ibsimu.message( 1 ) << "GMG_Precond: built " << _level.size() << " level(s), sizes:";
    for( size_t l = 0; l < _level.size(); l++ )
        ibsimu.message( 1 ) << " " << _level[l].n;
    ibsimu.message( 1 ) << "\n";
    ibsimu.message( 1 ) << "GMG_Precond: coarsest level solved by "
                        << (_lu.empty() ? "relaxation" : "LU") << "\n";
    if( nfull > 0 )
        ibsimu.message( 1 ) << "GMG_Precond: " << (nfull-A.rows()) << " of " << nfull
                            << " mesh nodes (" << (100.0*(nfull-A.rows())/nfull)
                            << "%) are fixed/eliminated at level 0\n";
    ibsimu.flush();

    _prepared = true;
}


/* Values-only refresh of the hierarchy that prepare() built. Only the
 * diagonal of A may differ from the matrix prepare() saw (the plasma
 * and compensation terms of the Jacobian); any other change needs
 * clear() first.
 */
void GMG_Precond::construct( const CRowMatrix &A )
{
    if( !_prepared ) {
        prepare( A );
        return;
    }

    const int nt = nthreads();

    refresh_level0_values( A );
    update_inv_diag( _level[0] );

    if( _level.size() >= 2 ) {
        // Level 1 = its prepare() values + R*diag(delta)*P, see
        // build_delta_tables(). Each entry is summed by one thread in
        // a fixed order, so the result does not depend on threading.
        const Level &lev0 = _level[0];
        Level &lev1 = _level[1];

        for_rows( lev0.n, nt, [&]( int m ) {
            int di = lev0.diag_idx[m];
            _lvl0_delta[m] = (di >= 0 ? lev0.val[di] - _lvl0_diag_ref[m] : 0.0);
        } );

        for_rows( (int)lev1.val.size(), nt, [&]( int p ) {
            double acc = 0.0;
            for( size_t k = _delta_ptr[p]; k < _delta_ptr[p+1]; k++ )
                acc += _delta_coef[k]*_lvl0_delta[ _delta_src[k] ];
            lev1.val[p] = _lvl1_lin_val[p] + acc;
        } );

        update_inv_diag( lev1 );
    }

    for( size_t l = 2; l < _level.size(); l++ ) {
        galerkin_values( _level[l-1], _level[l] );
        update_inv_diag( _level[l] );
    }

    factor_coarsest();
}


/* *****************************************************************************
 * V-cycle
 */

/* SOR update of the rows of one colour. Rows of a colour do not read
 * each other (see color_graph()), so the result does not depend on
 * threading.
 */
void GMG_Precond::relax( const Level &lev, const std::vector<uint32_t> &rows,
                          std::vector<double> &x, const std::vector<double> &b, double w ) const
{
    const int *ptr = lev.ptr.data();
    const int *col = lev.col.data();
    const double *val = lev.val.data();
    const double *inv_diag = lev.inv_diag.data();
    const uint32_t *row = rows.data();
    const double *bp = b.data();
    double *xp = x.data();
    const double w2 = 1.0-w;

    for_rows( (int)rows.size(), nthreads(), [=]( int k ) {
        int i = row[k];
        double inv = inv_diag[i];
        if( inv == 0.0 )
            return;

        double sum = bp[i];
        for( int p = ptr[i]; p < ptr[i+1]; p++ ) {
            int j = col[p];
            if( j != i )
                sum -= val[p]*xp[j];
        }
        xp[i] = w*sum*inv + w2*xp[i];
    } );
}


void GMG_Precond::smooth( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                           uint32_t nsweeps, double w ) const
{
    for( uint32_t s = 0; s < nsweeps; s++ ) {
        for( size_t c = 0; c < lev.color.size(); c++ )
            relax( lev, lev.color[c], x, b, w );
    }
}


/* r = b - A*x
 */
void GMG_Precond::residual( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                             std::vector<double> &r ) const
{
    const int *ptr = lev.ptr.data();
    const int *col = lev.col.data();
    const double *val = lev.val.data();
    const double *xp = x.data();
    const double *bp = b.data();
    double *rp = r.data();

    for_rows( lev.n, nthreads(), [=]( int i ) {
        double sum = bp[i];
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum -= val[p]*xp[col[p]];
        rp[i] = sum;
    } );
}


/* rcoarse = R*rfine
 */
void GMG_Precond::restrict_vec( const Level &coarse, const std::vector<double> &rfine,
                                 std::vector<double> &rcoarse ) const
{
    const int *ptr = coarse.r_ptr.data();
    const int *col = coarse.r_col.data();
    const double *val = coarse.r_val.data();
    const double *fp = rfine.data();
    double *cp = rcoarse.data();

    for_rows( coarse.n, nthreads(), [=]( int i ) {
        double sum = 0.0;
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum += val[p]*fp[col[p]];
        cp[i] = sum;
    } );
}


/* xfine += P*ecoarse
 */
void GMG_Precond::prolong_add( const Level &coarse, const std::vector<double> &ecoarse,
                                std::vector<double> &xfine ) const
{
    const int *ptr = coarse.p_ptr.data();
    const int *col = coarse.p_col.data();
    const double *val = coarse.p_val.data();
    const double *cp = ecoarse.data();
    double *fp = xfine.data();

    for_rows( (int)xfine.size(), nthreads(), [=]( int i ) {
        double sum = 0.0;
        for( int p = ptr[i]; p < ptr[i+1]; p++ )
            sum += val[p]*cp[col[p]];
        fp[i] += sum;
    } );
}


/* Solves the coarsest level: exactly with the LU factors of
 * factor_coarsest() if there are any, otherwise by relaxing until the
 * residual has dropped by the factor _coarse_rtol or _ncoarse_sweeps
 * sweeps are done.
 */
void GMG_Precond::coarse_solve( std::vector<double> &x, const std::vector<double> &b ) const
{
    const Level &lev = _level.back();
    const int n = lev.n;
    _last_coarse_sweeps = 0;

    if( !_lu.empty() ) {
        for( int i = 0; i < n; i++ )
            x[i] = (lev.inv_diag[i] != 0.0 ? b[i] : 0.0);
        for( int k = 0; k < n; k++ )
            std::swap( x[k], x[ _lu_piv[k] ] );
        for( int i = 1; i < n; i++ ) {
            double sum = x[i];
            for( int j = 0; j < i; j++ )
                sum -= _lu[i*n+j]*x[j];
            x[i] = sum;
        }
        for( int i = n-1; i >= 0; i-- ) {
            double sum = x[i];
            for( int j = i+1; j < n; j++ )
                sum -= _lu[i*n+j]*x[j];
            x[i] = sum/_lu[i*n+i];
        }
        return;
    }

    std::vector<double> &r = _r.back();
    double rn2 = 0.0;
    residual( lev, x, b, r );
    for( int i = 0; i < n; i++ )
        rn2 += r[i]*r[i];
    const double target = _coarse_rtol*_coarse_rtol*rn2;

    while( rn2 > target && _last_coarse_sweeps < _ncoarse_sweeps ) {
        smooth( lev, x, b, 1, _coarse_w );
        _last_coarse_sweeps++;

        rn2 = 0.0;
        residual( lev, x, b, r );
        for( int i = 0; i < n; i++ )
            rn2 += r[i]*r[i];
    }
}


void GMG_Precond::vcycle( uint32_t level, std::vector<double> &x, const std::vector<double> &b ) const
{
    if( level+1 == _level.size() ) {
        coarse_solve( x, b );
        return;
    }

    const Level &lev = _level[level];
    const Level &clev = _level[level+1];
    std::vector<double> &r = _r[level];
    std::vector<double> &bc = _b[level+1];
    std::vector<double> &xc = _x[level+1];

    smooth( lev, x, b, _npre, _w );

    // Coarse correction: solve A_coarse*xc = R*(b - A*x) from xc = 0
    residual( lev, x, b, r );
    restrict_vec( clev, r, bc );
    std::fill( xc.begin(), xc.end(), 0.0 );
    vcycle( level+1, xc, bc );
    prolong_add( clev, xc, x );

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

    // The V-cycle works on the full mesh, starting from zero. An
    // eliminated node has rhs 0 and stays at 0.
    const int nt = nthreads();
    const int32_t *node_map = _node_map.data();
    const int32_t *row_to_node = _row_to_node.data();
    double *x0 = _x[0].data();
    double *b0 = _b[0].data();
    const double *bp = b.get_data();
    double *xp = x.get_data();

    for_rows( _level[0].n, nt, [=]( int m ) {
        int32_t r = node_map[m];
        b0[m] = (r >= 0 ? bp[r] : 0.0);
        x0[m] = 0.0;
    } );

    for( uint32_t c = 0; c < _ncycles; c++ )
        vcycle( 0, _x[0], _b[0] );

    for_rows( dof, nt, [=]( int d ) {
        xp[d] = x0[ row_to_node[d] ];
    } );
}
