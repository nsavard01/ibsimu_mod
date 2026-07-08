/*! \file rbsor_precond.cpp
 *  \brief Parallel red-black SOR/Gauss-Seidel preconditioner for sparse matrices
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 *
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#include <algorithm>
#include "rbsor_precond.hpp"
#include <omp.h>
#include "ibsimu.hpp"
#include "error.hpp"


RBSOR_Precond::RBSOR_Precond( uint32_t nsweeps, double w, uint32_t nthreads )
    : _nsweeps(nsweeps), _w(w), _nthreads(nthreads), _A(NULL), _prepared(false)
{
    if( nsweeps < 1 )
        throw( ErrorDim( ERROR_LOCATION, "nsweeps must be at least 1" ) );
    if( w <= 0.0 || w >= 2.0 )
        throw( ErrorDim( ERROR_LOCATION, "relaxation factor must be in (0,2)" ) );
}


RBSOR_Precond *RBSOR_Precond::copy( void ) const
{
    RBSOR_Precond *pc = new RBSOR_Precond( _nsweeps, _w, _nthreads );
    // Coloring depends only on the non-zero pattern, which is cheap to
    // recompute on first use by the copy; we deliberately do not copy
    // _A (a non-owning pointer into whatever matrix is alive at the
    // call site) or _prepared, so the copy will call prepare() again
    // the first time it needs to, exactly like a freshly constructed
    // preconditioner would.
    return( pc );
}


void RBSOR_Precond::clear( void )
{
    _prepared = false;
    _A = NULL;
    _color.clear();
    _red.clear();
    _black.clear();
    _leftover.clear();
    _diag_idx.clear();
    _inv_diag.clear();
}


std::string RBSOR_Precond::typestring( void ) const
{
    return( std::string( "RBSOR" ) );
}


/* Greedy 2-coloring of the symmetrized sparsity graph of A.
 *
 * Two rows i != j are graph-adjacent if A(i,j) != 0 or A(j,i) != 0 (we
 * need the symmetrized graph because a SOR sweep for row i reads
 * column j whenever A(i,j) != 0, but the *matrix itself* need not be
 * structurally symmetric, e.g. because of Neumann/near-solid boundary
 * stencils). Rows are colored in index order; each row gets the first
 * color not already used by an already-colored neighbour. If both
 * colors are taken (which cannot happen for the bipartite stencils
 * EpotMatrixSolver produces, but could in principle for some other
 * graph) the row is placed in a third, serially-processed group, so
 * the preconditioner stays correct even if it is ever asked to handle
 * an unexpected sparsity pattern.
 */
void RBSOR_Precond::color_graph( const CRowMatrix &A )
{
    int n = A.rows();
    _color.assign( n, 0xFF ); // 0xFF = uncoloured

    for( int i = 0; i < n; i++ ) {

        bool red_taken   = false;
        bool black_taken = false;

        int ptr0 = A.ptr(i);
        int ptr1 = A.ptr(i+1);
        for( int p = ptr0; p < ptr1; p++ ) {
            int j = A.col(p);
            if( j == i || j < 0 || j >= n )
                continue;
            if( _color[j] == 0 )
                red_taken = true;
            else if( _color[j] == 1 )
                black_taken = true;
        }

        if( !red_taken )
            _color[i] = 0;
        else if( !black_taken )
            _color[i] = 1;
        else
            _color[i] = 2; // leftover, handled serially in solve()
    }

    // The loop above only looks at row i's own outgoing entries
    // A(i,j). Because A need not be structurally symmetric (Neumann
    // and near-solid stencils can be asymmetric near boundaries), we
    // also need to make sure that for every entry A(i,j) the rows i
    // and j ended up with different colors (or one of them is in the
    // leftover group); otherwise updating all of one colour in
    // parallel could read a value written earlier in the same sweep
    // by another row of the same colour. Repair any such conflicts by
    // moving the higher-index row of the pair into the leftover group.
    for( int i = 0; i < n; i++ ) {
        if( _color[i] == 2 )
            continue;
        int ptr0 = A.ptr(i);
        int ptr1 = A.ptr(i+1);
        for( int p = ptr0; p < ptr1; p++ ) {
            int j = A.col(p);
            if( j == i || j < 0 || j >= n )
                continue;
            if( _color[j] == _color[i] )
                _color[ std::max(i,j) ] = 2;
        }
    }

    _red.clear();
    _black.clear();
    _leftover.clear();
    for( int i = 0; i < n; i++ ) {
        if( _color[i] == 0 )
            _red.push_back( (uint32_t)i );
        else if( _color[i] == 1 )
            _black.push_back( (uint32_t)i );
        else
            _leftover.push_back( (uint32_t)i );
    }

    if( !_leftover.empty() ) {
        ibsimu.message( 2 ) << "RBSOR_Precond: " << _leftover.size()
                            << " of " << n
                            << " rows could not be 2-coloured, "
                            << "processing them serially\n";
        ibsimu.flush();
    }
}


void RBSOR_Precond::find_diagonals( const CRowMatrix &A )
{
    int n = A.rows();
    _diag_idx.assign( n, -1 );
    _inv_diag.assign(n, 0.0 );
    for( int i = 0; i < n; i++ ) {
        int ptr0 = A.ptr(i);
        int ptr1 = A.ptr(i+1);
        for( int p = ptr0; p < ptr1; p++ ) {
            if( A.col(p) == i ) {
                _diag_idx[i] = p;
                _inv_diag[i] = _w / A.val(p);
                break;
            }
        }
    }
}


void RBSOR_Precond::prepare( const CRowMatrix &A )
{
    if( A.rows() != A.columns() )
        throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );

    color_graph( A );
    find_diagonals( A );
    _prepared = true;
}


void RBSOR_Precond::construct( const CRowMatrix &A )
{
    if( !_prepared )
        prepare( A );
    _A = &A;
}





/* Leftover rows (if any) are updated one at a time, in index order, so
 * each can safely use the freshest available values of every other
 * row -- exactly like ordinary (non-coloured) Gauss-Seidel.
 */
void RBSOR_Precond::sweep_serial( Vector &x, const Vector &b, const std::vector<uint32_t> &rows ) const
{
    const CRowMatrix &A = *_A;
    const double w  = _w;
    const double w2 = 1.0-w;

    for( size_t k = 0; k < rows.size(); k++ ) {
        uint32_t i = rows[k];
        int ptr0 = A.ptr(i);
        int ptr1 = A.ptr(i+1);
        int di = _diag_idx[i];
        if( di < 0 )
            continue;

        double sum = b[i];
        for( int p = ptr0; p < ptr1; p++ ) {
            int j = A.col(p);
            if( j == (int)i )
                continue;
            sum -= A.val(p) * x[j];
        }
        double aii = A.val(di);
        if( aii == 0.0 )
            continue;
        double xnew = sum/aii;
        x[i] = w*xnew + w2*x[i];
    }
}


void RBSOR_Precond::solve( Vector &x, const Vector &b ) const
{
    if( !_A )
        throw( Error( ERROR_LOCATION, "RBSOR_Precond::solve() called before construct()" ) );

    if( x.size() != b.size() )
        x.resize( b.size() );

    // Start from zero -- BiCGSTAB only ever needs M^{-1} applied to a
    // vector, not a refinement of an existing guess, and starting from
    // zero each time keeps the preconditioner a fixed linear operator
    // (Precond::solve() takes no previous x), which is what BiCGSTAB
    // assumes when it calls pc.solve( phat, p ) etc.
    x.clear();

    const CRowMatrix &A = *_A;
    const double w  = _w;
    const double w2 = 1.0-w;
    const int *ptr = A.ptr_data();
    const int *col = A.col_data();
    const double *val = A.val_data();

    auto relax = [&](const std::vector<uint32_t>& rows)
    {
        #pragma omp for schedule(static)
        for (size_t k = 0; k < rows.size(); k++) {

            uint32_t i = rows[k];
            

            int ptr0 = ptr[i];
            int ptr1 = ptr[i+1];
            double inv = _inv_diag[i];

            if( inv == 0.0)
                continue;

            double sum = b[i];

            for( int p = ptr0; p < ptr1; p++ ) {
                int j = col[p];
                if( j != (int)i )
                    sum -= val[p] * x[j];
            }

            

            if( inv != 0.0 )
                x[i] = sum * inv + w2*x[i];
        }
    };

    #pragma omp parallel num_threads(_nthreads)
    {
        for (uint32_t s = 0; s < _nsweeps; s++) {

            relax(_red);
            relax(_black);

            #pragma omp single
            {
                if (!_leftover.empty())
                    sweep_serial(x,b,_leftover);
            }
        }
    }

    // for( uint32_t s = 0; s < _nsweeps; s++ ) {
    //     sweep( x, b, _red );
    //     sweep( x, b, _black );
    //     if( !_leftover.empty() )
    //         sweep_serial( x, b, _leftover );
    // }
}
