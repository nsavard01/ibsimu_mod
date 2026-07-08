/*! \file rbsor_precond.hpp
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

#ifndef RBSOR_PRECOND_HPP
#define RBSOR_PRECOND_HPP 1


#include <vector>
#include <stdint.h>
#include "crowmatrix.hpp"
#include "precond.hpp"


/*! \brief Parallel red-black SOR/Gauss-Seidel preconditioner.
 *
 *  Implements the preconditioner \a M as a fixed number of red-black
 *  ordered SOR sweeps of the matrix \a A, i.e. for each sweep
 *
 *    x_i <- (1-w)*x_i + w/A_ii * ( b_i - sum_{j!=i} A_ij*x_j )
 *
 *  applied first to all "red" rows and then to all "black" rows, using
 *  in each half-sweep only values that are guaranteed not to be
 *  written by another row of the same half-sweep. This is exactly the
 *  IBSimu multigrid solver's structured-mesh red-black Gauss-Seidel
 *  smoother (see EpotMGSubSolver::rbgs_loop_3d() and friends in
 *  epot_mgsubsolver.cpp), generalized to work directly on the sparse
 *  CRowMatrix that EpotMatrixSolver assembles, rather than on the
 *  structured mesh.
 *
 *  Operating on the matrix graph (instead of reusing the i,j,k mesh
 *  coloring that the multigrid solver uses) has two advantages:
 *   - It works for the matrix actually handed to a Precond by
 *     EpotMatrixSolver-derived classes (e.g. EpotBiCGSTABSolver),
 *     whatever node numbering scheme was used to build it, without
 *     assuming the matrix row order follows Geometry's i,j,k indexing.
 *   - It is solver-agnostic: the same class can be used to
 *     precondition any sufficiently diagonally-dominant CRowMatrix
 *     coming from a different discretization, not just the specific
 *     stencils implemented in EpotMGSubSolver.
 *
 *  The coloring is a greedy graph 2-coloring of A's symmetrized
 *  sparsity pattern, built once in prepare() (called only when the
 *  non-zero pattern changes, e.g. once per geometry) and reused by
 *  every later construct()/solve() call, so its O(nnz) cost is paid
 *  only rarely, not on every nonlinear (Newton) or BiCGSTAB iteration.
 *
 *  NOTE: a greedy 2-coloring is not guaranteed to succeed (general
 *  graphs may need 3+ colors). When it fails for some rows, those rows
 *  are placed in an extra "leftover" group that is processed serially
 *  (still correct, just not parallelized). For the 5/7/9-point stencils
 *  that EpotMatrixSolver produces on a structured mesh, the standard
 *  even/odd checkerboard coloring always succeeds, so in practice the
 *  leftover group is empty -- this fallback exists purely so the class
 *  stays correct if it is ever pointed at a different matrix.
 */
class RBSOR_Precond : public Precond {

    uint32_t _nsweeps;     //!< \brief Number of red-black sweep pairs to do in solve()
    double   _w;           //!< \brief SOR relaxation factor (1.0 = plain Gauss-Seidel)
    uint32_t _nthreads;    //!< \brief Number of threads to use, 0 = ask ibsimu

    const CRowMatrix *_A;  //!< \brief Matrix most recently passed to construct()

    bool _prepared;        //!< \brief True after a successful prepare()

    // Greedy 2-coloring of A's symmetrized sparsity graph, computed in
    // prepare(). _color[i] is 0 (red), 1 (black) or 2 (leftover/serial)
    // for row i.
    std::vector<uint8_t>  _color;
    std::vector<uint32_t> _red;
    std::vector<uint32_t> _black;
    std::vector<uint32_t> _leftover;

    // Cached diagonal index (offset into A's col()/val() arrays) for
    // each row, so solve() does not have to scan each row for its
    // diagonal element every sweep. Built in prepare(); refreshed if
    // construct() is ever called with a matrix whose pattern changed
    // without a matching prepare() (defensive; should not happen in
    // normal use).
    std::vector<int> _diag_idx;
    std::vector<double> _inv_diag;

    RBSOR_Precond( const RBSOR_Precond &pc );
    const RBSOR_Precond &operator=( const RBSOR_Precond &pc );

    void color_graph( const CRowMatrix &A );
    void find_diagonals( const CRowMatrix &A );
    void sweep( Vector &x, const Vector &b, const std::vector<uint32_t> &rows ) const;
    void sweep_serial( Vector &x, const Vector &b, const std::vector<uint32_t> &rows ) const;
    

public:

    /*! \brief Constructor.
     *
     *  \param nsweeps  Number of red+black sweep pairs done in solve().
     *                  More sweeps give a more accurate preconditioner
     *                  solve (closer to the true SOR fixed point) at
     *                  the cost of more work per BiCGSTAB iteration.
     *                  1-2 is a reasonable starting point; this plays
     *                  the same role that the single ILU0 application
     *                  played in EpotBiCGSTABSolver's previous default
     *                  preconditioner.
     *  \param w        SOR relaxation factor. w = 1 gives plain
     *                  Gauss-Seidel. For the kind of close-to-Poisson
     *                  matrices EpotMatrixSolver produces, w slightly
     *                  above 1 (problem-size dependent, similar to
     *                  EpotMGSubSolver's SOR factor) usually improves
     *                  convergence; 1.0 is always a safe default.
     *  \param nthreads Number of threads to use for each colour sweep.
     *                  0 (default) asks ibsimu.get_thread_count().
     */
    RBSOR_Precond( uint32_t nsweeps = 4, double w = 1.4, uint32_t nthreads = 0 );

    ~RBSOR_Precond() {}

    RBSOR_Precond *copy( void ) const;

    void prepare( const CRowMatrix &A );

    void construct( const CRowMatrix &A );

    void clear( void );

    bool is_prepared( void ) const { return( _prepared ); }

    std::string typestring( void ) const;

    void solve( Vector &x, const Vector &b ) const;

    /*! \brief Set number of red-black sweep pairs done in solve().
     */
    void set_nsweeps( uint32_t nsweeps ) { _nsweeps = nsweeps; }

    /*! \brief Set SOR relaxation factor.
     */
    void set_relaxation( double w ) { _w = w; }
};


#endif
