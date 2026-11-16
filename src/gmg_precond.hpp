/*! \file gmg_precond.hpp
 *  \brief Parallel geometric multigrid preconditioner for sparse matrices
 *         arising from structured-grid finite-difference discretizations
 *         (companion to rbsor_precond.hpp).
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
 * (same license header as rbsor_precond.cpp -- reproduced there in full)
 */

#ifndef GMG_PRECOND_HPP
#define GMG_PRECOND_HPP

#include <cstdint>
#include <vector>
#include <string>

#include "precond.hpp"
#include "crowmatrix.hpp"


/*! \brief Parallel geometric multigrid (GMG) preconditioner.
 *
 * Implements the Precond interface for BiCGSTAB like RBSOR_Precond
 * does. solve() applies one or more V-cycles from x=0, which is a
 * fixed linear operator as BiCGSTAB requires.
 *
 * The caller gives the size of the finest structured grid
 * (nx,ny,nz,dim), with mesh node index m = k*nx*ny + j*nx + i. If the
 * matrix has no row for some mesh nodes (EpotMatrixSolver eliminates
 * the fixed nodes), set_node_map() must tell where each row sits on
 * the mesh before the first prepare()/construct(). The hierarchy is
 * built on the full mesh, with the eliminated nodes as inert identity
 * rows:
 *
 *   - Grids: every active direction is halved, size = (size+1)/2, for
 *     as long as the sizes are odd, as in
 *     EpotMGSolver::prepare_mg_geom(), or until \a nlevels is reached.
 *
 *   - Transfer: multilinear prolongation P and restriction
 *     R = (1/2^dim)*P^T, the stencils of EpotMGSolver::prolong_3d() /
 *     restrict_3d().
 *
 *   - Coarse operators: the Galerkin product A_{l+1} = R_l*A_l*P_l, so
 *     the near-solid, Neumann and dielectric stencils of the fine
 *     matrix carry over to every level without any knowledge of the
 *     geometry.
 *
 *   - Smoother: multi-colour SOR, the rows of one colour in parallel.
 *
 *   - Coarsest level: dense LU if it is small enough, otherwise SOR
 *     sweeps until converged.
 *
 * prepare() builds everything that depends on the sparsity pattern.
 * construct() then only refreshes values, and assumes that nothing
 * but the diagonal of the matrix has changed since prepare(); after
 * any other change call clear() (or set_node_map() with a new map)
 * first. The result of solve() does not depend on the number of
 * threads.
 */
class GMG_Precond : public Precond {
public:
    /*! \brief Constructor.
     *
     * \param nx,ny,nz        Finest grid size. Set ny=nz=1 for 1D,
     *                        nz=1 for 2D/cylindrical problems.
     * \param dim             Number of active grid directions (1,2,3).
     * \param nlevels         Maximum number of levels, including the
     *                        finest. 0 = as many as the grid allows.
     * \param npre,npost      Number of SOR sweeps before/after the
     *                        coarse correction on every level except
     *                        the coarsest.
     * \param ncoarse_sweeps  Maximum number of SOR sweeps on the
     *                        coarsest level. Used only if that level
     *                        is too large for the direct solve.
     * \param coarse_rtol     The coarsest-level sweeps stop once the
     *                        residual has dropped by this factor, in
     *                        (0,1).
     * \param w               SOR relaxation factor of the smoother, in
     *                        (0,2). 1.0 (Gauss-Seidel) damps the
     *                        high-frequency error best.
     * \param coarse_w        SOR relaxation factor of the
     *                        coarsest-level sweeps, in (0,2), where
     *                        over-relaxation helps convergence. 0 =
     *                        same as \a w.
     * \param ncycles         Number of V-cycles per solve().
     * \param nthreads        Number of threads. 0 =
     *                        ibsimu.get_thread_count().
     */
    GMG_Precond( uint32_t nx, uint32_t ny, uint32_t nz, uint32_t dim,
                 uint32_t nlevels = 0,
                 uint32_t npre = 2,
                 uint32_t npost = 2,
                 uint32_t ncoarse_sweeps = 50,
                 double coarse_rtol = 1e-2,
                 double w = 1.0,
                 double coarse_w = 0.0,
                 uint32_t ncycles = 1,
                 uint32_t nthreads = 0 );

    virtual ~GMG_Precond();

    virtual GMG_Precond *copy( void ) const;
    virtual void clear( void );
    virtual std::string typestring( void ) const;

    /*! \brief Build the hierarchy for matrix \a A: grids, transfer
     *  operators, coarse operators and colouring.
     */
    virtual void prepare( const CRowMatrix &A );

    /*! \brief Refresh the values of the hierarchy from matrix \a A,
     *  which may differ from the matrix given to prepare() only on
     *  its diagonal. Calls prepare() on first use.
     */
    virtual void construct( const CRowMatrix &A );

    /*! \brief Apply the V-cycles as an approximate solve of A*x=b,
     *  starting from x=0.
     */
    virtual void solve( Vector &x, const Vector &b ) const;

    /*! \brief Returns true once prepare() has built the hierarchy. */
    virtual bool is_prepared( void ) const { return( _prepared ); }

    /*! \brief Set the mesh node -> matrix row map.
     *
     *  \a node_row has nx*ny*nz entries in mesh order. Entry \a m is
     *  the row of mesh node \a m in the matrix, or -1 if the node has
     *  no row. Every row must appear exactly once; prepare() checks
     *  this. Without a map every mesh node is taken to have its own
     *  row.
     *
     *  A map equal to the current one changes nothing, so this can be
     *  called before every solve. A different map makes the next
     *  construct() rebuild the hierarchy.
     */
    void set_node_map( const std::vector<int32_t> &node_row ) {
	if( _prepared && node_row == _node_map )
	    return;
	_node_map = node_row;
	_prepared = false;
    }

    void set_ncycles( uint32_t ncyc ) { _ncycles = ncyc; }
    uint32_t levels( void ) const { return( (uint32_t)_level.size() ); }

    /*! \brief Number of SOR sweeps the coarsest level took in the
     *  most recent V-cycle. 0 if it is solved directly.
     */
    uint32_t last_coarse_sweeps( void ) const { return( _last_coarse_sweeps ); }

private:

    /* One level of the hierarchy. Level 0 is the finest. */
    struct Level {
        uint32_t nx, ny, nz;   // grid size
        int      n;            // number of rows = nx*ny*nz

        // Operator in CSR form: the matrix embedded in the full mesh
        // for level 0, the Galerkin operator R*A*P below it.
        std::vector<int>    ptr;
        std::vector<int>    col;
        std::vector<double> val;

        std::vector<int>    diag_idx;    // position of the diagonal in val, -1 if none
        std::vector<double> inv_diag;    // 1/diagonal, see update_inv_diag()

        // Rows of each colour, see color_graph().
        std::vector< std::vector<uint32_t> > color;

        // Prolongation P from this level to its parent (parent.n rows
        // x n columns) and restriction R = (1/2^dim)*P^T from the
        // parent to this level. Empty for level 0.
        std::vector<int>    p_ptr, p_col;
        std::vector<double> p_val;
        std::vector<int>    r_ptr, r_col;
        std::vector<double> r_val;

        // (parent operator)*P, the intermediate product of the
        // Galerkin operator, see galerkin_pattern().
        std::vector<int>    ap_ptr, ap_col;
        std::vector<double> ap_val;
    };

    uint32_t _nx0, _ny0, _nz0, _dim;
    uint32_t _nlevels_req;
    uint32_t _npre, _npost, _ncoarse_sweeps;
    double   _coarse_rtol;
    double   _w;
    double   _coarse_w;
    uint32_t _ncycles;
    uint32_t _nthreads;

    mutable uint32_t _last_coarse_sweeps;

    bool _prepared;
    std::vector<Level> _level;

    // Mesh node -> matrix row (-1 = no row) and its inverse, see
    // set_node_map().
    std::vector<int32_t> _node_map;
    std::vector<int32_t> _row_to_node;

    // Work vectors of solve(), one set per level.
    mutable std::vector< std::vector<double> > _x, _b, _r;

    // Update of level 1 from the change of the level 0 diagonal, see
    // build_delta_tables().
    std::vector<double>  _lvl0_diag_ref;  // level 0 diagonal at prepare()
    std::vector<double>  _lvl0_delta;     // level 0 diagonal now - _lvl0_diag_ref
    std::vector<double>  _lvl1_lin_val;   // level 1 values at prepare()
    std::vector<size_t>  _delta_ptr;      // per level 1 entry, its range in _delta_src/_delta_coef
    std::vector<int32_t> _delta_src;      // level 0 row
    std::vector<double>  _delta_coef;     // level 1 entry += _delta_coef * _lvl0_delta[_delta_src]

    // LU factors of the coarsest operator (n*n, row-major) and the
    // row swaps, see factor_coarsest(). Empty if it is not factorised.
    std::vector<double>  _lu;
    std::vector<int>     _lu_piv;

    int nthreads( void ) const;

    void build_level0( const CRowMatrix &A );
    void refresh_level0_values( const CRowMatrix &A );
    bool next_level_size( const Level &fine, uint32_t &cnx, uint32_t &cny, uint32_t &cnz ) const;
    void build_transfer_operators( const Level &fine, Level &coarse,
                                    const std::vector<int32_t> *fine_node_map ) const;
    void galerkin_pattern( const Level &fine, Level &coarse ) const;
    void galerkin_values( const Level &fine, Level &coarse ) const;
    void color_graph( Level &lev, const std::vector<int32_t> *node_map ) const;
    void find_diagonals( Level &lev ) const;
    void update_inv_diag( Level &lev ) const;
    void build_delta_tables( void );
    void factor_coarsest( void );

    void relax( const Level &lev, const std::vector<uint32_t> &rows,
                std::vector<double> &x, const std::vector<double> &b, double w ) const;
    void smooth( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                 uint32_t nsweeps, double w ) const;
    void residual( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                   std::vector<double> &r ) const;
    void restrict_vec( const Level &coarse, const std::vector<double> &rfine,
                        std::vector<double> &rcoarse ) const;
    void prolong_add( const Level &coarse, const std::vector<double> &ecoarse,
                       std::vector<double> &xfine ) const;
    void coarse_solve( std::vector<double> &x, const std::vector<double> &b ) const;
    void vcycle( uint32_t level, std::vector<double> &x, const std::vector<double> &b ) const;
};

#endif