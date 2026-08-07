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

/* Confirmed against the actual ibsimu headers: Precond declares
 * prepare()/construct()/solve()/clear()/copy()/is_prepared()/typestring()
 * as pure virtual, and CRowMatrix/Vector come in transitively through
 * precond.hpp (which includes crowmatrix.hpp and mvector.hpp). Including
 * crowmatrix.hpp explicitly here too is harmless (include-guarded) and
 * keeps this header self-sufficient if ever used without precond.hpp.
 */
#include "precond.hpp"
#include "crowmatrix.hpp"


/*! \brief Parallel geometric multigrid (GMG) preconditioner.
 *
 * Meant as a drop-in alternative to RBSOR_Precond for BiCGSTAB: it
 * implements the same Precond interface (construct()/solve()), and
 * internally performs one or more V-cycles as its approximate M^-1
 * application, always starting from x=0 (a fixed linear operator, as
 * required by BiCGSTAB's use of the preconditioner).
 *
 * The caller supplies the size of the *finest* structured grid
 * (nx,ny,nz,dim) matching the row ordering a = k*nx*ny + j*nx + i.
 *
 * IMPORTANT -- eliminated (fixed/Dirichlet) nodes: matrix assemblers
 * such as EpotMatrixSolver do not give every mesh node a row in A --
 * fixed nodes are eliminated and folded into the right-hand side, so
 * A.rows() is the number of *free* nodes ("dof"), which is normally
 * smaller than nx*ny*nz, and row index r in A has no direct relation
 * to mesh index a. Geometric coarsening needs to know where every
 * matrix row actually sits on the structured mesh, so if your matrix
 * eliminates nodes this way you must call set_node_map() once (see
 * below) before the first construct()/prepare() call. If
 * set_node_map() is never called, GMG_Precond falls back to assuming
 * a 1:1 identity mapping (A.rows() == nx*ny*nz, no elimination),
 * which is only correct for assemblers that keep every mesh node as
 * an explicit row.
 *
 * Internally, eliminated nodes are embedded into the full-mesh-sized
 * multigrid hierarchy as trivial decoupled identity rows (always
 * x=0): this lets the geometric transfer operators and Galerkin
 * coarsening operate on a proper structured grid at every level
 * without needing to know anything about which nodes were
 * eliminated beyond the map itself. solve() transparently embeds the
 * dof-sized b into the full-mesh-sized working vector and extracts
 * the dof-sized x back out, so the public interface is unaffected.
 * From that, the class builds:
 *
 *   - A hierarchy of coarser grids, halving each active dimension
 *     (size[d] = (size[d]+1)/2) as long as the parity allows, exactly
 *     like EpotMGSolver::prepare_mg_geom() does. Coarsening stops
 *     automatically once further halving is not exactly representable
 *     (an even dimension) or a user cap (nlevels) is hit.
 *
 *   - Geometric prolongation/restriction operators using the same
 *     multilinear stencil weights as EpotMGSolver::prolong_3d() /
 *     restrict_3d() (generalized to 1D/2D/3D), stored explicitly as
 *     sparse matrices so they can also be used to build the coarse
 *     operators.
 *
 *   - Coarse-grid operators via the Galerkin product
 *         A_{l+1} = R_l * A_l * P_l ,   R_l = (1/2^dim) * P_l^T
 *     computed once in construct(). This makes the method correct for
 *     whatever near-solid / Neumann-modified stencil is present in
 *     the fine matrix, without needing to know anything about
 *     Geometry or plasma boundary conditions on the coarser levels --
 *     it is inherited automatically through the Galerkin projection.
 *
 *   - A red-black(-leftover) coloring of every level's operator, so
 *     that the smoother (plain SOR / Gauss-Seidel, matching
 *     RBSOR_Precond) can run in parallel exactly like
 *     RBSOR_Precond::solve() does.
 *
 * Rows that are not connected to the structured grid in the expected
 * way (all boundary/fixed rows for instance) are handled the same as
 * in RBSOR_Precond: rows with zero diagonal are simply skipped by the
 * smoother (a no-op relaxation), which is safe because the Galerkin
 * coarsening naturally keeps such rows decoupled.
 */
class GMG_Precond : public Precond {
public:
    /*! \brief Constructor.
     *
     * \param nx,ny,nz        Finest grid size. Set ny=nz=1 for 1D,
     *                        nz=1 for 2D/cylindrical problems.
     * \param dim              Number of active grid directions (1,2,3).
     * \param nlevels          Maximum number of multigrid levels,
     *                        including the finest. 0 = automatic (coarsen
     *                        for as long as the mesh parity allows).
     * \param npre,npost       Number of red-black SOR sweeps used for
     *                        pre-/post-smoothing on every level except
     *                        the coarsest.
     * \param ncoarse_sweeps   Maximum number of SOR sweeps used as an
     *                        (approximate) direct solve on the coarsest
     *                        level -- a safety cap, not a target: the
     *                        coarse solve stops early once its residual
     *                        has dropped by a factor of coarse_rtol
     *                        relative to where it started (see
     *                        coarse_rtol below and last_coarse_sweeps()).
     * \param coarse_rtol      Relative residual reduction target for the
     *                        coarsest-level solve, in (0,1). E.g. 1e-2
     *                        stops once the coarse residual has dropped
     *                        two orders of magnitude. Since this is only
     *                        ever used as a preconditioner (an
     *                        approximate M^-1, not an exact solve is
     *                        required), a modest reduction is normally
     *                        enough; tightening it trades more coarse
     *                        sweeps for a slightly better V-cycle.
     * \param w                SOR relaxation factor used for pre-/post-
     *                        smoothing, in (0,2). 1.0 (plain
     *                        Gauss-Seidel) is a safe, robust default:
     *                        as a multigrid *smoother* what matters is
     *                        damping high-frequency error efficiently,
     *                        which w=1 already does close to optimally
     *                        -- over-relaxing here trades that away for
     *                        faster overall convergence, which isn't
     *                        the smoother's job.
     * \param coarse_w         SOR relaxation factor used only for the
     *                        coarsest-level solve, in (0,2). Unlike the
     *                        smoothing w above, the coarsest level is
     *                        meant to actually converge, so classical
     *                        over-relaxation (coarse_w > 1) is often
     *                        beneficial there: for a plain discrete
     *                        Laplacian the optimal factor is
     *                        approximately 2/(1+sin(pi/(N+1))), N being
     *                        the coarsest grid's linear size, which
     *                        tends toward the upper end of (1,2) for
     *                        small N -- e.g. ~1.7 for N=19. That
     *                        formula assumes an idealized full-domain
     *                        stencil though, so treat it as a starting
     *                        point and tune from there using
     *                        last_coarse_sweeps() as feedback (a value
     *                        that's too high can also slow convergence
     *                        or destabilize it, so search rather than
     *                        just cranking it up). Defaults to the same
     *                        value as w if left at 0.
     * \param ncycles          Number of V-cycles performed per call to
     *                        solve(). 1 is normally sufficient since
     *                        this is used as a preconditioner inside
     *                        BiCGSTAB rather than as a standalone
     *                        solver.
     * \param nthreads         OpenMP threads to use. 0 = let OpenMP
     *                        decide (omp_get_max_threads()).
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

    /*! \brief Build the grid hierarchy, transfer operators and Galerkin
     *  coarse operators from matrix A. Required by the Precond
     *  interface; construct() calls this automatically the first
     *  time it is needed, exactly like RBSOR_Precond does.
     */
    virtual void prepare( const CRowMatrix &A );

    /*! \brief Build the preconditioner for matrix A. Calls prepare()
     *  on first use (mirrors RBSOR_Precond::construct()).
     */
    virtual void construct( const CRowMatrix &A );

    /*! \brief Apply _ncycles V-cycles as an approximate solve of A*x=b,
     *  always starting from x=0.
     */
    virtual void solve( Vector &x, const Vector &b ) const;

    /*! \brief Returns true once prepare() has built the hierarchy. */
    virtual bool is_prepared( void ) const { return( _prepared ); }

    /*! \brief Supply the mesh-node -> matrix-row map for assemblers
     *  that eliminate fixed/Dirichlet nodes (see class documentation).
     *
     *  \a node_row must have exactly nx*ny*nz entries, in the same
     *  a = k*nx*ny + j*nx + i order as the mesh. Entry \a m gives the
     *  row index of mesh node \a m in the matrix A that will be
     *  passed to construct(), or -1 if that node is eliminated
     *  (fixed/Dirichlet, no row of its own). Every row of A must be
     *  referenced by exactly one entry -- prepare() checks this and
     *  throws if not.
     *
     *  Safe (and expected) to call this again with an unchanged map --
     *  e.g. a caller that doesn't itself track whether the fixed-node
     *  set changed since the last solve, and so calls this once before
     *  every solve just to be safe. If \a node_row is identical to
     *  whatever was already prepared, this is a no-op: it does *not*
     *  invalidate the hierarchy, so the next construct() takes the
     *  cheap incremental-update path instead of a full prepare()
     *  rebuild. Only a genuine change forces a rebuild. If the mesh or
     *  the set of fixed nodes *does* change, call this (together with
     *  clear(), if you want to force starting over rather than let the
     *  next prepare() reuse scratch buffers) before the next
     *  construct()/prepare() call. Not needed at all if the matrix
     *  assembler keeps one row per mesh node (A.rows() == nx*ny*nz).
     */
    void set_node_map( const std::vector<int32_t> &node_row ) {
	if( _prepared && node_row == _node_map )
	    return; // unchanged since the last prepare() -- nothing to invalidate
	_node_map = node_row;
	_prepared = false;
    }

    void set_ncycles( uint32_t ncyc ) { _ncycles = ncyc; }
    uint32_t levels( void ) const { return( (uint32_t)_level.size() ); }

    /*! \brief Number of sweeps actually used by the coarsest-level
     *  solve on the most recent solve() call (see coarse_rtol in the
     *  constructor). Useful for tuning ncoarse_sweeps/coarse_rtol:
     *  consistently hitting the ncoarse_sweeps cap means it's too low
     *  (or coarse_rtol too tight) for this problem; consistently
     *  finishing in a handful of sweeps means there's room to lower
     *  ncoarse_sweeps without any loss of quality.
     */
    uint32_t last_coarse_sweeps( void ) const { return( _last_coarse_sweeps ); }

private:

    /* One level of the multigrid hierarchy. Level 0 is the finest
     * (user-supplied) matrix, level size()-1 is the coarsest. */
    struct Level {
        uint32_t nx, ny, nz;   // grid size at this level
        int      n;            // number of unknowns = nx*ny*nz

        // CSR copy of the operator at this level (straight copy of A
        // for level 0, Galerkin operator R*A*P for level>0).
        std::vector<int>    ptr;
        std::vector<int>    col;
        std::vector<double> val;

        std::vector<int>    diag_idx;
        std::vector<double> inv_diag;    // w/aii, precomputed for the smoother

        // Greedy 2-coloring of the symmetrized sparsity graph (same
        // algorithm as RBSOR_Precond::color_graph()).
        std::vector<uint32_t> red, black, leftover;

        // Prolongation P: parent(fine).n rows x this(coarse).n cols,
        // i.e. maps a correction on *this* (coarser) level up to the
        // parent (finer) level. Empty for the finest level.
        std::vector<int>    p_ptr, p_col;
        std::vector<double> p_val;

        // Restriction R = (1/2^dim) * P^T: this(coarse).n rows x
        // parent(fine).n cols, i.e. maps a fine-level residual down to
        // this (coarser) level. Empty for the finest level.
        std::vector<int>    r_ptr, r_col;
        std::vector<double> r_val;
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

    // Mesh-node (m = k*nx*ny+j*nx+i, full nx0*ny0*nz0 grid) -> matrix
    // row map and its inverse, see set_node_map(). _node_map has one
    // entry per mesh node (-1 = eliminated/fixed); _row_to_node has
    // one entry per matrix row and is built from _node_map in
    // prepare(). Both are empty/identity-filled automatically if
    // set_node_map() is never called.
    std::vector<int32_t> _node_map;
    std::vector<int32_t> _row_to_node;

    // Scratch vectors reused by solve() on every level, sized once in
    // prepare() and reused across calls (avoids per-call heap churn).
    mutable std::vector< std::vector<double> > _x, _b, _r;

    // Cached state for the cheap level0->level1 Galerkin update (see
    // build_delta_update_tables()/construct() in gmg_precond.cpp for the
    // full derivation). Only level 0's diagonal ever changes between
    // construct() calls (the off-diagonal/geometric part is fixed for
    // the lifetime of this preconditioner, inherited from
    // EpotMatrixSolver's linear/nonlinear split), so the resulting
    // change to level 1's Galerkin-coarsened operator is exactly
    // R0*diag(delta)*P0 -- a cheap, precomputed scatter-accumulate,
    // rather than redoing the general (and, for the finest hop,
    // dominant-cost) sparse-sparse Galerkin product every iteration.
    std::vector<double>  _lvl0_diag_ref;  // level 0 diagonal as of prepare(), indexed by fine row m
    std::vector<double>  _lvl1_lin_val;   // level 1 operator values with delta==0, i.e. as built by galerkin_coarsen() in prepare(); same layout as _level[1].val
    std::vector<int32_t> _delta_src;      // fine row m for each precomputed triple
    std::vector<int32_t> _delta_dst;      // matching index into _level[1].val
    std::vector<double>  _delta_coef;     // contribution is _delta_coef[k] * (current diag[m] - _lvl0_diag_ref[m])

    void build_delta_update_tables( void );

    void build_level0( const CRowMatrix &A );
    void refresh_level0_values( const CRowMatrix &A );
    void color_graph( Level &lev, const std::vector<int32_t> *node_map = NULL ) const;
    void find_diagonals( Level &lev ) const;
    bool next_level_size( const Level &fine, uint32_t &cnx, uint32_t &cny, uint32_t &cnz ) const;
    void build_transfer_operators( const Level &fine, Level &coarse,
                                    const std::vector<int32_t> *fine_node_map ) const;
    void galerkin_coarsen( const Level &fine, Level &coarse ) const;

    void relax( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                const std::vector<uint32_t> &rows, double w ) const;
    void relax_serial( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                        const std::vector<uint32_t> &rows, double w ) const;
    void smooth( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                 uint32_t nsweeps, double w ) const;
    void smooth_to_convergence( const Level &lev, std::vector<double> &x, const std::vector<double> &b,
                                 uint32_t max_sweeps, double rtol, double w ) const;
    void residual( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                   std::vector<double> &r ) const;
    double residual_norm2( const Level &lev, const std::vector<double> &x, const std::vector<double> &b,
                            std::vector<double> &r ) const;
    void restrict_vec( const Level &coarse, const std::vector<double> &rfine,
                        std::vector<double> &rcoarse ) const;
    void prolong_add( const Level &coarse, const std::vector<double> &ecoarse,
                       std::vector<double> &xfine ) const;

    void vcycle( uint32_t level, std::vector<double> &x, const std::vector<double> &b ) const;
};

#endif