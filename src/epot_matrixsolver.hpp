/*! \file epot_matrixsolver.hpp
 *  \brief Matrix solver for electric potential problem
 */

/* Copyright (c) 2011,2012,2017,2021 Taneli Kalvas. All rights reserved.
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


#ifndef EPOT_MATRIXSOLVER_HPP
#define EPOT_MATRIXSOLVER_HPP 1


#include "epot_solver.hpp"
#include "crowmatrix.hpp"
#include "mvector.hpp"
#include <vector>
#include <utility>


#define N2D_TYPE_MASK   0x80000000 // 100...
#define N2D_TYPE_FIXED  0x80000000 // 100...
#define N2D_TYPE_FREE   0x00000000 // 000...

#define N2D_INDEX_MASK  0x7FFFFFFF // 011...


/*! \brief Parent class for Matrix-based solvers for Electric potential problem.
 */
class EpotMatrixSolver : public EpotSolver {

protected:

    /*! \brief Class nodes to degrees of freedom mapping.
     *
     *  Uses running numbers starting from 0 to point to vectors with
     *  free variables (electric potential and matrices during solving
     *  for example). If a node is solid interior point or dirichlet
     *  boundary, a fixed negative number corresponding to the solid
     *  is used, i.e. from -1 to -6 for boundaries and starting from
     *  -7 for electrodes.
     */
    class Node2DoF {
        Int3D         _size;          /*!< \brief Size of mesh */
        uint32_t     *_n2d;           /*!< \brief Nodes to degrees of freedom array. */
        
    public:
        
        Node2DoF();
        Node2DoF( Int3D size );
        ~Node2DoF();
	void clear( void );
        void resize( Int3D size );
	
        uint32_t &operator()( int i ) 
            { return( _n2d[i] ); }
        uint32_t &operator()( int i, int j ) 
            { return( _n2d[i+j*_size[0]] ); }
        uint32_t &operator()( int i, int j, int k ) 
            { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
        
        const uint32_t &operator()( int i ) const
            { return( _n2d[i] ); }
        const uint32_t &operator()( int i, int j )  const 
            { return( _n2d[i+j*_size[0]] ); }
        const uint32_t &operator()( int i, int j, int k ) const 
            { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
        
        /*! \brief Print debugging information to os.
         */
        void debug_print( std::ostream &os ) const;
    };

    uint32_t               _dof;           /*!< \brief Degrees of freedom. */
    Node2DoF               _n2d;           /*!< \brief Nodes to degrees of freedom map. */
    CRowMatrix            *_fd_mat;        /*!< \brief Finite Difference matrix. */
    Vector                *_fd_vec;        /*!< \brief Finite Difference vector. */
    Vector                *_d_vec;         /*!< \brief Derivative vector for nonlinear solution. */
    const Vector          *_sol;           /*!< \brief Current solution vector. */

    MeshScalarField       *_epot;
    const MeshScalarField *_scharge;

    /*! \brief Per-row buffers used to make build_mat_vec() thread-safe.
     *
     *  CRowMatrix::construct_add() appends each entry at a single
     *  shared cursor (the next free slot is only known once every
     *  earlier row has finished writing), so it cannot be called
     *  directly from multiple threads working on different rows at
     *  once. Instead, set_link() buffers each row's (column, value)
     *  pairs here -- indexed by dof/row so each row's slot is written
     *  by exactly one thread -- and build_mat_vec() copies them into
     *  _fd_mat via construct_add() afterwards, in a cheap serial
     *  compaction pass that does no stencil math.
     */
    std::vector<std::vector<std::pair<int32_t,double> > > _row_entries;

    /*! \brief True once the linear (geometry-only) part of the matrix
     *  and rhs has been built and cached for the current preprocess()
     *  .. postprocess() bracket.
     *
     *  The matrix stencil (neighbour coefficients, boundary/near-solid
     *  weighting, Dirichlet-neighbour rhs contributions) and the
     *  constant space-charge rhs term depend only on the mesh,
     *  boundary conditions and \a scharge -- none of which change
     *  between Newton iterations of a single subsolve() call. Only the
     *  plasma-dependent rhs term and diagonal (_d_vec) actually depend
     *  on the current solution guess \a X, and those are recomputed
     *  every call by update_nonlinear_node(). Set back to false by
     *  reset_matrix(), so it is naturally rebuilt once per solve.
     */
    bool                    _linear_built;

    /*! \brief Cached constant (epot-independent) part of the rhs
     *  vector, snapshotted right after the linear part is built.
     *  Restored into _fd_vec at the start of every later build_mat_vec()
     *  call so the nonlinear pass only has to add the plasma term on
     *  top, instead of re-deriving it from scratch.
     */
    Vector                 *_fd_vec_base;

    /*! \brief Cached pure-linear diagonal of _fd_mat (J0(a,a) for each
     *  row a), captured right after the linear part is built.
     *
     *  get_resjac() forms the true Jacobian diagonal each Newton
     *  iteration as J0(a,a) - D(a) (see build_mat_vec()'s doc comment
     *  at the top of the file for the R = J0*X - B(X), J = J0 + I*D(X)
     *  convention). With _fd_mat now a cached, persistent matrix
     *  instead of one rebuilt from scratch every call, get_resjac()
     *  can no longer just subtract D(a) from whatever is currently in
     *  the diagonal -- that would accumulate corrections across
     *  iterations instead of re-deriving them from the fixed linear
     *  diagonal. This cache is what lets it set() the diagonal fresh
     *  from J0(a,a) each time instead.
     */
    std::vector<double>     _fd_mat_diag0;

    /*! \brief Cached self_material_at() result for *every* mesh node
     *  (free or fixed alike) -- dielectric solid number, or 0 for
     *  plain vacuum/near-solid/Neumann/conductor -- indexed by the
     *  full flat mesh index (k*size(1)+j)*size(0)+i, not by dof row.
     *
     *  Two call sites need this on a hot path: update_nonlinear_node()
     *  (every Newton iteration and every step-size backtracking
     *  re-evaluation, for every free node's own material) and
     *  vacuum_face_coefficient() (every stencil face of every
     *  vacuum/near-solid/Neumann node, for its *neighbour's* material,
     *  to detect a dielectric interface). self_material_at()'s fast
     *  path only recognizes a node whose raw tag directly carries a
     *  material number -- for every other node (which is to say
     *  essentially all of them, dielectric or not) it falls back to
     *  Geometry::inside(), a real geometric containment query against
     *  every solid. Calling that on demand at either site, repeatedly,
     *  is what caused two severe slowdowns in turn: first
     *  update_nonlinear_node() doing it for ~all DOF on every Newton
     *  iteration, then (after that was cached per major cycle)
     *  vacuum_face_coefficient() doing it for ~all DOF times up to 6
     *  faces, redone from scratch every major cycle even though the
     *  geometry -- and therefore every node's material -- never
     *  changes after the mesh is built.
     *
     *  So unlike _linear_built's cache (which legitimately must be
     *  redone every major cycle, since the rhs depends on the current
     *  space charge), this one is populated exactly once, ever, the
     *  first time build_mat_vec() runs -- see _dielectric_mat_built --
     *  and never invalidated afterward, since nothing here depends on
     *  scharge, epot, or the FIXED bit (material_number() already
     *  handles the transiently-fixed _FIX tag variants equivalently to
     *  their unfixed counterparts, so the cached value is correct
     *  whether or not the node happens to be force-pinned when read).
     */
    std::vector<uint32_t>   _dielectric_mat;

    /*! \brief True once _dielectric_mat has been populated. Deliberately
     *  never reset by reset_matrix() -- see _dielectric_mat's doc
     *  comment for why it has a longer lifetime than _linear_built.
     */
    bool                    _dielectric_mat_built;

    /*! \brief Accumulated time spent in _dielectric_mat's one-time (ever)
     *  population pass -- see build_mat_vec(). IS reset by
     *  reset_build_timing() like the other two, even though the
     *  underlying _dielectric_mat cache itself is never rebuilt: the
     *  timed block is gated by _dielectric_mat_built, so once that
     *  flag is set the block simply never executes again and this
     *  accumulator harmlessly stays at 0 for the rest of the program.
     *  (An earlier version of this code deliberately excluded this
     *  variable from reset_build_timing(), on the reasoning that
     *  "it only happens once, so don't bother resetting it" -- but
     *  that left it holding its one nonzero value forever, so every
     *  subsolve() after the first misleadingly reprinted the
     *  original one-time cost as if it had been paid again. Resetting
     *  it unconditionally, every call, is what actually makes it read
     *  0 after the first subsolve(), which was always the intent.)
     */
    double                  _time_dielectric_cache;

    /*! \brief Accumulated time spent in build_mat_vec()'s one-time-per-
     *  major-cycle linear (geometry+rhs) build -- the part gated by
     *  _linear_built. Reset by reset_build_timing().
     */
    double                  _time_linbuild;

    /*! \brief Accumulated time spent in build_mat_vec()'s nonlinear
     *  (plasma term) update pass, summed across every call within the
     *  current subsolve() (i.e. every Newton round and step-size
     *  backtracking re-evaluation). Reset by reset_build_timing().
     */
    double                  _time_nonlin;

    /*! \brief Constructor.
     */
    EpotMatrixSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotMatrixSolver( Geometry &geom, std::istream &s );

    /*! \brief Return const pointers to the matrix \a A and vector \a
     *  B of the linear problem.
     */
    void get_vecmat( const CRowMatrix **A, const Vector **B );

    /*! \brief Return const pointers to jacobian matrix and residual
     *  vector of the problem to \a J and \a R at \a X.
     */
    void get_resjac( const CRowMatrix **J, const Vector **R, const Vector &X );

    /*! \brief Return true if problem is linear.
     */
    //bool linear( void ) const;

    /*! \brief Load initial solution vector from electric potential.
     */
    void set_initial_guess( const MeshScalarField &epot, Vector &X ) const;

    /*! \brief Load electric potential from solution vector.
     *
     *  Only free nodes are set. Fixed nodes are set during preprocess.
     */
    void set_solution( MeshScalarField &epot, const Vector &X ) const;

    /*! \brief Preprocess.
     *
     *  Modify solid mesh suitable for solver, build n2d map and
     *  linear matrix. Make right-hand-side.
     */
    void preprocess( MeshScalarField &epot, const MeshScalarField &scharge );

    /*! \brief Postprocess.
     *
     *  Return solid mesh back to original state and remove temporary
     *  variables.
     */
    void postprocess( void );

    /*! \brief Reset matrix representation.
     */
    void reset_matrix( void );

    /*! \brief Zero the per-solve build_mat_vec() timing accumulators
     *  (_time_linbuild, _time_nonlin, _time_dielectric_cache). Called
     *  once at the start of each EpotBiCGSTABSolver::subsolve() so the
     *  values reported afterward reflect that solve only, not an
     *  accumulation across the whole program run. _time_dielectric_cache
     *  is included here too even though the cache it measures is only
     *  ever (re)built once: resetting the timer is what makes it read
     *  0 in every subsolve() after the one where the build actually
     *  happened, instead of forever reprinting the original one-time
     *  cost. See _time_dielectric_cache's own doc comment.
     */
    void reset_build_timing( void ) {
	_time_linbuild = 0.0;
	_time_nonlin = 0.0;
	_time_dielectric_cache = 0.0;
    }

    /*! \brief Time (seconds) spent in _dielectric_mat's one-time-ever
     *  population pass. 0 for every subsolve() after the first.
     */
    double time_dielectric_cache( void ) const { return( _time_dielectric_cache ); }

    /*! \brief Time (seconds) spent in this solve's one-time linear
     *  (geometry+rhs) matrix build.
     */
    double time_linbuild( void ) const { return( _time_linbuild ); }

    /*! \brief Time (seconds) spent in this solve's nonlinear (plasma
     *  term) update passes, summed over every Newton round/backtrack.
     */
    double time_nonlin( void ) const { return( _time_nonlin ); }

private:

    /*! \brief Build linear matrix and right-hand-side.
     */
    void build_mat_vec( void );

    void set_link( uint32_t a, uint32_t b, double val );

    /*! \brief Add the epot-dependent plasma contribution (rhs term and
     *  diagonal derivative) for free node \a a at mesh location
     *  (i,j,k)/x. Called once per Newton iteration for every free node,
     *  after the cached linear part has been restored into _fd_vec.
     *  This is the part of the old per-node add_*_node() bodies that
     *  actually depends on the current solution guess (via _sol) and
     *  so cannot be cached across Newton iterations. No-op if the
     *  solver is not using a plasma model.
     */
    void update_nonlinear_node( uint32_t a, uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );

    /*! \brief Return the dielectric solid number a mesh node belongs
     *  to, or 0 if it is not inside any dielectric (real vacuum, or any
     *  non-PURE_VACUUM tag: near-solid, Neumann, Dirichlet, fine-
     *  boundary).
     *
     *  A dielectric solid's interior is tagged SMESH_NODE_ID_PURE_VACUUM
     *  with its solid number (>=7) in the lower bits, instead of
     *  SMESH_NODE_ID_DIRICHLET (see
     *  Geometry::build_mesh_parallel_thread_3d()) -- this recovers that
     *  solid number, or 0 for plain vacuum/anything else. Used both to
     *  get a node's own permittivity (node_epsilon_r()) and to compare
     *  two adjacent nodes' materials (vacuum_face_coefficient()).
     */
    static uint32_t material_number( uint32_t mesh_value );

    /*! \brief Like material_number(), but also correct for a Neumann or
     *  Dirichlet box-wall node whose own mesh tag has been overridden
     *  and so no longer carries its material directly.
     *
     *  A dielectric solid that reaches the simulation box edge has
     *  that boundary layer's tag overwritten with the box's own
     *  condition (see Geometry::override_dielectric_box_boundary_3d())
     *  -- necessarily so, since the box edge's own Dirichlet/Neumann
     *  condition takes priority there, not the dielectric's bulk
     *  equation. But add_neumann_node() still needs to know that
     *  node's *permittivity* to weight its own links correctly
     *  whenever the boundary it sits on isn't the only thing that
     *  varies nearby (e.g. a Neumann side wall that a dielectric slab
     *  varying along a different axis happens to touch) -- material_number()
     *  alone would silently and incorrectly say "vacuum" there, since
     *  the tag itself no longer encodes it. This falls back to asking
     *  Geometry directly, by real-world position \a x, which solid (if
     *  any) the point is inside, bypassing the (possibly-overridden)
     *  mesh tag entirely. Returns 0 for plain vacuum, for a point
     *  outside the geometry, or for a non-dielectric (conductor) solid
     *  -- the last case would indicate a genuinely conflicting/
     *  pathological geometry (a conductor reaching exactly to a
     *  Neumann-tagged wall) that this function does not attempt to
     *  resolve.
     */
    uint32_t self_material_at( uint32_t mesh_value, const Vec3D &x ) const;

    /*! \brief Return the relative permittivity of whatever mesh node
     *  \a mesh_value describes.
     *
     *  1.0 (vacuum) unless material_number() is nonzero, in which case
     *  the solid's Bound (BOUND_DIELECTRIC) supplies epsilon_r.
     */
    double node_epsilon_r( uint32_t mesh_value ) const;

    /*! \brief Effective coefficient for the face linking free node
     *  (\a i, \a j, \a k) -- of permittivity \a eps_self and dielectric
     *  material number \a self_material (0 if plain vacuum) -- to its
     *  neighbour at (\a ni, \a nj, \a nk), reached in direction \a sign
     *  (+1/-1) along axis \a coord (0,1,2) and whose raw mesh tag is
     *  \a neighbor_mesh.
     *
     *  Three cases:
     *
     *  - Neighbour is Dirichlet and is the simulation box edge itself
     *    (boundary number 1-6, not a user-defined solid): there is no
     *    material information at a fixed node -- its lower bits
     *    identify *which* edge, not a permittivity -- so \a eps_self is
     *    used. This is exact, not an approximation: a dielectric node
     *    next to the box edge is, by construction, always exactly one
     *    full cell away from it (see
     *    Geometry::override_dielectric_box_boundary_3d()), so there is
     *    no sub-cell position to resolve.
     *
     *  - Neighbour is Dirichlet and is a real user-defined solid
     *    (boundary/solid number >=7, e.g. an STL electrode): unlike the
     *    box edge, this can sit at any sub-cell distance from \a self.
     *    Bisected against that solid's own geometry
     *    (Geometry::solid_face_frac()) to find the true distance alpha,
     *    then eps_self/alpha -- a single resistor of length alpha*h
     *    through self's own medium, since beyond the conductor surface
     *    is a known constant (its fixed voltage), not another variable
     *    material. Reduces to the old flush-contact assumption exactly
     *    at alpha=1 (a full cell away). This is what makes a dielectric
     *    touching a conductor at an arbitrary (e.g. STL-to-STL) position
     *    work correctly, not just when the two happen to be flush.
     *
     *  - Neighbour's material (material_number() of \a neighbor_mesh,
     *    treating near-solid/Neumann/fine-boundary as plain vacuum, 0)
     *    matches \a self_material: same medium on both sides of this
     *    face, no correction needed, coefficient is simply \a eps_self.
     *
     *  - Neighbour's material (material_number() of \a neighbor_mesh,
     *    treating near-solid/Neumann/fine-boundary as plain vacuum, 0)
     *    matches \a self_material: same medium on both sides of this
     *    face, no correction needed, coefficient is simply \a eps_self.
     *
     *  - Otherwise, a genuine material interface crosses this face.
     *    Its exact sub-cell position is found by bisecting against
     *    whichever side is an actual dielectric solid (self's own
     *    material if self is dielectric, else the neighbour's) via
     *    Geometry::solid_face_frac() -- this is what makes the result
     *    correct for an arbitrary (e.g. STL-imported) surface, not just
     *    one that happens to sit exactly halfway between two mesh
     *    nodes. The two half-cell permittivities are then combined in
     *    series, weighted by that fractional distance (alpha):
     *    1/(alpha/eps_self + (1-alpha)/eps_neighbor); this reduces
     *    exactly to the plain harmonic mean at alpha=0.5 (the
     *    grid-aligned case) and to eps_self (=eps_neighbor) whenever
     *    there is in fact no permittivity jump, regardless of alpha.
     *    Two different dielectric solids directly touching each other
     *    (both self_material and the neighbour's material nonzero and
     *    different) are not fully/rigorously handled -- self's own
     *    material is used for the bisection as a reasonable fallback,
     *    a documented limitation rather than a crash risk.
     *
     *    Geometry::solid_face_frac() always bisects assuming its
     *    (i,j,k) argument sits outside the target solid; when self is
     *    in fact the dielectric being bisected against, the bisection
     *    is done from the neighbour's side instead (walking back with
     *    -sign) and complemented (1-alpha) -- see the .cpp for why the
     *    naive self-side call silently returns a nonsense fraction in
     *    that case.
     */
    double vacuum_face_coefficient( int32_t i, int32_t j, int32_t k,
				     int32_t ni, int32_t nj, int32_t nk,
				     int sign, int coord,
				     uint32_t neighbor_mesh, uint32_t self_material,
				     double eps_self ) const;

    /*! \brief Effective permittivity in series across a face at
     *  fractional distance \a alpha (0,1] from the \a eps_self side,
     *  1/(alpha/eps_self + (1-alpha)/eps_neighbor). See
     *  vacuum_face_coefficient() for the derivation and reduction
     *  checks (alpha=0.5 gives the plain harmonic mean; eps_self ==
     *  eps_neighbor gives that common value regardless of alpha).
     */
    static double face_epsilon_alpha( double eps_self, double eps_neighbor, double alpha );

    /*! \brief O(1) lookup into _dielectric_mat by full mesh coordinates
     *  (i,j,k) -- 0 if that node is not inside a dielectric.
     *
     *  Used by add_near_solid_node_1d/2d/cyl/3d to detect a conductor/
     *  dielectric/vacuum triple point: a near-solid axis with a cached
     *  conductor distance on one side whose *other* side (not visible
     *  to is_solid(), which never recognizes a dielectric) turns out to
     *  be a dielectric-interior cell rather than plain vacuum. See
     *  add_near_solid_node_3d()'s doc comment on that branch for the
     *  full rationale.
     */
    uint32_t dielectric_material_at( uint32_t i, uint32_t j, uint32_t k ) const
        { return( _dielectric_mat[ (k*_geom.size(1)+j)*_geom.size(0)+i ] ); }

    void add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_near_solid_node_1d( uint32_t i, const Vec3D &x );
    void add_near_solid_node_2d( uint32_t i, uint32_t j, const Vec3D &x );
    void add_near_solid_node_cyl( uint32_t i, uint32_t j, const Vec3D &x );
    void add_near_solid_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_near_solid_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_neumann_node_1d( uint32_t i, const Vec3D &x );
    void add_neumann_node_2d( uint32_t i, uint32_t j, const Vec3D &x );
    void add_neumann_node_cyl( uint32_t i, uint32_t j, const Vec3D &x );
    void add_neumann_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_neumann_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );

public:

    /*! \brief Destructor.
     */
    virtual ~EpotMatrixSolver();

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;

    // epot_matrixsolver.hpp, public section
    void build_node_map( std::vector<int32_t> &map ) const;
};


#endif
