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

    /*! \brief Return the relative permittivity of whatever mesh node
     *  \a mesh_value describes.
     *
     *  Returns 1.0 (vacuum) unless the node is tagged
     *  SMESH_NODE_ID_PURE_VACUUM with a nonzero solid number in its
     *  lower bits, i.e. it is inside a dielectric solid (see
     *  Geometry::build_mesh_parallel_thread_3d() for how such nodes
     *  get tagged this way instead of as Dirichlet) -- in which case
     *  the solid's Bound (BOUND_DIELECTRIC) supplies epsilon_r. Any
     *  other node type (near-solid, Neumann, Dirichlet) is, by
     *  construction, never inside a dielectric in the current
     *  (grid-aligned dielectric surface) implementation, so it is
     *  always treated as vacuum here -- correct because a real vacuum
     *  node adjacent to a genuine conductor is always reclassified
     *  near-solid rather than staying pure vacuum, and a Dirichlet
     *  boundary's fixed value does not depend on permittivity anyway.
     */
    double node_epsilon_r( uint32_t mesh_value ) const;

    /*! \brief Effective relative permittivity to use for a face
     *  linking a free node (of permittivity \a eps_self) to whatever
     *  \a mesh_value describes.
     *
     *  If the neighbour is itself a plain/dielectric node
     *  (SMESH_NODE_ID_PURE_VACUUM[_FIX]), its own permittivity from
     *  node_epsilon_r() is used, same as any interior face. If the
     *  neighbour is Dirichlet (a fixed-voltage conductor, or a
     *  simulation box edge), there is no material information to be
     *  had from that node -- a Dirichlet node's lower bits identify
     *  *which conductor/edge* it is, not a permittivity -- so \a
     *  eps_self is used instead, i.e. the medium \a self sits in is
     *  assumed to extend right up to the fixed-voltage surface with no
     *  other material intervening. This is what makes a dielectric
     *  flush against a conductor (the "ceramic sitting on a biased
     *  electrode" case) work correctly: using 1.0 (vacuum) here
     *  unconditionally, as an earlier version of this function did,
     *  silently gave the wrong answer whenever a dielectric touched a
     *  conductor -- caught by the 1D two-media capacitor verification
     *  test, not by inspection. Any other neighbour type (near-solid,
     *  Neumann, fine-boundary) falls back to vacuum, consistent with
     *  node_epsilon_r() -- see its doc comment for why that is safe in
     *  the current grid-aligned-dielectric implementation.
     */
    double neighbor_epsilon_r( uint32_t mesh_value, double eps_self ) const;

    /*! \brief Effective face permittivity between two adjacent cells
     *  of relative permittivity \a eps_a and \a eps_b, i.e. the
     *  harmonic mean 2*eps_a*eps_b/(eps_a+eps_b).
     *
     *  This is the standard finite-volume treatment for a coefficient
     *  discontinuity: it is what falls out of modelling the two
     *  half-cells adjacent to the shared face as two permittivities in
     *  series. Reduces to eps_a (=eps_b) when both sides match, so
     *  using it unconditionally on every face -- not just ones that
     *  cross a material boundary -- changes nothing for a
     *  uniform-permittivity region (including plain vacuum, eps=1
     *  throughout, which reproduces every coefficient in
     *  add_vacuum_node() exactly as before dielectric support existed).
     */
    static double face_epsilon( double eps_a, double eps_b );

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
