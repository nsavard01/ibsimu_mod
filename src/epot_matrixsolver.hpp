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

    /*! \brief Inverse of _n2d for free nodes: row -> mesh node index.
     *
     *  _n2d holds a row number for a free node and a mesh node index for
     *  a fixed one, so the two must not be mixed in index arithmetic.
     */
    std::vector<int32_t>   _dof2node;
    CRowMatrix            *_fd_mat;        /*!< \brief Finite Difference matrix. */
    Vector                *_fd_vec;        /*!< \brief Finite Difference vector. */

    /*! \brief Electron confinement mask of the compensation model's
     *  ascent variant: 1 where an electron starting at rest is confined, 0
     *  where it reaches a solid or leaves the compensation region.
     *
     *  Covers the node window _cm_i0.._cm_i0+_cm_ni-1 etc. only, the
     *  bounding box of the compensation region. Nodes outside are not
     *  confined. Made by update_comp_mask().
     */
    std::vector<uint8_t>   _comp_mask;
    uint32_t               _cm_i0, _cm_j0, _cm_k0;
    uint32_t               _cm_ni, _cm_nj, _cm_nk;

    /*! \brief Make _comp_mask for the current solution.
     *
     *  From every node the potential is followed uphill, which is where a
     *  cold electron goes. The node is confined if the path ends at a
     *  potential maximum inside the compensation region.
     */
    void update_comp_mask( void );

public:

    /*! \brief Copy the confinement mask to \a out, 1.0 = confined.
     *
     *  This is the mask that the last solve used, so a diagnostic of the
     *  electron density has to use it too. Empty before the first solve
     *  with the ascent variant.
     */
    void get_comp_mask( MeshScalarField &out ) const;

protected:

    /*! \brief Solution dependent part of the jacobian diagonal, the D(X)
     *  in J = J0 - I*D(X). Written by update_nonlinear_node().
     */
    Vector                 _d_vec;
    const Vector          *_sol;           /*!< \brief Current solution vector. */

    MeshScalarField       *_epot;
    const MeshScalarField *_scharge;

    /*! \brief Linear system, kept from one solve to the next.
     *
     *  The node numbering (_n2d, _dof2node), the matrix J0 (_fd_mat with
     *  the diagonal _fd_mat_diag0) and the boundary part of the right hand
     *  side (_fd_vec_bc) depend neither on the solution nor on the space
     *  charge. build_linear() makes them, and preprocess() calls it only
     *  when linear_inputs_changed() reports a change or after
     *  reset_matrix().
     *
     *  Per solve, preprocess() adds the space charge term to get
     *  _fd_vec_base. Per Newton iteration, build_mat_vec() adds the plasma
     *  terms on top of that and get_resjac() writes the jacobian diagonal
     *  J0(a,a) - D(a) into _fd_mat.
     */
    Vector                  _fd_vec_bc;       /*!< \brief Boundary part of rhs. */
    Vector                  _fd_vec_base;     /*!< \brief Boundary and space charge part of rhs. */
    std::vector<uint8_t>    _row_scharge;     /*!< \brief Does the row take a space charge term? */
    std::vector<double>     _fd_mat_diag0;    /*!< \brief Diagonal of J0. */
    std::vector<int>        _fd_mat_diag_idx; /*!< \brief Index of each row's diagonal in the value array of _fd_mat. */

    /*! \brief State the linear system was built from, see
     *  linear_inputs_changed().
     */
    uint32_t                _key_revision;    /*!< \brief Geometry revision. */
    std::vector<uint32_t>   _key_mesh;        /*!< \brief Node tags. */
    std::vector<double>     _key_epot;        /*!< \brief Potential of fixed nodes. */

    /*! \brief Row entries collected by set_link() during build_linear(),
     *  empty otherwise.
     */
    std::vector<std::vector<std::pair<int32_t,double> > > _row_entries;

    double                  _time_linbuild;   /*!< \brief Time spent in build_linear(). */
    double                  _time_nonlin;     /*!< \brief Time spent in the plasma term update. */

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
     *  Modify solid mesh suitable for solver, build the linear system if
     *  it is not current and make the constant right-hand-side. Returns
     *  true if the linear system was rebuilt, in which case anything
     *  derived from the matrix (a preconditioner) is stale.
     */
    bool preprocess( MeshScalarField &epot, const MeshScalarField &scharge );

    /*! \brief Postprocess.
     *
     *  Return solid mesh back to original state.
     */
    void postprocess( void );

    /*! \brief Drop the linear system. It is rebuilt at the next solve.
     */
    void reset_matrix( void );

    /*! \brief Zero the timers below. Called at the start of a solve.
     */
    void reset_build_timing( void ) {
	_time_linbuild = 0.0;
	_time_nonlin = 0.0;
    }

    /*! \brief Time in seconds spent in build_linear() in this solve.
     */
    double time_linbuild( void ) const { return( _time_linbuild ); }

    /*! \brief Time in seconds spent in the plasma term updates of this
     *  solve.
     */
    double time_nonlin( void ) const { return( _time_nonlin ); }

private:

    /*! \brief Has the geometry, a node tag or a fixed potential changed
     *  since the previous call?
     */
    bool linear_inputs_changed( void );

    /*! \brief Build node numbering, linear matrix and boundary part of
     *  right-hand-side.
     */
    void build_linear( void );

    /*! \brief Set matrix to J0 and build right-hand-side at the current
     *  solution.
     */
    void build_mat_vec( void );

    /*! \brief Add the term \a val * phi(\a b) to the equation of free node \a a,
     *  where \a b is an entry of _n2d. A link into a Neumann mask is
     *  redirected to the mirror node.
     */
    void set_link( uint32_t a, uint32_t b, double val );

    /*! \brief Neumann-mask solids are implemented here, via set_link(). */
    virtual bool supports_neumann_mask( void ) const { return( true ); }

    /*! \brief Add the plasma and compensation terms of free node \a a
     *  at the current solution.
     */
    void update_nonlinear_node( uint32_t a, uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );

    /*! \brief Relative permittivity at node (\a i, \a j, \a k): that
     *  of the dielectric solid the node is in, otherwise 1.
     */
    double node_epsilon_r( uint32_t i, uint32_t j, uint32_t k ) const;

    /*! \brief Weight of the link from node \a i to its neighbour in
     *  direction \a sign (-1/+1) on axis \a coord.
     *
     *  The node has relative permittivity \a eps_self and is inside
     *  dielectric solid \a self_material (0 = vacuum). The weight is the
     *  permittivity of the face between the nodes, in units where a link
     *  between two vacuum nodes is 1:
     *
     *  - same medium at both nodes, or the neighbour on a box face or in
     *    a Neumann mask: eps_self.
     *  - neighbour inside an electrode whose surface is at alpha*h:
     *    eps_self/alpha.
     *  - a material interface at alpha*h: the two media in series,
     *    face_epsilon_alpha().
     *
     *  The surface positions come from Geometry::solid_face_frac(). The
     *  1/r factors of a cylindrical geometry are not included.
     */
    double vacuum_face_coefficient( const int32_t i[3], int coord, int sign,
				     uint32_t self_material, double eps_self ) const;

    /*! \brief Permittivity of two media in series,
     *  1/(alpha/eps_self + (1-alpha)/eps_neighbor), the interface being
     *  at fraction \a alpha of the way from this node to the neighbour.
     */
    static double face_epsilon_alpha( double eps_self, double eps_neighbor, double alpha );

    /* Stencil of a free node (i,j,k) at x, by node type. Each adds the
     * node's matrix row and the boundary part of its right hand side.
     */
    void add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_near_solid_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    void add_neumann_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x );
    double near_solid_axis( uint32_t a, const int32_t i[3], int coord,
			    uint8_t bnd, uint8_t cut, double alpha, double beta, const Vec3D &x );
    double near_solid_radial( uint32_t a, const int32_t i[3],
			      uint8_t bnd, uint8_t cut, double alpha, double beta, const Vec3D &x );

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

    /*! \brief Make the mesh node index -> row map, -1 for fixed nodes.
     */
    void build_node_map( std::vector<int32_t> &map ) const;
};


#endif
