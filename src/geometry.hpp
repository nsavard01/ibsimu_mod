/*! \file geometry.hpp
 *  \brief %Geometry definition
 */

/* Copyright (c) 2005-2013,2019 Taneli Kalvas. All rights reserved.
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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1


#include <stdint.h>
#include <vector>
#include <iostream>
#include "file.hpp"
#include "vec3d.hpp"
#include "vtriangle.hpp"
#include "solid.hpp"
#include "mesh.hpp"
#include "types.hpp"
#include "callback.hpp"


/*! \brief Boundary condition definition class.
 *
 *  Contains boundary condition type and numerical boundary value or a
 *  pointer to a callback functor providing the boundary value as a
 *  function of coordinates \a (x,y,z).  Dirichlet here means fixed,
 *  preset potential at the boundary \f[ \phi = \phi_0. \f] Neumann
 *  here means that the first derivative of the potential with respect
 *  to the unit outward normal (out of solid into the vacuum) of the
 *  surface is preset \f[ - \frac{\partial \phi}{\partial \vec{n}} = -
 *  \sum_i n_i \frac{\partial \phi}{\partial x_i} = q_0. \f]
 *
 *  Dielectric (solids only, n >= 7) is different in kind from the
 *  other two: it is not a boundary *value* at all, so \a value() is
 *  repurposed to carry the solid's relative permittivity \f$
 *  \epsilon_r \f$ instead of a potential. The solid's interior is not
 *  eliminated from the linear system the way a Dirichlet solid's is --
 *  it remains a free region of the mesh satisfying Laplace's equation,
 *  and the permittivity only enters the discretization at the
 *  boundary between differing materials (see EpotMatrixSolver). Must
 *  be a constant value (is_constant() == true); a spatially-varying
 *  permittivity within a single solid is not supported.
 */
class Bound 
{
    bound_e                   _type;    /*!< \brief Boundary type. */
    double                    _value;   /*!< \brief Boundary value if constant. */
    const CallbackFunctorD_V *_functor; /*!< \brief Value functor, NULL if constant. */

    /*! \brief Out-of-line cold paths for the inline value() accessors
     *  above, so that the common constant-value case inlines down to a
     *  null test plus a load without dragging either the exception
     *  machinery or the CallbackFunctorD_V definition into every
     *  translation unit that merely reads a permittivity.
     */
    void throw_nonconstant( void ) const;
    double functor_value( const Vec3D &x ) const;

public:

    /*! \brief Constructor for constant value boundary.
     */
    Bound( bound_e type, double value );

    /*! \brief Constructor for varying value boundary.
     */
    Bound( bound_e type, const CallbackFunctorD_V *functor );

    /*! \brief Constructor for loading boundary condition from a file.
     */
    Bound( std::istream &is );

    /*! \brief Return boundary type.
     *
     *  Inline: Geometry::dielectric_material_at() calls this on every
     *  lookup, and that in turn runs per free node per Newton iteration
     *  and per stencil face during the matrix build.
     */
    bound_e type( void ) const { return( _type ); }

    /*! \brief Set constant boundary value.
     */
    void set_value( double value );

    /*! \brief Return constant boundary value.
     *
     *  This function works only if boundary value is
     *  constant. Otherwise throws an error.
     */
    double value( void ) const {
	if( _functor )
	    throw_nonconstant();
	return( _value );
    }

    /*! \brief Return boundary value at \a x.
     */
    double value( const Vec3D &x ) const {
	if( _functor )
	    return( functor_value( x ) );
	return( _value );
    }

    /*! \brief Return if boundary value is constant.
     *
     *  Returns true if boundary value is constant and false if it is
     *  a function of location.
     */
    bool is_constant() const { return( !_functor ); }

    /*! \brief Saves data to stream \a os.
     */
    void save( std::ostream &os ) const;

    /*! \brief Outputting to stream.
     */
    friend std::ostream &operator<<( std::ostream &os, const Bound &b );
};


#define SMESH_NODE_ID_MASK             0xE0000000 // 111...

#define SMESH_NODE_ID_PURE_VACUUM      0x00000000 // 000...
#define SMESH_NODE_ID_NEAR_SOLID       0x20000000 // 001...
#define SMESH_NODE_ID_NEUMANN          0x40000000 // 010...
#define SMESH_NODE_ID_ROUGH_BOUNDARY   0x60000000 // 011...

#define SMESH_NODE_ID_PURE_VACUUM_FIX  0x80000000 // 100...
#define SMESH_NODE_ID_NEAR_SOLID_FIX   0xA0000000 // 101...
#define SMESH_NODE_ID_DIRICHLET        0xC0000000 // 110...
#define SMESH_NODE_ID_FINE_BOUNDARY    0xE0000000 // 111...

/*! \brief Node inside a Neumann-mask solid: removed from the solve, with
 *  a zero-flux (homogeneous Neumann) condition on every face separating it
 *  from a live node.
 *
 *  Takes the FINE_BOUNDARY slot, which is declared above but referenced
 *  nowhere else in the library -- so this costs no re-encoding. Sharing the
 *  bit pattern is safe precisely because nothing ever tests for
 *  FINE_BOUNDARY; if that ever changes, the two must be split and one of
 *  them given a different id.
 *
 *  The SMESH_NODE_FIXED bit is set (0xE0000000 & 0x80000000 != 0), so
 *  EpotMatrixSolver::preprocess() eliminates these nodes from the matrix
 *  without any change there. What makes them a Neumann mask rather than a
 *  Dirichlet block is EpotMatrixSolver::set_link(), which redirects a link
 *  pointing AT one of these onto the row's own diagonal instead of moving
 *  it to the right-hand side -- see the comment there for why that is
 *  exactly a zero-flux face.
 *
 *  A masked solid is created by giving a solid BOUND_NEUMANN, which
 *  set_boundary() accepts only for n >= 7.
 */
#define SMESH_NODE_ID_NEUMANN_MASK     0xE0000000 // 111...

/*! \brief True if this node is inside a Neumann-mask solid. */
#define SMESH_NODE_IS_NEUMANN_MASK(node) \
    ( ((node) & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN_MASK )

#define SMESH_NODE_FIXED               0x80000000 // 100...


#define SMESH_BOUNDARY_NUMBER_MASK     0x000000FF // limit to 0-255
#define SMESH_NEAR_SOLID_INDEX_MASK    0x1FFFFFFF // limit to 0-2^29 (5.4e8)


/*! \brief True if raw smesh value \a node is physically solid material
 *  -- impenetrable to particles and part of the visible/triangulated
 *  surface -- regardless of whether the potential solver eliminates it
 *  (a Dirichlet conductor) or keeps it as a free node satisfying
 *  Laplace's equation (a dielectric, see BOUND_DIELECTRIC in
 *  types.hpp).
 *
 *  This is deliberately a *different* question from is_solid()/
 *  is_near_solid(): those exist for the potential solver's mesh
 *  classification and must stay Dirichlet-only (a dielectric's
 *  interior is intentionally not eliminated, and its vacuum-facing
 *  neighbours must stay ordinary PURE_VACUUM nodes rather than being
 *  reclassified near-solid -- see EpotMatrixSolver::add_vacuum_node()).
 *  Surface triangulation (mc_case(), surface_cell_face_case_2d(),
 *  surface_inside_solid_number()) and particle collision detection
 *  (ParticleIterator::get_solid()/is_solid()) ask a geometric question
 *  instead -- "is there solid material here a particle should be
 *  absorbed by / a surface should be drawn at" -- which dielectrics
 *  must answer the same way conductors do.
 */
#define SMESH_NODE_IS_SOLID_MATERIAL(node) \
    ( ( ((node) & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET && \
	((node) & SMESH_BOUNDARY_NUMBER_MASK) >= 7 ) || \
      ( ((node) & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM && \
	((node) & SMESH_NEAR_SOLID_INDEX_MASK) >= 7 ) )


/*! \brief %Geometry defining class.
 *
 *  %Geometry class holds the definitions of the geometry
 *  dimensionality, mesh size and electrode configuration. Also it
 *  contains a signed char array for information about the type of
 *  each node. This array is known as the solid mesh.
 *
 *  The integer numbers in the solid mesh have the following meanings:
 *  The solid number 0 is reserved for vacuum and solid numbers from 1
 *  to 6 are reserved for Dirichlet type boundaries of the bounding
 *  box. %Solid numbers from -1 to -6 are reserved for Neumann type
 *  boundaries of the bounding box. Negative solid numbers starting
 *  from -7 are used for marking electrode edges and the positive
 *  numbers starting from 7 are used to mark the interior points of
 *  the electrodes.
 *
 *  The mesh nodes are marked using the following logic: First nodes,
 *  which are inside electrodes are marked solid (>=7). If a point is
 *  inside several solids, the highest solid number is marked. Other
 *  points are left as vacuum nodes (0). As the next step the solid
 *  nodes are mapped to find nodes which have vacuum as closest
 *  neighbour along any of the axes (not diagonal). These nodes are
 *  marked as solid edges (<=-7). As the last step, the vacuum nodes
 *  at the simulation box boundary are marked either as Neumann (<0
 *  and >-7) or Dirichlet (>0 and <7).
 *
 *  Starting from 1.0.3: A. Mark solids, B. Mark Neumann and Dirichlet
 *  boundaries, C. Mark edges taking in account that Neumann = Vacuum
 *  and Dirichlet != Vacuum.
 *
 *  Bounding box edges are numbered in order xmin, xmax, ymin, ymax,
 *  zmin, xmax.
 */
class Geometry : public Mesh
{
    uint32_t                   _n;         /*!< \brief Number of solids */
    std::vector<const Solid*>  _sdata;     /*!< \brief Array of solid definitions, size \a _n */
    std::vector<Bound>         _bound;     /*!< \brief Array of boundary conditions, size \a _n+6 */

    bool                       _built;     /*!< \brief Is solid mesh array built? */
    uint32_t                  *_smesh;     /*!< \brief Solid mesh array. */
    std::vector<uint8_t>       _nearsolid; /*!< \brief Near solid data. */

    /*! \brief Real solid number (0 = vacuum, else the solid -- conductor
     *  or dielectric alike -- this node's location truly belongs to),
     *  for *every* node, independent of and never touched by any
     *  stencil/tag reclassification (near-solid, box-edge overrides,
     *  the _FIX bit, etc.).
     *
     *  _smesh's own tag is fundamentally a single mutually-exclusive
     *  "what stencil applies here" category (plain free node, near a
     *  conductor, box Neumann, box Dirichlet...), and only *sometimes*
     *  also carries a material number in its unused low bits, for
     *  whichever category happened to claim the node first. That's
     *  fine as long as a node only ever needs one of those meanings at
     *  once, but it breaks at a genuine triple point -- a dielectric
     *  solid that reaches the simulation box edge *and* has a real
     *  conductor immediately adjacent there needs to be NEAR_SOLID-
     *  tagged (to carry that conductor's cached bracket-distance data)
     *  while *also* remembering its own dielectric material for its
     *  own permittivity -- two payloads the same node needs at once,
     *  which no single category tag can hold no matter which spare bit
     *  pattern is used for it.
     *
     *  This array sidesteps that by never being reclassified at all:
     *  it's populated once, directly from solid membership, at the
     *  same point build_mesh_parallel_thread_3d()/_2d()/_1d() already
     *  tests each node against every solid to build the initial tag --
     *  so filling it in costs one more array write per hit, no extra
     *  geometry queries. override_dielectric_box_boundary_3d()/_2d()/
     *  _1d() and the near-solid/box-edge classification passes still
     *  run exactly as before (they still have their own, unrelated job:
     *  deciding which *stencil* applies at a node) but none of them
     *  touch this array, so a node's real material is always readable
     *  here regardless of what its _smesh tag ends up being reclassified
     *  to. See material()/material_check().
     *
     *  Serialized directly by Geometry::save()/the istream constructor,
     *  rather than being re-derived from _sdata on load: save_solids is
     *  false in every real caller in this codebase, so a loaded Geometry's
     *  _sdata is all NULL placeholders and there is no live solid geometry
     *  left to re-query via inside() by the time a stream constructor
     *  would need to rebuild this array.
     *
     *  This used to be duplicated ad hoc wherever a caller needed
     *  material identity independent of tag reclassification --
     *  EpotMatrixSolver's own _dielectric_mat cache (populated
     *  reactively via a Geometry::inside() fallback query per node, the
     *  first time it was needed) being the most complete prior version.
     *  Promoting the same idea here, populated for free at mesh-build
     *  time instead of reactively, lets every consumer (the potential
     *  solver, scharge cleanup, surface triangulation, particle
     *  collision) share one already-correct answer instead of each
     *  needing its own cache or living with the tag's blind spots.
     */
    std::vector<uint8_t>       _material;

    double                     _surface_eps; /*!< \brief Vertec matching tolerance. */
    VTriangleSurface           _surface;   /*!< \brief Triangulated surface. */
    std::vector<int32_t>       _triptr;    /*!< \brief Pointer from mesh cell to first triangle. */
    
    /*! \brief Check if node is solid (n>=7).
     *
     *  Returns 0 if node is not solid (vacuum, outside mesh or
     *  simulation box boundary). If node is solid, the solid number
     *  >= 7 is returned.
     */
    uint32_t is_solid( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Add an entry to near solid data.
     *
     *  Appends a new entry to near solid data and updated \a
     *  near_solid_index to point to the next byte after the end of
     *  data.
     */
    void add_near_solid_entry( uint32_t &near_solid_index, int32_t i, int32_t j, int32_t k );

    /*! \brief Bracket solid surface between node points
     *
     *  Finds surface of solid \a solid between the outer node (\a i,
     *  \a j, \a k) and the inner node in \a sign direction (+/-1) of
     *  of coordinate axis \a coord (0,1 or 2). Returns parametric
     *  distance (between 1 and 255) from the outer node to the inner
     *  node.
     */
    uint8_t bracket_ndist( int32_t i, int32_t j, int32_t k, int32_t solid, int sign, int coord ) const;

    /*! \brief Check mesh definition validity.
     */
    void check_definition();

    Vec3D surface_normal_2d( const Vec3D &x ) const;
    Vec3D surface_normal_3d( const Vec3D &x ) const;


    void build_mesh_parallel_near_solid( uint32_t ind, int32_t i, int32_t j, int32_t k );
    void build_mesh_parallel_prepare_near_solid( uint32_t &near_solid_index,
						 int32_t i, int32_t j, int32_t k );
    void build_mesh_parallel_prepare_3d( void );
    void build_mesh_parallel_prepare_2d( void );
    void build_mesh_parallel_prepare_1d( void );

    /*! \brief Reset any of the six (3d)/four (2d)/two (1d) boundary
     *  faces that a dielectric solid's interior has claimed back to
     *  "unclaimed", so the box's own boundary condition (or a real
     *  conductor solid, if one happens to be adjacent right there) can
     *  be applied to it correctly.
     *
     *  A dielectric solid's interior is tagged SMESH_NODE_ID_PURE_VACUUM
     *  (see build_mesh_parallel_thread_3d()), so if it happens to reach
     *  the simulation box edge, build_mesh_parallel_prepare_3d()'s
     *  box-edge classification pass -- which only classifies still-
     *  unclaimed (mesh==0) nodes -- skips it, leaving that boundary
     *  layer tagged as ordinary free dielectric interior instead of the
     *  box's actual Dirichlet/Neumann condition. That is both physically
     *  wrong (the box condition, not the dielectric's bulk equation,
     *  applies exactly at that surface) and a crash: add_vacuum_node()
     *  would read one mesh node past the domain edge trying to reach a
     *  neighbour that doesn't exist.
     *
     *  Called once, *before* build_mesh_parallel_prepare_*(), this walks
     *  only the boundary faces (cheap, O(surface area) not O(volume))
     *  and resets any face node still tagged as dielectric back to
     *  mesh==0, then leaves the actual reclassification to
     *  build_mesh_parallel_prepare_*()'s own decision tree. That tree
     *  already checks is_near_solid() before any box-edge condition, in
     *  the correct xmin/xmax/ymin/ymax/zmin/zmax priority order, so a
     *  reset node is reclassified exactly as an ordinary free node in
     *  that position would be -- including the case where a *real
     *  conductor* solid is also immediately adjacent to it, which must
     *  become a properly cached NEAR_SOLID node, not a plain box-edge
     *  tag. An earlier version of this function ran *after* prepare_*()
     *  and duplicated its box-edge decision tree directly, which bypassed
     *  is_near_solid() entirely and left such a node neither NEAR_SOLID-
     *  tagged nor dielectric-tagged -- Geometry::solid_dist() would then
     *  throw "not a near solid node" the first time surface
     *  triangulation crossed that edge (hit on a real 5-solid STL
     *  geometry where a dielectric and a conductor solid both reached
     *  the same box-edge cell). Interior dielectric nodes are untouched.
     */
    void override_dielectric_box_boundary_3d( void );
    void override_dielectric_box_boundary_2d( void );
    void override_dielectric_box_boundary_1d( void );

    /*! \brief Narrow [imin,imax] (inclusive) along one mesh axis to
     *  the node range that can overlap world-space bounds [bmin,bmax],
     *  clamped to the existing mesh axis. Leaves imin/imax untouched
     *  (full axis) if bmin/bmax indicate the axis is unrestricted --
     *  see Solid::bbox_from_local_box() for that convention. If the
     *  bounds don't overlap the mesh at all along this axis, sets
     *  imin/imax so the range is empty (imin > imax).
     */
    void restrict_axis_range( double bmin, double bmax, double origo, int32_t size,
			      int32_t &imin, int32_t &imax ) const;

    /*! \brief Compute the inclusive node index range [imin,imax] x
     *  [jmin,jmax] x [kmin,kmax] that a solid's bounding box can
     *  possibly touch, defaulting to the whole mesh on any axis (or
     *  entirely) where Solid::get_bbox() doesn't provide a usable
     *  bound. Used to skip calling solid->inside() on mesh nodes that
     *  provably cannot be inside the solid.
     */
    void solid_node_range( const Solid *solid,
			   int32_t &imin, int32_t &imax,
			   int32_t &jmin, int32_t &jmax,
			   int32_t &kmin, int32_t &kmax ) const;

    void build_mesh_parallel_thread_3d( void );
    void build_mesh_parallel_thread_2d( void );
    void build_mesh_parallel_thread_1d( void );

    void build_mesh_parallel( void );

    static const int32_t mc_faces[15*256];

    Vec3D mc_surface( int32_t i, int32_t j, int32_t k, uint8_t cn, int32_t ei ) const;
    uint8_t mc_case( int32_t i, int32_t j, int32_t k ) const;
    uint32_t mc_trianglec( uint8_t cn ) const;
    int32_t mc_add_vertex_try( const Vec3D &x, int32_t i, int32_t j, int32_t k ) const;
    uint32_t mc_add_vertex( const Vec3D &x, int32_t i, int32_t j, int32_t k, uint32_t firstv );
    void mc_triangulate( int32_t i, int32_t j, int32_t k );

    uint8_t surface_cell_face_case_2d( const int32_t i[3], const int32_t vb[3] ) const;
    double surface_cell_face_dist( const int32_t i[3], const int32_t vb[3], 
				   int32_t dx, int32_t dy, 
				   int32_t dir ) const;
    uint32_t surface_cell_face_add_vertex( VTriangleSurfaceSolid &solid,
					   const int32_t i[3], const int32_t vb[3], 
					   double dx, double dy ) const;
    void surface_cell_face_add_triangles( VTriangleSurfaceSolid &solid,
					  int32_t i0, int32_t i1, int32_t i2, 
					  int32_t vb0, int32_t vb1, int32_t vb2 ) const;
    void surface_cell_face_add_triangles( VTriangleSurfaceSolid &solid,
					  const int32_t i[3], const int32_t vb[3] ) const;
    uint32_t surface_inside_solid_number( int32_t i, int32_t j,int32_t k ) const;

public:

    /*! \brief Constructor for geometry class.
     *
     *  Sets geometry mode, mesh cell size \a h, mesh size \a size and
     *  origo \a origo.
     */
    Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h );

    /*! \brief Constructor for loading geometry from a stream \a is.
     */
    explicit Geometry( std::istream &is );

    /*! \brief Destructor for geometry.
     */
    ~Geometry();

    /*! \brief Return number of solids.
     */
    uint32_t number_of_solids() const;

    /*! \brief Return number of boundaries.
     *
     *  Always returns at least \a n >= 6.
     */
    uint32_t number_of_boundaries() const;
    
    /*! \brief Sets solid number \a n to \a s.
     *
     *  Solids have to be defined in incresing order. %Solid number \a
     *  n should be >= 7. This function can also be used to overwrite
     *  a previous solid definition. Pointer to solid \a s is saved to
     *  geometry structure. %Solid will not be deleted when geometry is
     *  deleted. The newly defined defined solids default to Dirichlet
     *  boundary with potential zero.
     */
    void set_solid( uint32_t n, const Solid *s );

    /*! \brief Returns a const pointer to solid number \a n.
     *
     *  %Solid number \a n should be >= 7.
     */
    const Solid *get_solid( uint32_t n ) const;

    /*! \brief Sets boundary condition \a b for solid number \a n.
     *
     *  %Solid number \a n should be > 0 here. Boundary numbers from 1
     *  to 6 are the boundary conditions for the bounding box. Numbers
     *  starting from 7 are the user defined solids. All boundaries of
     *  the simulation box (n <= 6) default to Neumann boundary
     *  condition with derivative value zero. All defined solids (n >=
     *  7) default to Dirichlet boundary with potential zero.
     *
     *  In cylindrical geometry case the \a rmin boundary can be set
     *  to Dirichlet, which means that there is an infinitely thin
     *  wire with a fixed potential at the axis or to Neumann, which
     *  means that the natural boundary for cylindrical axis will be
     *  used.
     */
    void set_boundary( uint32_t n, const Bound &b );

    /*! \brief Returns boundary condition for solid number \a n.
     */
    Bound get_boundary( uint32_t n ) const;

    /*! \brief Returns a vector of boundary conditions.
     */
    std::vector<Bound> get_boundaries() const;

    /*! \brief Returns true if full solid data available.
     *
     *  Returns false if any defined solid does not include solid data.
     */
    bool have_solid_data( void ) const;

    /*! \brief Return if point is inside solids.
     *
     *  Returns 0 if point \a x is in vacuum or the number of
     *  solid of \a x is inside a defined solid. Returns a number from
     *  1 to 6 if point \a x is outside the defined geometry. If the
     *  point is inside several defined solids, the solid with the
     *  highest solid number is returned.
     *
     *  Uses solid data.
     */
    uint32_t inside( const Vec3D &x ) const;

    /*! \brief Returns true if point \a x is inside solid \a n.
     *
     *  Uses solid data.
     */
    bool inside( uint32_t n, const Vec3D &x ) const;

    /*! \brief Fractional distance (0,1] from mesh node (\a i, \a j, \a
     *  k) to the surface of solid \a solid, walking in direction \a
     *  sign (+1 or -1) along axis \a coord (0, 1 or 2). A return value
     *  of 1.0 means the surface is (at least) a full cell width \a h
     *  away in that direction -- i.e. no correction is needed -- and
     *  smaller values mean the surface sits closer to node (i,j,k) than
     *  a full cell width.
     *
     *  This is a public wrapper around the same bisection used by the
     *  conductor near-solid mechanism (see is_near_solid()/bracket_ndist()),
     *  exposed for EpotMatrixSolver::add_vacuum_node() to place a
     *  dielectric material interface at its true sub-cell position
     *  (needed for an arbitrary, e.g. STL-imported, surface) rather than
     *  assuming it always falls exactly halfway between two mesh nodes
     *  (which is only true for a grid-aligned interface). Unlike the
     *  conductor case, a dielectric interior stays an ordinary free
     *  (PURE_VACUUM-tagged) node rather than being reclassified
     *  near-solid, so this fractional distance is computed on demand
     *  here instead of being cached in _nearsolid.
     */
    double solid_face_frac( int32_t i, int32_t j, int32_t k, uint32_t solid, int sign, int coord ) const;

    /*! \brief Find solid \a n surface location by bracketing.
     *
     *  Searches for the solid \a n surface location on the line
     *  between points \a xin and \a xout by bracketing. Point \a xin
     *  should be inside the solid and point \a xout should be outside
     *  the solid. Function saves the coordinates of the surface to
     *  xsurf and returns parametrical distance (value from 0 to 1)
     *  from xin.
     *
     *  Uses solid data.
     */
    double bracket_surface( uint32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const;

    /*! \brief Find surface outward normal at location \a x.
     *
     *  Returns zero vector on failure.  Uses solid data. The location
     *  \a x should be within 1 mesh unit of a surface.
     */
    Vec3D surface_normal( const Vec3D &x ) const;

    /*! \brief Is the solid mesh built?
     */
    bool built( void ) const { return( _built ); }

    /*! \brief Builds (or rebuilds) the solid mesh and near solid data
     *  from solid definitions.
     */
    void build_mesh( void );

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i ) const { return( _smesh[i] ); }

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i, int32_t j ) const {
	return( _smesh[i + j*_size[0]] ); 
    }

    /*! \brief Returns a const reference to solid mesh array.
     */
    const uint32_t &mesh( int32_t i, int32_t j, int32_t k ) const {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i ) { return( _smesh[i] ); }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i, int32_t j ) {
	return( _smesh[i + j*_size[0]] );
    }

    /*! \brief Returns a reference to solid mesh array.
     */
    uint32_t &mesh( int32_t i, int32_t j, int32_t k ) {
	return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Returns the real solid number node \a i truly belongs to
     *  (0 = vacuum), independent of any stencil/tag reclassification.
     *  See _material's doc comment.
     */
    uint32_t material( int32_t i ) const { return( _material[i] ); }

    /*! \brief Returns the real solid number node (\a i, \a j) truly
     *  belongs to (0 = vacuum), independent of any stencil/tag
     *  reclassification. See _material's doc comment.
     */
    uint32_t material( int32_t i, int32_t j ) const {
	return( _material[i + j*_size[0]] );
    }

    /*! \brief Returns the real solid number node (\a i, \a j, \a k)
     *  truly belongs to (0 = vacuum), independent of any stencil/tag
     *  reclassification. See _material's doc comment.
     */
    uint32_t material( int32_t i, int32_t j, int32_t k ) const {
	return( _material[i + j*_size[0] + k*_size[0]*_size[1]] );
    }

    /*! \brief Like material(), but returns 0 (never solid) instead of
     *  throwing/reading out of bounds for a point outside the mesh --
     *  mirrors mesh_check()'s bounds handling, but since being outside
     *  the box is never "real solid material" there is no equivalent of
     *  mesh_check()'s synthetic Dirichlet boundary number to return.
     */
    uint32_t material_check( int32_t i ) const;

    /*! \brief See material_check( int32_t i ) const.
     */
    uint32_t material_check( int32_t i, int32_t j ) const;

    /*! \brief See material_check( int32_t i ) const.
     */
    uint32_t material_check( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Return the dielectric solid number node (\a i, \a j, \a k)
     *  truly belongs to, or 0 if it's vacuum or a conductor. Reads
     *  material() directly (see _material's doc comment) filtered to
     *  dielectrics only.
     *
     *  Public so every consumer that needs a node's material identity
     *  independent of any stencil/tag reclassification -- the potential
     *  solver (EpotMatrixSolver), surface triangulation (solid_dist(),
     *  below), scharge cleanup -- can share this one answer instead of
     *  each keeping its own copy or cache of the same logic.
     *
     *  Defined inline: EpotMatrixSolver::update_nonlinear_node() calls
     *  this once per free node on every Newton iteration, and
     *  vacuum_face_coefficient() calls it once per face (six times per
     *  node in 3D) during the matrix build, so an out-of-line call --
     *  uninlinable across translation units without LTO -- was pure
     *  overhead around what is otherwise two array reads and a compare.
     */
    uint32_t dielectric_material_at( int32_t i, int32_t j, int32_t k ) const {
	uint32_t solid_number = _material[i + j*_size[0] + k*_size[0]*_size[1]];
	if( solid_number == 0 || _bound[solid_number-1].type() != BOUND_DIELECTRIC )
	    return( 0 ); // vacuum, or a conductor -- not this function's concern
	return( solid_number );
    }

    /*! \brief Non-throwing, non-copying access to boundary \a n's
     *  definition -- the same object get_boundary() returns, but by
     *  const reference and without the range check.
     *
     *  get_boundary() returns Bound BY VALUE and validates \a n on every
     *  call, which is the right interface for user code but makes it a
     *  poor fit for the solver's inner loops (node_epsilon_r() and
     *  vacuum_face_coefficient() ask for a permittivity per node/face):
     *  the copy is pointless there, and the potentially-throwing check
     *  blocks inlining. Callers of this overload are responsible for
     *  passing an \a n they already know is valid -- in practice one
     *  that came out of material()/dielectric_material_at(), which can
     *  only ever return 0 or a real solid number by construction.
     */
    const Bound &boundary( uint32_t n ) const {
	return( _bound[n-1] );
    }

    /*! \brief Returns number from solid mesh array.
     *
     *  For 1D geometries. Returns number from solid mesh array at \a
     *  i or Dirichlet boundary number (1 or 2) if point is outside
     *  mesh.
     */
    uint32_t mesh_check( int32_t i ) const;

    /*! \brief Returns number from solid mesh array.
     *
     *  For 2D geometries. Returns number from solid mesh array at \a
     *  (i,j) or Dirichlet boundary number (1-4) if point is outside
     *  mesh.
     */
    uint32_t mesh_check( int32_t i, int32_t j ) const;

    /*! \brief Returns number from solid mesh array.
     *
     *  Returns number from solid mesh array at \a (i,j,k) or
     *  Dirichlet boundary number (1-6) if point is outside mesh.
     */
    uint32_t mesh_check( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Returns true if node is a potential near solid point.
     *
     *  Returns true if any of the neighbouring points is a solid
     *  point (Dirichlet with solid number >= 7). The state of node \a
     *  (i,j,k) is not checked.
     */
    bool is_near_solid( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Returns a const pointer to start of near-solid data for
     *  near-solid node \a index. 
     *
     *  The first byte contains the bit
     *  flags for the existance of neighbouring solids. From bit 0 to
     *  bit 5 the boolean flags are for directions: xmin, xmax, ymin,
     *  ymax, zmin, zmax. The next bytes contain the parametric
     *  distances of the solid surfaces from the node in each
     *  direction. Only the directions with set bit flag are saved to
     *  data. The distances are saved in the same order as the flags
     *  (from xmin to zmax). The distance information is an unsigned
     *  8-bit integer (0 to 255), where 0 means distance 0.0 and 255
     *  means 1.0.
     */
    const uint8_t *nearsolid_ptr( int32_t index ) const {
	return( &_nearsolid[index] );
    }

    /*! \brief Returns distance of solid boundary from point.
     *
     *  Returns the distance (0 to 255) of solid surface in direction
     *  \a dir from near solid point at (\a i, \a j, \a k). The
     *  direction \a dir is an integer from 0 to 5, with 0 meaning -x,
     *  1 meaning +x, 2 meaning -y, 3 meaning +y, 4 meaning -z and 5
     *  meaning +z. If the node at (\a i, \a j, \a k) is not a near
     *  solid node or if there is no solid nearby in the direction an
     *  error will be thrown.
     */
    uint8_t solid_dist( uint32_t i, uint32_t j, uint32_t k, uint32_t dir ) const;

    /*! \brief Returns distance of solid boundary from point.
     *
     *  Same as solid_dist() above, just using one dimensional index
     *  for mesh.
     */
    uint8_t solid_dist( uint32_t i, uint32_t dir ) const;

    /*! \brief Build surface triangulation data.
     */
    void build_surface( void );
    
    /*! \brief Is the solid surface representation built?
     */
    bool surface_built( void ) const { return( _triptr.size() ); }

    /*! \brief Finds if point is inside surface triangulation.
     *
     *  Returns 0 if point \a x is in vacuum or the number of
     *  solid of \a x is inside a defined solid. Returns a number from
     *  1 to 6 if point \a x is outside the defined geometry. Uses 
     *  solid mesh and surface triangulation data.
     *
     *  The surface triangulation does not separate solids which are
     *  touching (no vacuum node in between). Therefore in these cases
     *  this function can not separate these solids with accuracy
     *  higher than one grid cell. The vacuum to solid separation
     *  always has the maximum resolution.
     */
    uint32_t surface_inside( const Vec3D &x ) const;

    /*! \brief Return total surface vertex count.
     */
    uint32_t surface_vertexc( void ) const {
	return( _surface.vertexc() );
    }

    /*! \brief Return total surface triangle count.
     */
    uint32_t surface_trianglec( void ) const {
	return( _surface.trianglec() );
    }

    /*! \brief Return reference to surface vertex \a a.
     */
    const Vec3D &surface_vertex( int32_t a ) const {
	return( _surface.vertex(a) );
    }

    /*! \brief Return normal of surface triangle near point x.
     *
     *  Uses surface triangulation, returns NaN if surface not found.
     */
    Vec3D surface_triangle_normal( const Vec3D &x ) const;

    /*! \brief Return normal of surface triangle \a a.
     */
    Vec3D surface_triangle_normal( int32_t a ) const;

    /*! \brief Return reference to surface triangle \a a.
     */
    const VTriangle &surface_triangle( int32_t a ) const {
	return( _surface.triangle(a) );
    }

    /*! \brief Return index of first surface triangle at mesh cube \a (i,j,k).
     */
    uint32_t surface_triangle_ptr( int32_t i, int32_t j, int32_t k ) const {
	return( _triptr[(k*(_size[1]-1) + j)*(_size[0]-1) + i] );
    }

    /*! \brief Return surface triangle count at mesh cube \a (i,j,k).
     */
    int32_t surface_trianglec( int32_t i, int32_t j, int32_t k ) const;

    /*! \brief Saves data to a new file \a filename.
     *
     *  If \a save_solids is true, also the solid definitions are
     *  saved. Not all solid types are saveable. A warning is printed
     *  if trying to save such a solid.
     */
    void save( const std::string &filename, bool save_solids = false ) const;

    /*! \brief Saves data to stream \a os.
     *
     *  If \a save_solids is true, also the solid definitions are
     *  saved. Not all solid types are saveable. A warning is printed
     *  if trying to save such a solid.
     */
    void save( std::ostream &os, bool save_solids = false ) const;

    /*! \brief Print debugging information to stream \a os.
     */
    void debug_print( std::ostream &os ) const;
};


#endif

