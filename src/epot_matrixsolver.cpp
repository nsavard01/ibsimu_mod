/*! \file epot_matrixsolver.cpp
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


#include "epot_matrixsolver.hpp"
#include "constants.hpp"
#include "ibsimu.hpp"
#include <omp.h>
#include <chrono>

namespace {
    using Clock = std::chrono::steady_clock;
}


EpotMatrixSolver::Node2DoF::Node2DoF() 
  : _size(0), _n2d(0) 
{

}


EpotMatrixSolver::Node2DoF::Node2DoF( Int3D size ) 
  : _size(size) 
{
    _n2d = new uint32_t[_size[0]*_size[1]*_size[2]];
}


EpotMatrixSolver::Node2DoF::~Node2DoF() 
{
    if( _n2d )
	delete[] _n2d; 
}


void EpotMatrixSolver::Node2DoF::clear( void ) 
{
    if( _n2d )
	delete[] _n2d;
    _n2d = 0;
    _size = Int3D(0,0,0);
}


void EpotMatrixSolver::Node2DoF::resize( Int3D size ) 
{
    _size = size;
    if( _n2d )
	delete[] _n2d;
    _n2d = new uint32_t[_size[0]*_size[1]*_size[2]];
}


void EpotMatrixSolver::Node2DoF::debug_print( std::ostream &os ) const
{
    uint32_t a;
    uint32_t nc = _size[0]*_size[1]*_size[2];
    os << "**Node2DoF\n";
    os << "size = ("
        << _size[0] << ", "
        << _size[1] << ", "
        << _size[2] << ")\n";
    os << "n2d = (";
    if( _n2d ) {
	if( nc <= 10 ) {
	    for( a = 0; a < nc-1; a++ )
		os << _n2d[a] << ", ";
	    os << _n2d[a] << ")\n";
	} else {
	    for( a = 0; a < 10; a++ )
		os << _n2d[a] << ", ";
	    os << " ...)\n";
	}
    } else {
	os << ")\n";
    }
}


EpotMatrixSolver::EpotMatrixSolver( Geometry &geom )
    : EpotSolver(geom), _dof(0), _fd_mat(0), _fd_vec(0),
      _sol(0),
      _linear_built(false), _fd_vec_base(0),
      _time_linbuild(0.0), _time_nonlin(0.0)
{

}


EpotMatrixSolver::EpotMatrixSolver( Geometry &geom, std::istream &s )
    : EpotSolver(geom,s)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotMatrixSolver::~EpotMatrixSolver()
{
    if( _fd_mat )
        delete _fd_mat;
    if( _fd_vec )
        delete _fd_vec;
    if( _fd_vec_base )
        delete _fd_vec_base;
}


/* Set a link in matrix system.
 *
 * Makes the node a dependent on node b. If b is positive this means
 * that element (a,b) of matrix A is set to value val. If b is
 * negative, the node has a fixed potential and the dependence should
 * be added to the vector side of the system of equations on row
 * a. The value of potential at node b is stored in epot(-b).
 */
void EpotMatrixSolver::set_link( uint32_t a, uint32_t b, double val )
{
    //std::cout << "set_link( a = " << a << ", b = " << b << ", val = " << val << ")\n";

    if( (b & N2D_TYPE_MASK) == N2D_TYPE_FIXED ) {

	/* A link pointing INTO a Neumann-mask solid is a symmetry face:
	 * redirect it onto the MIRROR-IMAGE node instead of moving it to the
	 * right hand side.
	 *
	 * This matches what ibsimu already does on a Neumann BOX wall, which
	 * is the convention the rest of the code is built around and which
	 * should have been the starting point. add_neumann_node_2d() emits
	 *     set_link( a, _n2d(i,j-1), 2.0*eps_self );   cof += 2.0*eps_self
	 * -- weight 2*eps on the single INTERIOR neighbour. That is the
	 * standard five-point stencil with the ghost eliminated by
	 *     phi_(j+1) := phi_(j-1),
	 * i.e. a NODE-CENTRED mirror whose zero-derivative plane lies exactly
	 * ON the boundary node. The factor of two on one neighbour is its
	 * signature.
	 *
	 * The mirror image of neighbour b about node a is 2a - b in linear
	 * index space, for any axis and any geometry mode, because the mesh
	 * index is affine: a neighbour is a +-1 or +-nx or +-nx*ny offset, so
	 * negating the offset negates the index difference. Redirecting the
	 * face weight there gives, for a node with mask on +y,
	 *     w_m*phi_(j-1) + w_p*phi_(j-1) - (w_m+w_p)*phi_j
	 * which for uniform eps is exactly 2*eps*phi_(j-1) - 2*eps*phi_j --
	 * the box-wall stencil above. The duplicate-column merge in
	 * build_mat_vec() folds the two entries on column (j-1) into one.
	 *
	 * The PREVIOUS version of this sent the weight to the row's own
	 * diagonal, which imposes phi_ghost = phi_self. That is also a valid
	 * zero-flux condition, but a FACE-centred one: its symmetry plane
	 * sits half a cell OUTSIDE the last live node. Particles reflect off
	 * the analytic solid surface at the node, so field and particles were
	 * mirroring half a cell apart, leaving a half-cell the field treated
	 * as symmetric but no particle ever entered. It also did not match
	 * the *= 2.0 space-charge correction that scharge_finalize_*() applies
	 * at box walls, which exists precisely because a node-centred mirror
	 * node owns half a control volume.
	 *
	 * Doing it in set_link() still means it applies to every stencil
	 * builder at once -- add_vacuum_node(), add_neumann_node_1d/2d/cyl/3d(),
	 * add_near_solid_node_1d/2d/cyl/3d() -- with no per-site changes.
	 *
	 * The masked node's own index is a MESH node index here: preprocess()
	 * sets _n2d(a) = N2D_TYPE_FIXED | a for every fixed node, so the
	 * index survives unchanged and _geom.mesh() can be queried directly.
	 */
	if( SMESH_NODE_IS_NEUMANN_MASK( _geom.mesh( (int32_t)(b & N2D_INDEX_MASK) ) ) ) {

	    /* Mirror in MESH index space. The row index a is a DOF index and
	     * b (being FIXED) carries a MESH index, so 2*a - b would subtract
	     * one index space from the other and land on an arbitrary column.
	     * _dof2node maps the row back to its mesh node first. */
	    const int32_t dof_a  = (int32_t)(a & N2D_INDEX_MASK);
	    const int32_t node_a = _dof2node[dof_a];
	    const int32_t node_b = (int32_t)(b & N2D_INDEX_MASK);
	    const int32_t node_m = 2*node_a - node_b;   /* image of b about a */

	    if( node_m >= 0 && node_m < (int32_t)_geom.nodecount() ) {

		uint32_t n2d_m = _n2d(node_m);

		if( (n2d_m & N2D_TYPE_MASK) != N2D_TYPE_FIXED ) {
		    /* Image is a free node: it carries a DOF index, which is
		     * what a matrix column must be. */
		    _row_entries[dof_a].push_back(
			std::make_pair( (int32_t)(n2d_m & N2D_INDEX_MASK), val ) );
		    return;
		} else if( !SMESH_NODE_IS_NEUMANN_MASK( _geom.mesh( node_m ) ) ) {
		    /* Image is a Dirichlet node: its potential is known, so
		     * the term belongs on the right hand side, exactly as an
		     * ordinary fixed neighbour would. */
		    (*_fd_vec)(a) += -val * (*_epot)( node_m );
		    return;
		}
		/* Image is itself masked -- fall through. */
	    }

	    /* No usable image: off the mesh, or masked on both sides (a live
	     * region one node wide). Fold onto the diagonal instead, which is
	     * the face-centred mirror phi_ghost = phi_self. Less accurate,
	     * but local and diagonally dominant rather than referencing a
	     * node whose value means nothing. */
	    _row_entries[dof_a].push_back( std::make_pair( dof_a, val ) );
	    return;
	}

	//std::cout << "epot = " << (*_epot)(b & N2D_INDEX_MASK) << "\n";
        (*_fd_vec)(a) += -val * (*_epot)(b & N2D_INDEX_MASK);
    } else {
	//std::cout << "matrix construct\n";
	// Buffered instead of a direct CRowMatrix::construct_add() call
	// so build_mat_vec() can fill different rows from different
	// threads -- see _row_entries' doc comment in the header. Row
	// a's buffer is only ever touched by the single thread that
	// owns node/row a in build_mat_vec()'s parallel loop.
	_row_entries[a & N2D_INDEX_MASK].push_back(
	    std::make_pair( (int32_t)(b & N2D_INDEX_MASK), val ) );
    }
}


/* Epot-dependent plasma contribution to rhs/diagonal for free node a.
 * Consolidated out of add_vacuum_node()/add_near_solid_node()/
 * add_neumann_node(), which used to each carry an identical copy of
 * this block. Called once per Newton iteration (not cached), since it
 * is the one part of matrix/rhs construction that genuinely depends on
 * the current solution guess X via _sol.
 */
void EpotMatrixSolver::update_nonlinear_node( uint32_t a, uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    const bool have_plasma = ( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP ||
                               _plasma == PLASMA_SHIELD );
    if( !have_plasma && !_comp )
	return;

    // A dielectric solid's interior is a fixed, chargeless medium --
    // never actual plasma -- regardless of how close it sits to a real
    // plasma region. This function is called for every free node
    // (vacuum, near-solid, and dielectric interior alike, since a
    // dielectric's own row is otherwise an ordinary free node -- see
    // build_mesh_parallel_thread_3d()), so without this check the
    // nonlinear electron/Boltzmann source term below gets added to a
    // dielectric's own (correctly eps-weighted, charge-free) row too,
    // pulling its potential toward the plasma equilibrium instead of
    // solving the linear Laplace equation it should. This also catches
    // a dielectric node reclassified NEAR_SOLID/NEUMANN at a box edge
    // or next to a conductor, where the raw tag alone no longer says
    // "dielectric".
    //
    // Geometry::dielectric_material_at() is an O(1) read (see its doc
    // comment) that's already correct regardless of any such
    // reclassification, computed once for free during mesh build --
    // this function runs every Newton iteration (and every step-size
    // backtracking re-evaluation) for every free node, so it matters
    // that this isn't a Geometry::inside() query repeated on every call.
    const bool dielectric_free = ( _geom.dielectric_material_at( i, j, k ) == 0 );

    // The two models are ACCUMULATED, not selected between. A source plasma
    // and a compensated drift coexist in one mesh and a run needs both; they
    // are gated by SEPARATE region functions because they occupy different
    // parts of it.
    double fd = 0.0, d = 0.0;

    if( dielectric_free && have_plasma ) {
	bool inplasma = true;
	if( _plasma_calc_func )
	    inplasma = (*_plasma_calc_func)(x);
	if( inplasma ) {
	    double p = (*_sol)(a);
	    double rhst, drhst;
	    if( _plasma == PLASMA_PEXP ) {
		pexp_newton( rhst, drhst, p );
		fd += rhst;
		d  += drhst;
	    } else if( _plasma == PLASMA_NSIMP ) {
		nsimp_newton( rhst, drhst, p );
		fd += rhst;
		d  += drhst;
	    } else {
		shield_newton( rhst, drhst, p );
		double R = -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
		fd += R*rhst;
		d  += R*drhst;
	    }
	}
    }

    if( dielectric_free && _comp ) {
	bool incomp = true;
	if( _comp_region_func )
	    incomp = (*_comp_region_func)(x);
	if( incomp ) {
	    // Prefactor is the LOCAL ion density, exactly as the shield model
	    // takes its own -- so _scharge must hold ION charge only, with the
	    // compensating charge generated here rather than supplied.
	    double rho_i = (*_scharge)(i,j,k);
	    if( rho_i > _compRhoMin ) {
		double p = (*_sol)(a);
		double g, dg;
		comp_newton( g, dg, p );
		// rho_e = -rho_i*g contributes -rho_e*h^2/eps0 = +rho_i*g*h^2/eps0,
		// the same sign convention the other models use.
		double R = rho_i*_geom.h()*_geom.h()/EPSILON0;
		fd += R*g;
		d  += R*dg;
	    }
	}
    }

    (*_fd_vec)(a) += fd;
    _d_vec(a) = d;
}


/* See the doc comment on node_epsilon_r() in the header. */
double EpotMatrixSolver::node_epsilon_r( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint32_t solid_number = _geom.dielectric_material_at( i, j, k );
    if( solid_number == 0 )
	return( 1.0 ); // plain vacuum
    // boundary(), not get_boundary(): solid_number came straight from
    // dielectric_material_at() so it is valid by construction, and this
    // runs per node -- no need to copy a Bound or re-range-check it.
    return( _geom.boundary( solid_number ).value() );
}


double EpotMatrixSolver::face_epsilon_alpha( double eps_self, double eps_neighbor, double alpha )
{
    return( 1.0/( alpha/eps_self + (1.0-alpha)/eps_neighbor ) );
}


/* See the doc comment on vacuum_face_coefficient() in the header. */
double EpotMatrixSolver::vacuum_face_coefficient( int32_t i, int32_t j, int32_t k,
						   int32_t ni, int32_t nj, int32_t nk,
						   int sign, int coord,
						   uint32_t neighbor_mesh, uint32_t self_material,
						   double eps_self ) const
{
    uint32_t node_id = neighbor_mesh & SMESH_NODE_ID_MASK;

    /* Neumann mask on the far side: this face carries no flux at all, so
     * its coefficient is about to be cancelled by set_link()'s redirect
     * (it sends the value to the row's own diagonal, where it exactly
     * undoes the caller's cof += w). The returned number is therefore
     * arbitrary -- but computing it properly is not merely wasted work,
     * it is meaningless work: the code below would see nb_material == 0
     * (dielectric_material_at() returns 0 for a mask, since a mask is not
     * BOUND_DIELECTRIC), decide there is a genuine material interface, and
     * bisect against the dielectric with bracket_ndist() from a point that
     * may well be INSIDE it -- which violates that routine's stated
     * precondition and returns a meaningless fraction.
     *
     * eps_self keeps it finite and dimensionally sane for anyone reading a
     * matrix dump. It cannot affect the solution.
     */
    if( SMESH_NODE_IS_NEUMANN_MASK( neighbor_mesh ) )
	return( eps_self );

    if( node_id == SMESH_NODE_ID_DIRICHLET ) {
	uint32_t boundary_number = neighbor_mesh & SMESH_BOUNDARY_NUMBER_MASK;
	if( boundary_number < 7 )
	    // Simulation box edge, not a user-defined solid -- its
	    // position is exactly the mesh's own coordinate system (no
	    // Solid::inside() to bisect against, and no ambiguity: a
	    // dielectric node next to it is, by construction, always
	    // exactly one full cell away -- see
	    // Geometry::override_dielectric_box_boundary_3d()). Assume
	    // self's own medium extends right up to it, as before.
	    return( eps_self );

	// A real solid electrode (n>=7): unlike the box edge, this can
	// sit at any sub-cell distance from self, e.g. an STL electrode
	// touching an STL dielectric. Bisect against that solid's own
	// geometry to find the true distance, rather than assuming
	// flush contact (alpha=1, eps_self unweighted, as an earlier
	// version of this function did). self (dielectric or vacuum) is
	// always geometrically outside the conductor solid, so -- unlike
	// the dielectric-vs-dielectric case below -- this never needs
	// the reversed/complemented bisection.
	//
	// Physically: beyond the conductor surface is a known constant
	// (its fixed voltage), not another variable material, so the
	// flux from self to it is a single resistor of length
	// alpha*h through self's own medium: eps_self/alpha. This
	// reduces to the old eps_self exactly at alpha=1 (conductor a
	// full cell away), and grows without bound as alpha->0 (the
	// conductor sitting arbitrarily close), consistent with how the
	// existing near-solid conductor formulas also blow up in that
	// limit.
	double alpha = _geom.solid_face_frac( i, j, k, boundary_number, sign, coord );
	return( eps_self / alpha );
    }

    // The neighbour's raw tag alone is not enough here: a Neumann (or
    // near-solid/fine-boundary) neighbour's own tag may itself have
    // been overridden away from a dielectric tag that reached that
    // node -- e.g. a Neumann side wall at a z-level that a dielectric
    // slab varying along z actually occupies there.
    // Geometry::dielectric_material_at() recovers the true answer
    // directly (O(1), independent of any such reclassification -- see
    // its doc comment), so no on-demand Geometry::inside() fallback or
    // solver-private cache is needed here at all.
    uint32_t nb_material = _geom.dielectric_material_at( ni, nj, nk );
    if( nb_material == self_material )
	return( eps_self ); // same medium on both sides of this face -- no correction needed

    // Genuine material interface on this face -- find its true sub-cell
    // position by bisecting against whichever side is an actual
    // dielectric solid, rather than assuming it sits exactly halfway
    // (only true for a grid-aligned interface).
    //
    // Geometry::solid_face_frac() (bracket_ndist()) bisects assuming
    // its (i,j,k) argument is OUTSIDE the target solid, walking towards
    // it -- true whenever the *neighbour* is the dielectric (self is
    // plain vacuum here). But when *self* is the dielectric being
    // bisected against (self_material == bisect_solid), self is
    // INSIDE that solid, not outside, and calling it from self's side
    // would search assuming the wrong monotonicity and return
    // garbage (caught by a non-grid-aligned analytic verification
    // test, not by inspection) -- so in that case bisect from the
    // neighbour (guaranteed to be on the outside, since its material
    // differs) walking back with the opposite sign instead, and take
    // the complementary fraction.
    uint32_t bisect_solid = ( self_material != 0 ) ? self_material : nb_material;
    double eps_neighbor = ( nb_material != 0 ) ? _geom.boundary( nb_material ).value() : 1.0;
    double alpha;
    if( bisect_solid == self_material )
	alpha = 1.0 - _geom.solid_face_frac( ni, nj, nk, bisect_solid, -sign, coord );
    else
	alpha = _geom.solid_face_frac( i, j, k, bisect_solid, sign, coord );
    return( face_epsilon_alpha( eps_self, eps_neighbor, alpha ) );
}


void EpotMatrixSolver::add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    (void)x; // no longer used here -- see update_nonlinear_node()
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    // Relative permittivity of this node, and material number (0 =
    // plain vacuum, else the dielectric solid this node is inside).
    // Per direction, vacuum_face_coefficient() compares this against
    // the neighbour's own material: matching materials (including
    // vacuum-vacuum, the common case) reduce to a plain eps_self link,
    // exactly reproducing the original fixed 1.0/-2.0/-4.0/-6.0 stencil
    // whenever no dielectric is involved anywhere near this node. A
    // genuine material interface on a face gets a fractional-distance-
    // weighted coefficient found by bisecting against the true (e.g.
    // STL-imported) solid surface -- see vacuum_face_coefficient() --
    // so, unlike the initial grid-aligned-only implementation, this is
    // correct for an arbitrarily positioned dielectric-vacuum
    // interface. A dielectric next to a Dirichlet conductor (fixed
    // node) at an arbitrary sub-cell position is not yet corrected this
    // way -- see vacuum_face_coefficient()'s doc comment for that
    // scoped limitation.
    uint32_t self_material = _geom.dielectric_material_at( i, j, k );
    double eps_self = node_epsilon_r( i, j, k );
    double cof = 0.0;

    switch( _geom.geom_mode() ) {
    case MODE_1D: {
        double wm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), self_material, eps_self );
        double wp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), self_material, eps_self );
        set_link( a, _n2d(i-1,j,k), wm );
        set_link( a, _n2d(i+1,j,k), wp );
        cof = wm+wp;
        break; }
    case MODE_2D: {
        double wxm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), self_material, eps_self );
        double wxp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), self_material, eps_self );
        double wym = vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), self_material, eps_self );
        double wyp = vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), self_material, eps_self );
        set_link( a, _n2d(i,j-1,k), wym );
        set_link( a, _n2d(i-1,j,k), wxm );
        set_link( a, _n2d(i+1,j,k), wxp );
        set_link( a, _n2d(i,j+1,k), wyp );
        cof = wxm+wxp+wym+wyp;
        break; }
    case MODE_CYL: {
        double wxm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), self_material, eps_self );
        double wxp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), self_material, eps_self );
        double wym = (1.0-0.5/j)*vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), self_material, eps_self );
        double wyp = (1.0+0.5/j)*vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), self_material, eps_self );
        set_link( a, _n2d(i,j-1,k), wym );
        set_link( a, _n2d(i-1,j,k), wxm );
        set_link( a, _n2d(i+1,j,k), wxp );
        set_link( a, _n2d(i,j+1,k), wyp );
        cof = wxm+wxp+wym+wyp;
        break; }
    case MODE_3D: {
        double wxm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), self_material, eps_self );
        double wxp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), self_material, eps_self );
        double wym = vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), self_material, eps_self );
        double wyp = vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), self_material, eps_self );
        double wzm = vacuum_face_coefficient( i,j,k, i,j,k-1, -1,2, _geom.mesh(i,j,k-1), self_material, eps_self );
        double wzp = vacuum_face_coefficient( i,j,k, i,j,k+1, +1,2, _geom.mesh(i,j,k+1), self_material, eps_self );
        set_link( a, _n2d(i,j,k-1), wzm );
        set_link( a, _n2d(i,j-1,k), wym );
        set_link( a, _n2d(i-1,j,k), wxm );
        set_link( a, _n2d(i+1,j,k), wxp );
        set_link( a, _n2d(i,j+1,k), wyp );
        set_link( a, _n2d(i,j,k+1), wzp );
        cof = wxm+wxp+wym+wyp+wzm+wzp;
        break; }
    }

    set_link( a, _n2d(i,j,k), -cof );

    // Nonlinear (plasma) contribution to rhs/diagonal is epot-dependent
    // and is added separately, once per Newton iteration, by
    // update_nonlinear_node() -- see build_mat_vec(). Everything in
    // this function is geometry-only (now including permittivity,
    // which is a fixed material property) and gets cached after the
    // first build.

    // Right hand side. Left as-is (no epsilon_r factor) even inside a
    // dielectric: this is the flux-conservative form of
    // div(eps*grad(phi)) = -rho/eps0, so the permittivity belongs
    // entirely on the coefficients above, not the rhs.
    //
    // self_material == 0 guards this rather than trusting scharge to
    // already be zero inside a dielectric's bulk: PIC/linear deposition
    // has no notion of Geometry at all (see scharge_clear_solid_nodes()'s
    // doc comment in scharge.hpp), so without this guard, correctness
    // here would depend on every scharge-producing code path upstream
    // remembering to clean solid nodes out first -- exactly the kind of
    // cross-module invariant that's easy to silently break later. This
    // is a genuine dielectric interior node (this function is shared
    // between plain vacuum and dielectric interior -- see the dielectric
    // material/permittivity comment above), so it truly has no free
    // charge in its own bulk equation, unlike add_near_solid_node()'s and
    // add_neumann_node()'s unconditional reads of the same scharge term,
    // which are genuine vacuum locations that can legitimately have real
    // nearby charge and must keep reading it.
    if( _plasma != PLASMA_SHIELD && self_material == 0 )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::add_near_solid_node_1d( uint32_t i, const Vec3D &x )
{
    uint32_t a = _n2d(i) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis. See add_near_solid_node_3d() for the
    // rationale of the conductor/dielectric same-axis triple-point
    // handling -- unlike the 2d/cyl/3d variants, this function has no
    // separate "neither side near a conductor" branch (impossible here:
    // being in add_near_solid_node_1d at all means is_near_solid() was
    // true, and X is the only axis, so sflag always has at least one
    // bit set), so the triple-point check applies directly in the
    // final else below.
    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i), -2.0/(beta*beta) );
	set_link( a, _n2d(i+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a, _n2d(i), -2.0/(alpha*alpha) );
	set_link( a, _n2d(i-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	bool xm_conductor = sflag & 0x01;
	bool xp_conductor = sflag & 0x02;
	bool xm_dielectric = !xm_conductor && _geom.dielectric_material_at(i-1,0,0) != 0;
	bool xp_dielectric = !xp_conductor && _geom.dielectric_material_at(i+1,0,0) != 0;

	if( xm_dielectric || xp_dielectric ) {
	    double wm = xm_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,0,0, i-1,0,0, -1,0, _geom.mesh(i-1), 0, 1.0 );
	    double wp = xp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,0,0, i+1,0,0, +1,0, _geom.mesh(i+1), 0, 1.0 );
	    set_link( a, _n2d(i), -(wm+wp) );
	    set_link( a, _n2d(i-1), wm );
	    set_link( a, _n2d(i+1), wp );
	} else {
	    set_link( a, _n2d(i), -2.0/(alpha*beta) );
	    set_link( a, _n2d(i-1), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i+1), 2.0/((alpha+beta)*beta) );
	}
    }
}


void EpotMatrixSolver::add_near_solid_node_2d( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis. See add_near_solid_node_3d() for the
    // rationale of the sflag-gated third branch, including the
    // conductor/dielectric same-axis triple-point handling.
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else if( sflag & 0x03 ) {
	bool xm_conductor = sflag & 0x01;
	bool xp_conductor = sflag & 0x02;
	bool xm_dielectric = !xm_conductor && _geom.dielectric_material_at(i-1,j,0) != 0;
	bool xp_dielectric = !xp_conductor && _geom.dielectric_material_at(i+1,j,0) != 0;

	if( xm_dielectric || xp_dielectric ) {
	    double wm = xm_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), 0, 1.0 );
	    double wp = xp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), 0, 1.0 );
	    set_link( a, _n2d(i-1,j), wm );
	    set_link( a, _n2d(i+1,j), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
	}
    } else {
	double wm = vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), 0, 1.0 );
	set_link( a, _n2d(i-1,j), wm );
	set_link( a, _n2d(i+1,j), wp );
	cof += wm+wp;
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else if( sflag & 0x0c ) {
	bool ym_conductor = sflag & 0x04;
	bool yp_conductor = sflag & 0x08;
	bool ym_dielectric = !ym_conductor && _geom.dielectric_material_at(i,j-1,0) != 0;
	bool yp_dielectric = !yp_conductor && _geom.dielectric_material_at(i,j+1,0) != 0;

	if( ym_dielectric || yp_dielectric ) {
	    double wm = ym_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,0, i,j-1,0, -1,1, _geom.mesh(i,j-1), 0, 1.0 );
	    double wp = yp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,0, i,j+1,0, +1,1, _geom.mesh(i,j+1), 0, 1.0 );
	    set_link( a, _n2d(i,j-1), wm );
	    set_link( a, _n2d(i,j+1), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i,j-1), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i,j+1), 2.0/((alpha+beta)*beta) );
	}
    } else {
	double wm = vacuum_face_coefficient( i,j,0, i,j-1,0, -1,1, _geom.mesh(i,j-1), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,0, i,j+1,0, +1,1, _geom.mesh(i,j+1), 0, 1.0 );
	set_link( a, _n2d(i,j-1), wm );
	set_link( a, _n2d(i,j+1), wp );
	cof += wm+wp;
    }

    // Middle node
    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_near_solid_node_cyl( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis (axial). See add_near_solid_node_3d() for the
    // rationale of the sflag-gated third branch, including the
    // conductor/dielectric same-axis triple-point handling. (Y, the
    // radial axis below, is a separate, still-unresolved limitation --
    // see its own comment.)
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else if( sflag & 0x03 ) {
	bool xm_conductor = sflag & 0x01;
	bool xp_conductor = sflag & 0x02;
	bool xm_dielectric = !xm_conductor && _geom.dielectric_material_at(i-1,j,0) != 0;
	bool xp_dielectric = !xp_conductor && _geom.dielectric_material_at(i+1,j,0) != 0;

	if( xm_dielectric || xp_dielectric ) {
	    double wm = xm_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), 0, 1.0 );
	    double wp = xp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), 0, 1.0 );
	    set_link( a, _n2d(i-1,j), wm );
	    set_link( a, _n2d(i+1,j), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
	}
    } else {
	double wm = vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), 0, 1.0 );
	set_link( a, _n2d(i-1,j), wm );
	set_link( a, _n2d(i+1,j), wp );
	cof += wm+wp;
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis (radial). The cylindrical curvature terms are
    // baked into these coefficients, so unlike the other axes they are not
    // a flat Shortley-Weller form and the dielectric-aware version had to
    // be derived rather than copied.
    //
    // DERIVATION (conservative / finite-volume, which is what a material
    // interface requires -- flux must be continuous across it):
    //
    //   Control volume for node j spans the inner face at
    //   r_- = r_j - alpha*h/2 and the outer face at r_+ = r_j + beta*h/2.
    //   The radial flux through a face is eps_face * r_face * dphi/dr, so
    //   in this file's row scaling (a uniform vacuum node has cof = 4 and
    //   rhs = -rho h^2 / eps0) each face coefficient is
    //
    //       w = eps_face * (r_face / r_j) / (face distance in units of h)
    //
    //   giving  w_m = eps_m * (j - alpha/2) / (j*alpha)
    //           w_p = eps_p * (j + beta /2) / (j*beta )
    //
    // Checked in three limits:
    //   alpha=beta=1  -> eps_m*(1-1/(2j)), eps_p*(1+1/(2j)), i.e. exactly
    //                    add_vacuum_node()'s cylindrical form;
    //   j -> infinity -> eps_m/alpha, eps_p/beta, i.e. exactly the axial
    //                    dielectric branch above;
    //   eps = 1       -> (alpha+beta)/2 times the Taylor form kept below.
    //
    // That last ratio is not an error: it is the same relationship the
    // AXIAL axis already has between its plain branch (2/((alpha+beta)*alpha),
    // a Taylor expansion, more accurate for smooth coefficients) and its
    // dielectric branch (1/alpha, conservative, exact in flux). They agree
    // for uniform spacing and differ only in cut cells, and the
    // conservative one is the correct choice when eps jumps.
    //
    // Validated against the analytic coaxial two-layer dielectric solution
    // phi = C - (A/eps) ln r: this form converges at second order, whereas
    // ignoring eps on the radial faces leaves a ~24% error that does NOT
    // shrink with the mesh -- a wrong equation rather than a resolution
    // problem.
    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis. The factor 4 is the coordinate-singularity limit of the
	// radial operator, not an ordinary face, so -- as in
	// add_neumann_node_cyl() -- no cut-cell or interface treatment is
	// attempted for a dielectric sitting exactly on the axis.
	cof += 4.0;
	set_link( a, _n2d(i,j+1), 4.0 );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	// Node on the rmax face. Left as the vacuum form: reaching here
	// needs a node that is simultaneously on the outer box face and
	// within a cell of a conductor, and the box face is a boundary
	// condition rather than a material interface.
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += (1.0)/(2.0*alpha)*(2.0/alpha+1.0/j)*
	    2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else {
	bool ym_conductor = sflag & 0x04;
	bool yp_conductor = sflag & 0x08;
	bool ym_dielectric = !ym_conductor && _geom.dielectric_material_at(i,j-1,0) != 0;
	bool yp_dielectric = !yp_conductor && _geom.dielectric_material_at(i,j+1,0) != 0;

	if( ym_dielectric || yp_dielectric ) {

	    // A near-solid node is by definition a vacuum node next to a
	    // conductor, so self is vacuum (self_material 0, eps_self 1)
	    // and a conductor face contributes eps = 1 with its distance
	    // already carried by alpha/beta -- passing it through
	    // vacuum_face_coefficient() would apply that distance twice.
	    double rj = (double)j;
	    double em = ym_conductor ? 1.0 :
		vacuum_face_coefficient( i,j,0, i,j-1,0, -1,1, _geom.mesh(i,j-1), 0, 1.0 );
	    double ep = yp_conductor ? 1.0 :
		vacuum_face_coefficient( i,j,0, i,j+1,0, +1,1, _geom.mesh(i,j+1), 0, 1.0 );

	    double wm = em*(rj - 0.5*alpha)/(rj*alpha);
	    double wp = ep*(rj + 0.5*beta )/(rj*beta );

	    set_link( a, _n2d(i,j-1), wm );
	    set_link( a, _n2d(i,j+1), wp );
	    cof += wm+wp;

	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i,j-1), 1.0/(alpha+beta)*(2.0/alpha-1.0/j) );
	    set_link( a, _n2d(i,j+1), 1.0/(alpha+beta)*(2.0/beta+1.0/j) );
	}
    }

    // Middle node
    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_near_solid_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j,k);

    const uint8_t *nearsolid_ptr = _geom.nearsolid_ptr( _geom.mesh(i,j,k) & SMESH_NEAR_SOLID_INDEX_MASK );
    uint8_t sflag = nearsolid_ptr[0];
    uint8_t *ptr = (uint8_t *)&nearsolid_ptr[1];

    double cof = 0.0;

    // Xmin direction
    double alpha = 1.0;
    if( sflag & 0x01 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Xmax direction
    double beta = 1.0;
    if( sflag & 0x02 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j,k), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j,k), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else if( sflag & 0x03 ) {
	// Near a conductor on (at least) one side of this axis. If a
	// conductor is cached on *both* sides, self's own medium (vacuum,
	// eps=1) genuinely spans between the two conductors with no
	// material change anywhere on this axis, so the existing
	// Taylor-difference Shortley-Weller treatment applies exactly.
	// But is_solid() (and therefore this near-solid cache) never
	// recognizes a dielectric -- only a conductor -- so if just *one*
	// side is cached, the *other* side (defaulted to alpha/beta=1.0
	// above) might actually be a dielectric-interior cell rather than
	// plain vacuum: a true conductor/dielectric/vacuum triple point.
	// The joint Taylor formula has no way to fold in a permittivity
	// discontinuity (it assumes a well-defined second derivative
	// through both neighbours; a material interface only guarantees
	// flux continuity, i.e. a kink, not that). Detect that case via
	// _dielectric_mat and, only then, recompose this axis from two
	// independent per-face resistor conductances instead: the
	// conductor side reuses its already-bisected cached fractional
	// distance (eps_self/alpha or eps_self/beta, eps_self=1 here --
	// exactly vacuum_face_coefficient()'s own Dirichlet-branch
	// formula, just without re-bisecting), and the other side goes
	// through vacuum_face_coefficient() itself, which bisects/weights
	// correctly for the dielectric neighbour. When neither uncached
	// side turns out to be a dielectric, this reduces to the ordinary
	// case below and the existing Taylor formula is used unchanged.
	bool xm_conductor = sflag & 0x01;
	bool xp_conductor = sflag & 0x02;
	bool xm_dielectric = !xm_conductor && _geom.dielectric_material_at(i-1,j,k) != 0;
	bool xp_dielectric = !xp_conductor && _geom.dielectric_material_at(i+1,j,k) != 0;

	if( xm_dielectric || xp_dielectric ) {
	    double wm = xm_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), 0, 1.0 );
	    double wp = xp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), 0, 1.0 );
	    set_link( a, _n2d(i-1,j,k), wm );
	    set_link( a, _n2d(i+1,j,k), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i-1,j,k), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i+1,j,k), 2.0/((alpha+beta)*beta) );
	}
    } else {
	// Neither side of this axis is near a conductor (alpha=beta=1,
	// both defaulted, no _nearsolid data cached for either) -- this
	// near-solid node (near a conductor on a *different* axis) can
	// still be directly adjacent to a dielectric along *this* axis,
	// so use the same per-face dielectric-aware treatment as
	// add_vacuum_node()/add_neumann_node() rather than assuming
	// plain vacuum unconditionally. self_material=0, eps_self=1.0
	// always here (this function never runs for a dielectric-
	// interior self node).
	double wm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), 0, 1.0 );
	set_link( a, _n2d(i-1,j,k), wm );
	set_link( a, _n2d(i+1,j,k), wp );
	cof += wm+wp;
    }

    // Ymin direction
    alpha = 1.0;
    if( sflag & 0x04 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Ymax direction
    beta = 1.0;
    if( sflag & 0x08 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Y axis
    if( bindex & EPOT_SOLVER_BYMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j+1,k), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1,k), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else if( sflag & 0x0c ) {
	// See the X-axis branch above for the rationale (conductor and
	// dielectric on opposite sides of the *same* axis).
	bool ym_conductor = sflag & 0x04;
	bool yp_conductor = sflag & 0x08;
	bool ym_dielectric = !ym_conductor && _geom.dielectric_material_at(i,j-1,k) != 0;
	bool yp_dielectric = !yp_conductor && _geom.dielectric_material_at(i,j+1,k) != 0;

	if( ym_dielectric || yp_dielectric ) {
	    double wm = ym_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), 0, 1.0 );
	    double wp = yp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), 0, 1.0 );
	    set_link( a, _n2d(i,j-1,k), wm );
	    set_link( a, _n2d(i,j+1,k), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i,j-1,k), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i,j+1,k), 2.0/((alpha+beta)*beta) );
	}
    } else {
	double wm = vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), 0, 1.0 );
	set_link( a, _n2d(i,j-1,k), wm );
	set_link( a, _n2d(i,j+1,k), wp );
	cof += wm+wp;
    }

    // Zmin direction
    alpha = 1.0;
    if( sflag & 0x10 ) {
	alpha = *ptr/255.0;
	ptr++;
    }

    // Zmax direction
    beta = 1.0;
    if( sflag & 0x20 ) {
	beta = *ptr/255.0;
	ptr++;
    }

    // Factors for Z axis
    if( bindex & EPOT_SOLVER_BZMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i,j,k+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(5).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j,k-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(6).value(x) / alpha;
    } else if( sflag & 0x30 ) {
	// See the X-axis branch above.
	bool zm_conductor = sflag & 0x10;
	bool zp_conductor = sflag & 0x20;
	bool zm_dielectric = !zm_conductor && _geom.dielectric_material_at(i,j,k-1) != 0;
	bool zp_dielectric = !zp_conductor && _geom.dielectric_material_at(i,j,k+1) != 0;

	if( zm_dielectric || zp_dielectric ) {
	    double wm = zm_conductor ? 1.0/alpha :
		vacuum_face_coefficient( i,j,k, i,j,k-1, -1,2, _geom.mesh(i,j,k-1), 0, 1.0 );
	    double wp = zp_conductor ? 1.0/beta :
		vacuum_face_coefficient( i,j,k, i,j,k+1, +1,2, _geom.mesh(i,j,k+1), 0, 1.0 );
	    set_link( a, _n2d(i,j,k-1), wm );
	    set_link( a, _n2d(i,j,k+1), wp );
	    cof += wm+wp;
	} else {
	    cof += 2.0/(alpha*beta);
	    set_link( a, _n2d(i,j,k-1), 2.0/((alpha+beta)*alpha) );
	    set_link( a, _n2d(i,j,k+1), 2.0/((alpha+beta)*beta) );
	}
    } else {
	double wm = vacuum_face_coefficient( i,j,k, i,j,k-1, -1,2, _geom.mesh(i,j,k-1), 0, 1.0 );
	double wp = vacuum_face_coefficient( i,j,k, i,j,k+1, +1,2, _geom.mesh(i,j,k+1), 0, 1.0 );
	set_link( a, _n2d(i,j,k-1), wm );
	set_link( a, _n2d(i,j,k+1), wp );
	cof += wm+wp;
    }

    // Middle node
    set_link( a, _n2d(i,j,k), -cof );
}


void EpotMatrixSolver::add_near_solid_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	add_near_solid_node_1d(i,x);
	break;
    case MODE_2D:
	add_near_solid_node_2d(i,j,x);
	break;
    case MODE_CYL:
	add_near_solid_node_cyl(i,j,x);
	break;
    case MODE_3D:
	add_near_solid_node_3d(i,j,k,x);
	break;
    }

    // Nonlinear (plasma) contribution to rhs/diagonal is epot-dependent
    // and is added separately, once per Newton iteration, by
    // update_nonlinear_node() -- see build_mat_vec(). Everything in
    // this function is geometry-only and gets cached after the first
    // build.

    // Right hand side
    if( _plasma != PLASMA_SHIELD )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::add_neumann_node_1d( uint32_t i, const Vec3D &x )
{
    uint32_t a = _n2d(i) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i);

    // A Neumann box-wall node's own mesh tag generally cannot carry
    // material info (it may have been overridden away from a dielectric
    // tag that reached this edge), so its permittivity has to be
    // recovered from Geometry::dielectric_material_at() (O(1), correct
    // regardless of that reclassification) rather than the raw tag
    // alone. This matters whenever some *other* axis has a material
    // boundary near this node -- in 1D there is only one axis, so this
    // mainly keeps add_neumann_node_1d() consistent with the other geom
    // modes.
    uint32_t self_material = _geom.dielectric_material_at( i, 0, 0 );
    double eps_self = ( self_material == 0 ) ? 1.0 : _geom.get_boundary( self_material ).value();

    if( bindex & EPOT_SOLVER_BXMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i), -w );
	set_link( a, _n2d(i+1), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(1).value( x );
    } else {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i-1), w );
	set_link( a, _n2d(i), -w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(2).value( x );
    }
}


void EpotMatrixSolver::add_neumann_node_2d( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    // See add_neumann_node_1d() -- recover this node's true permittivity
    // even though its own tag may not carry it, since a material
    // boundary along the *other* axis can sit right next to a
    // Neumann-tagged node on this one.
    uint32_t self_material = _geom.dielectric_material_at( i, j, 0 );
    double eps_self = ( self_material == 0 ) ? 1.0 : _geom.get_boundary( self_material ).value();
    double cof = 0.0;

    if( bindex & EPOT_SOLVER_BXMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i+1,j), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(1).value( x );
	cof += w;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i-1,j), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(2).value( x );
	cof += w;
    } else {
	double wm = vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), self_material, eps_self );
	double wp = vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), self_material, eps_self );
	set_link( a, _n2d(i-1,j), wm );
	set_link( a, _n2d(i+1,j), wp );
	cof += wm+wp;
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j+1), w );
	(*_fd_vec)(a) += -w*_geom.h()*_geom.get_boundary(3).value( x );
	cof += w;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j-1), w );
	(*_fd_vec)(a) += -w*_geom.h()*_geom.get_boundary(4).value( x );
	cof += w;
    } else {
	double wm = vacuum_face_coefficient( i,j,0, i,j-1,0, -1,1, _geom.mesh(i,j-1), self_material, eps_self );
	double wp = vacuum_face_coefficient( i,j,0, i,j+1,0, +1,1, _geom.mesh(i,j+1), self_material, eps_self );
	set_link( a, _n2d(i,j-1), wm );
	set_link( a, _n2d(i,j+1), wp );
	cof += wm+wp;
    }

    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_neumann_node_cyl( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    // See add_neumann_node_1d().
    uint32_t self_material = _geom.dielectric_material_at( i, j, 0 );
    double eps_self = ( self_material == 0 ) ? 1.0 : _geom.get_boundary( self_material ).value();
    double cof = 0.0;

    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    double w = 2.0*eps_self;
	    set_link( a, _n2d(i+1,j), w );
	    (*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(1).value( x );
	    cof += w;
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    double w = 2.0*eps_self;
	    set_link( a, _n2d(i-1,j), w );
	    (*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(2).value( x );
	    cof += w;
	} else {
	    double wm = vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), self_material, eps_self );
	    double wp = vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), self_material, eps_self );
	    set_link( a, _n2d(i-1,j), wm );
	    set_link( a, _n2d(i+1,j), wp );
	    cof += wm+wp;
	}
	// On-axis cylindrical factor of 4 -- a coordinate-singularity
	// term, not an ordinary face, so (unlike the other links here)
	// this does not attempt a bisected cut-cell treatment if a
	// dielectric interface happens to sit exactly on the axis --
	// scoped the same way as a dielectric-Dirichlet interface is
	// elsewhere (see vacuum_face_coefficient()'s doc comment).
	double wr = 4.0*eps_self;
	set_link( a, _n2d(i,j+1), wr );
	cof += wr;
    } else {
	// Off-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    double w = 2.0*eps_self;
	    set_link( a, _n2d(i+1,j), w );
	    (*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(1).value( x );
	    cof += w;
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    double w = 2.0*eps_self;
	    set_link( a, _n2d(i-1,j), w );
	    (*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(2).value( x );
	    cof += w;
	} else {
	    double wm = vacuum_face_coefficient( i,j,0, i-1,j,0, -1,0, _geom.mesh(i-1,j), self_material, eps_self );
	    double wp = vacuum_face_coefficient( i,j,0, i+1,j,0, +1,0, _geom.mesh(i+1,j), self_material, eps_self );
	    set_link( a, _n2d(i-1,j), wm );
	    set_link( a, _n2d(i+1,j), wp );
	    cof += wm+wp;
	}

	if( bindex & EPOT_SOLVER_BYMAX ) {
	    double w = 2.0*eps_self;
	    set_link( a, _n2d(i,j-1), w );
	    (*_fd_vec)(a) += (1.0+0.5/j)*w*_geom.h()*_geom.get_boundary(2).value( x );
	    cof += w;
	} else {
	    double wm = (1.0-0.5/j)*vacuum_face_coefficient( i,j,0, i,j-1,0, -1,1, _geom.mesh(i,j-1), self_material, eps_self );
	    double wp = (1.0+0.5/j)*vacuum_face_coefficient( i,j,0, i,j+1,0, +1,1, _geom.mesh(i,j+1), self_material, eps_self );
	    set_link( a, _n2d(i,j-1), wm );
	    set_link( a, _n2d(i,j+1), wp );
	    cof += wm+wp;
	}
    }

    set_link( a, _n2d(i,j), -cof );
}


void EpotMatrixSolver::add_neumann_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j,k);

    // See add_neumann_node_1d().
    uint32_t self_material = _geom.dielectric_material_at( i, j, k );
    double eps_self = ( self_material == 0 ) ? 1.0 : _geom.get_boundary( self_material ).value();
    double cof = 0.0;

    if( bindex & EPOT_SOLVER_BXMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i+1,j,k), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(1).value( x );
	cof += w;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	double w = 2.0*eps_self;
	set_link( a,_n2d(i-1,j,k), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(2).value( x );
	cof += w;
    } else {
	double wm = vacuum_face_coefficient( i,j,k, i-1,j,k, -1,0, _geom.mesh(i-1,j,k), self_material, eps_self );
	double wp = vacuum_face_coefficient( i,j,k, i+1,j,k, +1,0, _geom.mesh(i+1,j,k), self_material, eps_self );
	set_link( a, _n2d(i-1,j,k), wm );
	set_link( a, _n2d(i+1,j,k), wp );
	cof += wm+wp;
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j+1,k), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(3).value( x );
	cof += w;
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j-1,k), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(4).value( x );
	cof += w;
    } else {
	double wm = vacuum_face_coefficient( i,j,k, i,j-1,k, -1,1, _geom.mesh(i,j-1,k), self_material, eps_self );
	double wp = vacuum_face_coefficient( i,j,k, i,j+1,k, +1,1, _geom.mesh(i,j+1,k), self_material, eps_self );
	set_link( a, _n2d(i,j-1,k), wm );
	set_link( a, _n2d(i,j+1,k), wp );
	cof += wm+wp;
    }

    if( bindex & EPOT_SOLVER_BZMIN ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j,k+1), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(5).value( x );
	cof += w;
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	double w = 2.0*eps_self;
	set_link( a, _n2d(i,j,k-1), w );
	(*_fd_vec)(a) += w*_geom.h()*_geom.get_boundary(6).value( x );
	cof += w;
    } else {
	double wm = vacuum_face_coefficient( i,j,k, i,j,k-1, -1,2, _geom.mesh(i,j,k-1), self_material, eps_self );
	double wp = vacuum_face_coefficient( i,j,k, i,j,k+1, +1,2, _geom.mesh(i,j,k+1), self_material, eps_self );
	set_link( a, _n2d(i,j,k-1), wm );
	set_link( a, _n2d(i,j,k+1), wp );
	cof += wm+wp;
    }

    set_link( a, _n2d(i,j,k), -cof );
}


void EpotMatrixSolver::add_neumann_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    switch( _geom.geom_mode() ) {
    case MODE_1D:
	add_neumann_node_1d( i, x );
        break;
    case MODE_2D:
	add_neumann_node_2d( i, j, x );
        break;
    case MODE_CYL:
	add_neumann_node_cyl( i, j, x );
	break;
    case MODE_3D:
	add_neumann_node_3d( i, j, k, x );
        break;
    }

    // Nonlinear (plasma) contribution to rhs/diagonal is epot-dependent
    // and is added separately, once per Newton iteration, by
    // update_nonlinear_node() -- see build_mat_vec(). Everything in
    // this function is geometry-only and gets cached after the first
    // build.

    // Right hand side
    if( _plasma != PLASMA_SHIELD )
	(*_fd_vec)(a) += -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
}


void EpotMatrixSolver::reset_matrix( void )
{
    if( _fd_mat )
        delete _fd_mat;
    if( _fd_vec )
        delete _fd_vec;
    if( _fd_vec_base )
        delete _fd_vec_base;
    _fd_mat = 0;
    _fd_vec = 0;
    _fd_vec_base = 0;
    _dof = 0;
    _linear_built = false;
    _fd_mat_diag0.clear();
    // Must be dropped together with the matrix these index into -- see
    // _fd_mat_diag_idx's doc comment.
    _fd_mat_diag_idx.clear();

    _n2d.clear();
}


void EpotMatrixSolver::set_initial_guess( const MeshScalarField &epot, Vector &X ) const
{
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    X.resize( _dof );
    for( int32_t a = 0; a < (int32_t)_geom.nodecount(); a++ ) {
	if( (_n2d(a) & N2D_TYPE_MASK) == N2D_TYPE_FREE )
	    X( _n2d(a) & N2D_INDEX_MASK ) = epot(a);
    }
}


void EpotMatrixSolver::set_solution( MeshScalarField &epot, const Vector &X ) const
{ 
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    for( uint32_t a = 0; a < _geom.nodecount(); a++ ) {
	uint32_t b = _n2d(a);
        if( (b & N2D_TYPE_MASK) == N2D_TYPE_FREE )
            epot(a) = X(b);
    }
}


void EpotMatrixSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    _epot = &epot;
    _scharge = &scharge;
    EpotSolver::preprocess( epot );

    reset_matrix();

    // Build n2d array and calculate degrees of freedom.
    //
    // NOTE the two index spaces stored here, which are NOT interchangeable:
    //   fixed node -> N2D_TYPE_FIXED | (MESH index)
    //   free  node -> N2D_TYPE_FREE  | (DOF  index)
    // set_link()'s row argument is a DOF index, its column argument for a
    // fixed neighbour is a MESH index. Mixing them silently produces valid-
    // looking but meaningless matrix columns.
    _n2d.resize( _geom.size() );
    _dof2node.clear();
    _dof2node.reserve( _geom.nodecount() );
    _dof = 0;
    for( int32_t a = 0; a < (int32_t)_geom.nodecount(); a++ ) {

	uint32_t mesh = _geom.mesh(a);
	bool fixed = mesh & SMESH_NODE_FIXED;

	if( fixed )
	    // Fixed node
	    _n2d(a) = N2D_TYPE_FIXED | a;
	else { 
	    // Free node
	    _n2d(a) = N2D_TYPE_FREE | _dof;
	    _dof2node.push_back( a );          // inverse map, see set_link()
	    _dof++;
	}
    }

    ibsimu.message( 1 ) << "dof = " << _dof << "\n";
    if( _dof == 0 )
        throw( Error( ERROR_LOCATION, "zero degrees of freedom" ) );
    
    // Allocate problem matrix and vector
    _fd_mat = new CRowMatrix( _dof, _dof );
    _fd_vec = new Vector( _dof );

    // Sized once here rather than per get_resjac() call; contents are
    // reset at the top of each of those calls.
    _d_vec.resize( _dof );
    _d_vec.clear();
}


void EpotMatrixSolver::build_mat_vec( void )
{
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    int nthreads = (int)ibsimu.get_thread_count();

    // Every mesh node's dielectric material used to be precomputed here,
    // once ever, into a solver-private _dielectric_mat cache (each entry
    // populated via self_material_at()'s Geometry::inside() fallback,
    // since the raw tag alone wasn't always enough -- see
    // Geometry::dielectric_material_at()'s doc comment). That's no
    // longer needed: Geometry now populates the same information into
    // its own material() array for free during mesh build, so
    // update_nonlinear_node()/vacuum_face_coefficient()/add_vacuum_node()
    // etc. just call Geometry::dielectric_material_at(i,j,k) directly.

    if( !_linear_built ) {

	auto time_t0 = Clock::now();

	// One-time build of the linear (geometry-only) part: matrix
	// stencil coefficients, Dirichlet-neighbour rhs contributions
	// and the constant space-charge rhs term. None of this depends
	// on the current solution guess X, only on the mesh, boundary
	// conditions and scharge, all of which are fixed between here
	// and the next preprocess()/reset_matrix() call -- so it is
	// built once and cached, instead of being redone on every
	// Newton iteration.
	_fd_mat->clear();
	_fd_vec->clear();

	// Reset per-row entry buffers (see header doc comment on
	// _row_entries for why these exist instead of writing straight
	// into _fd_mat from the parallel loop below). Reserve enough
	// capacity for the largest possible stencil (3D: 6 neighbours + 1
	// combined centre term = 7) so no row needs to reallocate while
	// another thread is concurrently working on a different row.
	if( _row_entries.size() != _dof )
	    _row_entries.resize( _dof );
	for( uint32_t a = 0; a < _dof; a++ ) {
	    _row_entries[a].clear();
	    _row_entries[a].reserve( 7 );
	}

	// Parallel over grid nodes: every node is handled by exactly
	// one thread, which only ever writes to that node's own
	// _row_entries[a] slot and its own _fd_vec(a) entry (a is
	// unique per node), so there is no cross-row data race. This
	// does mean any user-supplied Boundary::value() callback
	// reachable from add_vacuum_node()/add_near_solid_node()/
	// add_neumann_node() must itself be safe to call concurrently
	// from multiple threads -- this parallel loop calls it with no
	// locking.
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
	for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		    Vec3D x( _geom.origo(0) + _geom.h()*i,
			     _geom.origo(1) + _geom.h()*j,
			     _geom.origo(2) + _geom.h()*k );

		    uint32_t a = (k*_geom.size(1)+j)*_geom.size(0)+i;
		    uint32_t mesh = _geom.mesh(a);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    bool fixed = mesh & SMESH_NODE_FIXED;

		    if( fixed ) {
			//ibsimu.message( 1 ) << "Fixed\n";
			continue;
		    } else if( node_id == SMESH_NODE_ID_PURE_VACUUM ) {
			//ibsimu.message( 1 ) << "Vacuum\n";
			add_vacuum_node( i, j, k, x );
		    } else if( node_id == SMESH_NODE_ID_NEAR_SOLID ) {
			//ibsimu.message( 1 ) << "Near solid\n";
			add_near_solid_node( i, j, k, x );
		    } else if( node_id == SMESH_NODE_ID_NEUMANN ) {
			//ibsimu.message( 1 ) << "Neumann\n";
			add_neumann_node( i, j, k, x );
		    }
		}
	    }
	}

	// Compact the per-row buffers into _fd_mat via construct_add(), in
	// ascending row order, exactly as construct_add() requires. This
	// pass is pure data movement (no stencil math, no plasma
	// evaluation) over already-computed entries, so it stays cheap
	// even though it runs single-threaded.
	//
	// Duplicate columns within a row are MERGED here rather than being
	// handed to construct_add() twice. construct_add() appends
	// unconditionally -- it does not look for an existing entry -- so
	// duplicates would survive into the CSR arrays. A matrix-vector
	// product sums them and so stays correct, but everything that
	// SEARCHES for a particular entry sees only the first one:
	// _fd_mat_diag_idx below (which get_resjac() uses to apply the
	// nonlinear plasma term to the diagonal every Newton iteration),
	// GMG_Precond::find_diagonals(), and the ILU/RBSOR preconditioners'
	// diagonal handling would all then act on a fraction of the true
	// diagonal.
	//
	// Set_link()'s Neumann-mask redirect is what makes this reachable:
	// it emits an extra entry at (a,a) for every masked face, on top of
	// the diagonal the stencil builder emits at the end. Merging here
	// rather than relying on the builders always emitting the diagonal
	// last keeps that independent of call order.
	std::vector<double> rowacc;
	std::vector<int32_t> rowcol;
	for( uint32_t a = 0; a < _dof; a++ ) {
	    const std::vector<std::pair<int32_t,double> > &row = _row_entries[a];

	    rowacc.clear();
	    rowcol.clear();
	    for( size_t n = 0; n < row.size(); n++ ) {
		size_t m = 0;
		for( ; m < rowcol.size(); m++ ) {
		    if( rowcol[m] == row[n].first ) {
			rowacc[m] += row[n].second;
			break;
		    }
		}
		if( m == rowcol.size() ) {
		    rowcol.push_back( row[n].first );
		    rowacc.push_back( row[n].second );
		}
	    }

	    for( size_t m = 0; m < rowcol.size(); m++ )
		_fd_mat->construct_add( a, rowcol[m], rowacc[m] );
	}

	// Order matrix
	_fd_mat->order_ascending();

	// Snapshot the constant rhs so later calls can restore it
	// instead of rebuilding it from scratch.
	if( _fd_vec_base )
	    delete _fd_vec_base;
	_fd_vec_base = new Vector( *_fd_vec );

	// Snapshot the pure-linear diagonal J0(a,a) so get_resjac() can
	// re-derive the Jacobian diagonal fresh each Newton iteration
	// (J0(a,a) - D(a)) instead of repeatedly subtracting from
	// whatever the now-persistent _fd_mat's diagonal currently holds.
	//
	// Also resolve, once, where each diagonal element physically lives
	// in the matrix's value array, so get_resjac()'s two per-iteration
	// diagonal sweeps become direct indexed writes rather than
	// row-scanning set(a,a) calls -- see _fd_mat_diag_idx's doc
	// comment. Done after order_ascending() above, and the structure
	// is fixed from here until the next reset_matrix(), so these
	// indices stay valid for the whole solve.
	_fd_mat_diag0.resize( _dof );
	_fd_mat_diag_idx.resize( _dof );
	for( uint32_t a = 0; a < _dof; a++ ) {
	    int idx = _fd_mat->value_index( a, a );
	    if( idx < 0 )
		// Every free row gets an explicit centre-node link from
		// add_vacuum_node()/add_near_solid_node()/add_neumann_node(),
		// so a structurally missing diagonal means the stencil for
		// this row was never emitted -- a real bug upstream, and one
		// that would otherwise show up much later as a silently
		// wrong Jacobian rather than an error here.
		throw( Error( ERROR_LOCATION, "no diagonal element on matrix row "
			      + to_string(a) ) );
	    _fd_mat_diag_idx[a] = idx;
	    _fd_mat_diag0[a] = _fd_mat->value_ptr()[idx];
	}

	_linear_built = true;
	_time_linbuild += std::chrono::duration<double>( Clock::now()-time_t0 ).count();

    } else {

	// Matrix and constant rhs part are unchanged since the last
	// build (same mesh, boundary conditions and scharge) -- just
	// restore the cached constant rhs instead of rebuilding it.
	*_fd_vec = *_fd_vec_base;
    }

    // Nonlinear update: add the epot-dependent plasma rhs term and
    // diagonal derivative for every free node. This is the only part
    // of matrix/rhs construction that actually depends on X (via
    // _sol), so it is the only part that has to run every Newton
    // iteration. Skipped entirely when the solver isn't using a
    // plasma model, and cheap even when it does run since it does none
    // of the stencil/boundary bookkeeping the linear build does -- just
    // a mesh scan plus, for plasma-active nodes, the actual nonlinear
    // evaluation (pexp_newton/nsimp_newton/shield_newton).
    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP ||
	_plasma == PLASMA_SHIELD || _comp ) {

	auto time_t0 = Clock::now();

#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
	for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		    uint32_t a_full = (k*_geom.size(1)+j)*_geom.size(0)+i;
		    uint32_t mesh = _geom.mesh(a_full);
		    if( mesh & SMESH_NODE_FIXED )
			continue;

		    Vec3D x( _geom.origo(0) + _geom.h()*i,
			     _geom.origo(1) + _geom.h()*j,
			     _geom.origo(2) + _geom.h()*k );

		    uint32_t a = _n2d(a_full) & N2D_INDEX_MASK;
		    update_nonlinear_node( a, i, j, k, x );
		}
	    }
	}

	_time_nonlin += std::chrono::duration<double>( Clock::now()-time_t0 ).count();
    }
}


void EpotMatrixSolver::postprocess( void )
{ 
    reset_matrix();
    EpotSolver::postprocess();
}


void EpotMatrixSolver::get_vecmat( const CRowMatrix **A, const Vector **B )
{
    build_mat_vec();
    *A = _fd_mat;
    *B = _fd_vec;
}


void EpotMatrixSolver::get_resjac( const CRowMatrix **J, const Vector **R, const Vector &X )
{
    // Build residual and jacobian from linear matric and vector.
    // Calculate R = J0*X - B(X) and J = J0 + I*D(X)
    //
    // _d_vec is a member resized once in preprocess(), not a per-call
    // allocation. Cleared here so rows that update_nonlinear_node()
    // never visits (every row, when no plasma model is active) still
    // read back as exactly zero, exactly as a freshly calloc'd Vector
    // used to.
    _d_vec.clear();

    // Construct whole right hand side to _fd_vec, nonlinear component of diagonal
    // to _d_vec and linear part of jacobian to _fd_mat.
    _sol = &X;
    build_mat_vec();

    // _fd_mat is now a cached, persistent matrix (see build_mat_vec()):
    // whatever the *previous* get_resjac() call left in its diagonal is
    // the Jacobian diagonal J0(a,a) - D_prev(a), not the pure linear
    // J0(a,a). Restore the pure linear diagonal before using _fd_mat as
    // J0 for the w = J0*X multiply below, otherwise w would silently
    // pick up the previous iteration's nonlinear correction.
    //
    // Written straight into the value array via the diagonal positions
    // cached in _fd_mat_diag_idx rather than through set(a,a), which
    // would re-scan row a to find the diagonal on every one of these
    // 2*_dof accesses -- see _fd_mat_diag_idx's doc comment. Structure
    // is fixed here (no insertions since order_ascending()), so the
    // cached indices are valid; values_changed() below tells the matrix
    // its cached MKL handle is stale, once for the whole sweep instead
    // of once per element.
    double *mval = _fd_mat->value_ptr();
    for( uint32_t a = 0; a < _dof; a++ )
	mval[_fd_mat_diag_idx[a]] = _fd_mat_diag0[a];
    _fd_mat->values_changed();

    // Linear part
    Vector w = (*_fd_mat) * X;

    for( uint32_t a = 0; a < _dof; a++ ) {

	(*_fd_vec)(a) = w(a) - (*_fd_vec)(a);
	// Set (not subtract from) the diagonal: it must be re-derived
	// fresh from the pure-linear J0(a,a) each call, not adjusted
	// relative to whatever a previous iteration left there.
	mval[_fd_mat_diag_idx[a]] = _fd_mat_diag0[a] - _d_vec(a);
    }
    _fd_mat->values_changed();

    *J = _fd_mat;
    *R = _fd_vec;
}


/*
bool EpotMatrixSolver::linear( void ) const
{
    if( _plasma == PLASMA_NONE || _plasma == PLASMA_PEXP_INITIAL || _plasma == PLASMA_NSIMP_INITIAL )
        return( true );
    else
        return( false );
}
*/


void EpotMatrixSolver::debug_print( std::ostream &os ) const
{
    EpotSolver::debug_print( os );
    os << "**EpotMatrixSolver\n";
    os << "dof = " << _dof << "\n";
    _n2d.debug_print( os );
    if( _fd_mat )
	os << "fd_mat = \n" << *_fd_mat << "\n";
    else
	os << "fd_mat = NULL\n";
    if( _fd_vec )
	os << "fd_vec = \n" << *_fd_vec << "\n";
    else
	os << "fd_vec = NULL\n";
}

// epot_matrixsolver.cpp
void EpotMatrixSolver::build_node_map( std::vector<int32_t> &map ) const
{
    map.resize( _geom.nodecount() );
    for( uint32_t a = 0; a < _geom.nodecount(); a++ ) {
        uint32_t b = _n2d(a);
        if( (b & N2D_TYPE_MASK) == N2D_TYPE_FREE )
            map[a] = (int32_t)(b & N2D_INDEX_MASK);
        else
            map[a] = -1;
    }
}


void EpotMatrixSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}
