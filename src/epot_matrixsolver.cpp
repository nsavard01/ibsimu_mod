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
      _cm_i0(0), _cm_j0(0), _cm_k0(0), _cm_ni(0), _cm_nj(0), _cm_nk(0),
      _sol(0),
      _key_revision(0),
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
}


/* Add the term val*phi(b) to the equation of free node a.
 *
 * a is the row and b an entry of _n2d:
 *   - a free node: val is added to matrix element (a,b).
 *   - a fixed node of known potential: the term goes to the right hand
 *     side.
 *   - a node inside a Neumann mask: the mask surface is a symmetry plane
 *     through node a, as a Neumann face of the simulation box is, so the
 *     term goes to the mirror image of b about a. If there is no image
 *     (off the mesh, or masked as well) it goes to the diagonal, which
 *     puts the symmetry plane half a cell further out.
 */
void EpotMatrixSolver::set_link( uint32_t a, uint32_t b, double val )
{
    const int32_t row = (int32_t)(a & N2D_INDEX_MASK);
    const int32_t nb  = (int32_t)(b & N2D_INDEX_MASK);

    // The entries of a row are collected in _row_entries, see build_linear().
    if( (b & N2D_TYPE_MASK) != N2D_TYPE_FIXED ) {
	_row_entries[row].push_back( std::make_pair( nb, val ) );
	return;
    }

    // For a fixed node nb is a mesh node index, not a row.
    if( !SMESH_NODE_IS_NEUMANN_MASK( _geom.mesh( nb ) ) ) {
	(*_fd_vec)(row) += -val * (*_epot)( nb );
	return;
    }

    // Mirror image of b about a, which must be on the same mesh line.
    const int32_t node_a = _dof2node[row];
    const int32_t node_m = 2*node_a - nb;
    const int32_t nx = _geom.size(0), nxy = nx*_geom.size(1);
    const int32_t d = ( nb > node_a ? nb-node_a : node_a-nb );
    if( node_m >= 0 && node_m < (int32_t)_geom.nodecount() &&
	( d == 1 ? node_m/nx == node_a/nx : d == nx ? node_m/nxy == node_a/nxy : true ) ) {
	const uint32_t n2d_m = _n2d(node_m);
	if( (n2d_m & N2D_TYPE_MASK) != N2D_TYPE_FIXED ) {
	    _row_entries[row].push_back( std::make_pair( (int32_t)(n2d_m & N2D_INDEX_MASK), val ) );
	    return;
	} else if( !SMESH_NODE_IS_NEUMANN_MASK( _geom.mesh( node_m ) ) ) {
	    (*_fd_vec)(row) += -val * (*_epot)( node_m );
	    return;
	}
    }
    _row_entries[row].push_back( std::make_pair( row, val ) );
}


void EpotMatrixSolver::get_comp_mask( MeshScalarField &out ) const
{
    const uint32_t nx = _geom.size(0), ny = _geom.size(1), nz = _geom.size(2);
    for( uint32_t k = 0; k < nz; k++ )
	for( uint32_t j = 0; j < ny; j++ )
	    for( uint32_t i = 0; i < nx; i++ ) {
		double v = 0.0;
		if( i >= _cm_i0 && i < _cm_i0+_cm_ni && j >= _cm_j0 && j < _cm_j0+_cm_nj &&
		    k >= _cm_k0 && k < _cm_k0+_cm_nk ) {
		    size_t a = ((size_t)(k-_cm_k0)*_cm_nj + (j-_cm_j0))*_cm_ni + (i-_cm_i0);
		    if( a < _comp_mask.size() && _comp_mask[a] ) v = 1.0;
		}
		out( i, j, k ) = v;
	    }
}


void EpotMatrixSolver::update_comp_mask( void )
{
    const uint32_t nx = _geom.size(0), ny = _geom.size(1), nz = _geom.size(2);
    const double h = _geom.h();

    // Node window: the bounding box of the compensation region if set,
    // otherwise the whole mesh. The work arrays below cover the window only.
    int64_t i0 = 0, j0 = 0, k0 = 0, i1 = (int64_t)nx-1, j1 = (int64_t)ny-1, k1 = (int64_t)nz-1;
    if( _comp.bbox_set ) {
	const uint32_t n[3] = { nx, ny, nz };
	int64_t lo[3], hi[3];
	for( int a = 0; a < 3; a++ ) {
	    lo[a] = (int64_t)std::floor( (_comp.bmin[a] - _geom.origo(a)) / h );
	    hi[a] = (int64_t)std::ceil ( (_comp.bmax[a] - _geom.origo(a)) / h );
	    if( lo[a] < 0 ) lo[a] = 0;
	    if( hi[a] > (int64_t)n[a]-1 ) hi[a] = (int64_t)n[a]-1;
	}
	i0=lo[0]; i1=hi[0]; j0=lo[1]; j1=hi[1]; k0=lo[2]; k1=hi[2];
	if( i1 < i0 || j1 < j0 || k1 < k0 ) { _comp_mask.clear(); _cm_ni=_cm_nj=_cm_nk=0; return; }
    }
    _cm_i0 = (uint32_t)i0; _cm_j0 = (uint32_t)j0; _cm_k0 = (uint32_t)k0;
    _cm_ni = (uint32_t)(i1-i0+1); _cm_nj = (uint32_t)(j1-j0+1); _cm_nk = (uint32_t)(k1-k0+1);
    const size_t N = (size_t)_cm_ni*_cm_nj*_cm_nk;
    _comp_mask.assign( N, 0 );

    auto loc = [&]( int64_t i, int64_t j, int64_t k ) -> size_t {
	return ((size_t)(k-k0)*_cm_nj + (size_t)(j-j0))*_cm_ni + (size_t)(i-i0);
    };
    auto inbox = [&]( int64_t i, int64_t j, int64_t k ) -> bool {
	return( i>=i0 && i<=i1 && j>=j0 && j<=j1 && k>=k0 && k<=k1 );
    };
    auto phi_at = [&]( int64_t i, int64_t j, int64_t k ) -> double {
	size_t a_full = ((size_t)k*ny + j)*nx + i;
	uint32_t mesh = _geom.mesh( a_full );
	if( mesh & SMESH_NODE_FIXED )
	    return (*_epot)( a_full );
	return (*_sol)( _n2d( a_full ) & N2D_INDEX_MASK );
    };
    auto in_region = [&]( int64_t i, int64_t j, int64_t k ) -> bool {
	if( !_comp.region ) return( true );
	Vec3D x( _geom.origo(0) + h*i, _geom.origo(1) + h*j, _geom.origo(2) + h*k );
	return( (*_comp.region)(x) );
    };

    // next[] = window index of the neighbour of highest potential;
    // -1 if the node is a maximum, -2 if that neighbour is in a solid or
    // outside the window or the region.
    std::vector<int64_t> next( N, -1 );
    const int di[6] = {-1,1,0,0,0,0}, dj[6] = {0,0,-1,1,0,0}, dk[6] = {0,0,0,0,-1,1};
#pragma omp parallel for num_threads(ibsimu.get_thread_count()) collapse(3) schedule(static)
    for( int64_t k = k0; k <= k1; k++ ) {
	for( int64_t j = j0; j <= j1; j++ ) {
	    for( int64_t i = i0; i <= i1; i++ ) {
		size_t a = loc(i,j,k);
		if( _geom.material(i,j,k) != 0 || !in_region(i,j,k) ) { next[a] = -2; continue; }
		double best = phi_at(i,j,k);
		int64_t bi = -1;
		for( int d = 0; d < 6; d++ ) {
		    int64_t ii=i+di[d], jj=j+dj[d], kk=k+dk[d];
		    if( ii<0 || ii>=(int64_t)nx || jj<0 || jj>=(int64_t)ny ||
			kk<0 || kk>=(int64_t)nz ) continue;
		    double pb = phi_at(ii,jj,kk);
		    if( pb <= best ) continue;
		    best = pb;
		    bi = ( !inbox(ii,jj,kk) || _geom.material(ii,jj,kk) != 0 ||
			   !in_region(ii,jj,kk) ) ? -2 : (int64_t)loc(ii,jj,kk);
		}
		next[a] = bi;
	    }
	}
    }

    // Follow next[] from every node. Nodes whose path ends at a maximum
    // are confined.
    std::vector<uint8_t> done( N, 0 );
    std::vector<size_t> chain;
    for( size_t a0 = 0; a0 < N; a0++ ) {
	if( done[a0] || next[a0] == -2 ) { done[a0] = 1; continue; }
	chain.clear();
	size_t a = a0;
	uint8_t res;
	for(;;) {
	    if( done[a] ) { res = _comp_mask[a]; break; }
	    chain.push_back( a ); done[a] = 1;
	    int64_t n = next[a];
	    if( n == -1 ) { res = 1; break; }
	    if( n == -2 ) { res = 0; break; }
	    a = (size_t)n;
	}
	for( size_t c : chain ) _comp_mask[c] = res;
    }
}


/* Plasma and compensation terms of free node a, mesh node (i,j,k) at x,
 * at the current solution: the charge term is added to the right hand side
 * and its derivative is written to _d_vec. Called for every free node on
 * every Newton iteration.
 */
void EpotMatrixSolver::update_nonlinear_node( uint32_t a, uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    const bool have_plasma = ( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP ||
                               _plasma == PLASMA_SHIELD );
    if( !have_plasma && !_comp.enabled )
	return;

    // There is no plasma inside a dielectric.
    const bool dielectric_free = ( _geom.dielectric_material_at( i, j, k ) == 0 );

    // The plasma model and the compensation model are added up. Each has
    // its own region function, and one mesh can hold both a source plasma
    // and a compensated drift.
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

    if( dielectric_free && _comp.enabled ) {
	bool incomp = true;
	if( _comp.region )
	    incomp = (*_comp.region)(x);
	if( incomp && _comp.ascent ) {
	    // Only where update_comp_mask() found the electrons confined.
	    if( i < _cm_i0 || i >= _cm_i0+_cm_ni || j < _cm_j0 || j >= _cm_j0+_cm_nj ||
		k < _cm_k0 || k >= _cm_k0+_cm_nk )
		incomp = false;
	    else {
		size_t a_loc = ((size_t)(k-_cm_k0)*_cm_nj + (j-_cm_j0))*_cm_ni + (i-_cm_i0);
		if( a_loc >= _comp_mask.size() || !_comp_mask[a_loc] )
		    incomp = false;
	    }
	}
	if( incomp ) {
	    // The electron density is rho_e = -rho_i*g with the local ion
	    // density rho_i, so _scharge must hold the ion charge only.
	    double rho_i = (*_scharge)(i,j,k);
	    if( rho_i > _comp.rho_min ) {
		double p = (*_sol)(a);
		double g, dg;
		comp_newton( g, dg, p );
		double R = rho_i*_geom.h()*_geom.h()/EPSILON0;
		fd += R*g;
		d  += R*dg;
	    }
	}
    }

    (*_fd_vec)(a) += fd;
    _d_vec(a) = d;
}


double EpotMatrixSolver::node_epsilon_r( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint32_t solid_number = _geom.dielectric_material_at( i, j, k );
    if( solid_number == 0 )
	return( 1.0 );
    return( _geom.boundary( solid_number ).value() );
}


double EpotMatrixSolver::face_epsilon_alpha( double eps_self, double eps_neighbor, double alpha )
{
    return( 1.0/( alpha/eps_self + (1.0-alpha)/eps_neighbor ) );
}


double EpotMatrixSolver::vacuum_face_coefficient( const int32_t i[3], int coord, int sign,
						   uint32_t self_material, double eps_self ) const
{
    int32_t n[3] = { i[0], i[1], i[2] };
    n[coord] += sign;
    const uint32_t neighbor_mesh = _geom.mesh( n[0], n[1], n[2] );

    // Neumann mask: set_link() sends this link to the mirror node, in the
    // medium of this node.
    if( SMESH_NODE_IS_NEUMANN_MASK( neighbor_mesh ) )
	return( eps_self );

    if( (neighbor_mesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET ) {
	// A face of the simulation box is always a full cell away.
	const uint32_t boundary_number = neighbor_mesh & SMESH_BOUNDARY_NUMBER_MASK;
	if( boundary_number < 7 )
	    return( eps_self );

	// Electrode surface at alpha*h, reached through the medium of this
	// node.
	const double alpha = _geom.solid_face_frac( i[0], i[1], i[2], boundary_number, sign, coord );
	return( eps_self / alpha );
    }

    // dielectric_material_at() is right also for a neighbour whose node tag
    // says something else (Neumann, near-solid).
    const uint32_t nb_material = _geom.dielectric_material_at( n[0], n[1], n[2] );
    if( nb_material == self_material )
	return( eps_self );

    // Material interface at alpha*h: the two media in series.
    // solid_face_frac() has to start from a node outside the solid, so if
    // it is this node that is inside, the search runs from the neighbour
    // and the fraction is complemented. Where two dielectrics touch, the
    // surface of this node's solid is used.
    const double eps_neighbor = ( nb_material != 0 ) ? _geom.boundary( nb_material ).value() : 1.0;
    double alpha;
    if( self_material != 0 )
	alpha = 1.0 - _geom.solid_face_frac( n[0], n[1], n[2], self_material, -sign, coord );
    else
	alpha = _geom.solid_face_frac( i[0], i[1], i[2], nb_material, sign, coord );
    return( face_epsilon_alpha( eps_self, eps_neighbor, alpha ) );
}


/* Stencil of a node with no conductor surface and no box face within a
 * cell. The node is in vacuum or inside a dielectric.
 */
void EpotMatrixSolver::add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    (void)x;
    const uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    const int32_t ijk[3] = { (int32_t)i, (int32_t)j, (int32_t)k };
    const uint32_t dim = _geom.dim();
    const bool cyl = ( _geom.geom_mode() == MODE_CYL );
    const uint32_t self_material = _geom.dielectric_material_at( i, j, k );
    const double eps_self = node_epsilon_r( i, j, k );

    // Weight of the link to the lower/upper neighbour on each axis. With no
    // dielectric nearby this is the standard stencil, all weights 1 (times
    // 1-+1/(2j) for the radial axis of a cylindrical geometry).
    double w[3][2];
    double cof = 0.0;
    for( uint32_t c = 0; c < dim; c++ ) {
	w[c][0] = vacuum_face_coefficient( ijk, c, -1, self_material, eps_self );
	w[c][1] = vacuum_face_coefficient( ijk, c, +1, self_material, eps_self );
	if( cyl && c == 1 ) {
	    w[c][0] *= (1.0-0.5/j);
	    w[c][1] *= (1.0+0.5/j);
	}
	cof += w[c][0];
	cof += w[c][1];
    }

    // Links in ascending column order.
    for( uint32_t c = dim; c-- > 0; ) {
	int32_t m[3] = { ijk[0], ijk[1], ijk[2] };
	m[c]--;
	set_link( a, _n2d(m[0],m[1],m[2]), w[c][0] );
    }
    for( uint32_t c = 0; c < dim; c++ ) {
	int32_t p[3] = { ijk[0], ijk[1], ijk[2] };
	p[c]++;
	set_link( a, _n2d(p[0],p[1],p[2]), w[c][1] );
    }
    set_link( a, _n2d(i,j,k), -cof );

    // The permittivity is in the weights, so the space charge term added
    // in preprocess() has no epsilon_r in it. A node inside a dielectric
    // takes no space charge at all.
    if( self_material != 0 )
	_row_scharge[a] = 0;
}


/* Links of near-solid node a to its two neighbours along axis coord, for
 * every axis but the radial one of a cylindrical geometry. Returns the
 * weight that these links add to the diagonal.
 *
 * bnd: bit 0/1 set if the node is on the lower/upper face of the box.
 * cut: bit 0/1 set if a conductor surface cuts the link to the lower/upper
 *      neighbour, at distance alpha/beta (fraction of h, 1 if not cut).
 */
double EpotMatrixSolver::near_solid_axis( uint32_t a, const int32_t i[3], int coord,
					   uint8_t bnd, uint8_t cut, double alpha, double beta,
					   const Vec3D &x )
{
    int32_t m[3] = { i[0], i[1], i[2] };
    int32_t p[3] = { i[0], i[1], i[2] };
    m[coord]--;
    p[coord]++;

    // Neumann face of the box: the neighbour outside is the mirror image
    // of the one inside.
    if( bnd & 0x01 ) {
	const double w = 2.0/(beta*beta);
	set_link( a, _n2d(p[0],p[1],p[2]), w );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1+2*coord).value(x) / beta;
	return( w );
    } else if( bnd & 0x02 ) {
	const double w = 2.0/(alpha*alpha);
	set_link( a, _n2d(m[0],m[1],m[2]), w );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2+2*coord).value(x) / alpha;
	return( w );
    }

    const bool cm = cut & 0x01;
    const bool cp = cut & 0x02;
    double wm, wp, cof;
    if( (cm || cp) &&
	(cm || _geom.dielectric_material_at(m[0],m[1],m[2]) == 0) &&
	(cp || _geom.dielectric_material_at(p[0],p[1],p[2]) == 0) ) {
	// Vacuum all the way to the conductor(s): Shortley-Weller.
	wm = 2.0/((alpha+beta)*alpha);
	wp = 2.0/((alpha+beta)*beta);
	cof = 2.0/(alpha*beta);
    } else {
	// No conductor on this axis, or a conductor on one side and a
	// dielectric on the other. Shortley-Weller assumes a smooth
	// potential through both neighbours, which a permittivity jump
	// breaks, so each face gets its own conductance instead. The node
	// itself is in vacuum, so that is 1/distance towards a conductor.
	wm = cm ? 1.0/alpha : vacuum_face_coefficient( i, coord, -1, 0, 1.0 );
	wp = cp ? 1.0/beta  : vacuum_face_coefficient( i, coord, +1, 0, 1.0 );
	cof = wm+wp;
    }
    set_link( a, _n2d(m[0],m[1],m[2]), wm );
    set_link( a, _n2d(p[0],p[1],p[2]), wp );
    return( cof );
}


/* The same as near_solid_axis() for the radial axis of a cylindrical
 * geometry.
 */
double EpotMatrixSolver::near_solid_radial( uint32_t a, const int32_t i[3],
					     uint8_t bnd, uint8_t cut, double alpha, double beta,
					     const Vec3D &x )
{
    const int32_t j = i[1];
    if( bnd & 0x01 ) {
	// On axis: the limit of the radial operator at r=0. A dielectric
	// touching the axis gets no interface treatment here.
	set_link( a, _n2d(i[0],j+1), 4.0 );
	return( 4.0 );
    } else if( bnd & 0x02 ) {
	// Neumann face at rmax, in the vacuum form.
	const double w = 2.0/(alpha*alpha);
	set_link( a, _n2d(i[0],j-1), w );
	(*_fd_vec)(a) += (1.0)/(2.0*alpha)*(2.0/alpha+1.0/j)*
	    2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
	return( w );
    }

    const bool cm = cut & 0x01;
    const bool cp = cut & 0x02;
    double wm, wp, cof;
    if( (cm || _geom.dielectric_material_at(i[0],j-1,0) == 0) &&
	(cp || _geom.dielectric_material_at(i[0],j+1,0) == 0) ) {
	// No dielectric neighbour: Shortley-Weller with the 1/r term.
	wm = 1.0/(alpha+beta)*(2.0/alpha-1.0/j);
	wp = 1.0/(alpha+beta)*(2.0/beta+1.0/j);
	cof = 2.0/(alpha*beta);
    } else {
	// A dielectric neighbour: the flux through each face of the
	// control volume from r_j-alpha*h/2 to r_j+beta*h/2 is
	// eps_face*r_face*dphi/dr, which in the row scaling of this file
	// (vacuum node: diagonal 4) gives the face weight
	//     eps_face * (r_face/r_j) / (distance to the neighbour / h).
	// For alpha=beta=1 this is the stencil of add_vacuum_node(), and
	// for large j that of near_solid_axis(). A face towards a
	// conductor has eps 1, its distance is in alpha/beta already.
	const double rj = (double)j;
	const double em = cm ? 1.0 : vacuum_face_coefficient( i, 1, -1, 0, 1.0 );
	const double ep = cp ? 1.0 : vacuum_face_coefficient( i, 1, +1, 0, 1.0 );
	wm = em*(rj - 0.5*alpha)/(rj*alpha);
	wp = ep*(rj + 0.5*beta )/(rj*beta );
	cof = wm+wp;
    }
    set_link( a, _n2d(i[0],j-1), wm );
    set_link( a, _n2d(i[0],j+1), wp );
    return( cof );
}


/* Stencil of a vacuum node with a conductor surface within a cell.
 */
void EpotMatrixSolver::add_near_solid_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    const uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    const int32_t ijk[3] = { (int32_t)i, (int32_t)j, (int32_t)k };
    const bool cyl = ( _geom.geom_mode() == MODE_CYL );
    const uint8_t bindex = boundary_index_general(i,j,k);

    // Near-solid data of the node: a flag byte with one bit per direction
    // (xmin, xmax, ymin, ...) whose link is cut by a conductor surface,
    // followed by one byte per set bit, the distance to the surface in
    // units of h/255.
    const uint8_t *ptr = _geom.nearsolid_ptr( _geom.mesh(i,j,k) & SMESH_NEAR_SOLID_INDEX_MASK );
    const uint8_t sflag = *ptr++;

    double cof = 0.0;
    for( uint32_t c = 0; c < _geom.dim(); c++ ) {
	const uint8_t bnd = (bindex >> 2*c) & 0x03;
	const uint8_t cut = (sflag >> 2*c) & 0x03;
	const double alpha = (cut & 0x01) ? *ptr++/255.0 : 1.0;
	const double beta  = (cut & 0x02) ? *ptr++/255.0 : 1.0;
	if( cyl && c == 1 )
	    cof += near_solid_radial( a, ijk, bnd, cut, alpha, beta, x );
	else
	    cof += near_solid_axis( a, ijk, c, bnd, cut, alpha, beta, x );
    }
    set_link( a, _n2d(i,j,k), -cof );
}


/* Stencil of a node on a Neumann face of the simulation box with no
 * conductor surface within a cell. The node is in vacuum or inside a
 * dielectric. The boundary value is the derivative of the potential along
 * the inward normal of the face.
 */
void EpotMatrixSolver::add_neumann_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    const uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    const int32_t ijk[3] = { (int32_t)i, (int32_t)j, (int32_t)k };
    const bool cyl = ( _geom.geom_mode() == MODE_CYL );
    const uint8_t bindex = boundary_index_general(i,j,k);
    const uint32_t self_material = _geom.dielectric_material_at( i, j, k );
    const double eps_self = node_epsilon_r( i, j, k );

    double cof = 0.0;
    for( uint32_t c = 0; c < _geom.dim(); c++ ) {
	const uint8_t bnd = (bindex >> 2*c) & 0x03;
	const bool radial = ( cyl && c == 1 );
	int32_t m[3] = { ijk[0], ijk[1], ijk[2] };
	int32_t p[3] = { ijk[0], ijk[1], ijk[2] };
	m[c]--;
	p[c]++;

	if( radial && (bnd & 0x01) ) {
	    // On axis: the limit of the radial operator at r=0. A
	    // dielectric touching the axis gets no interface treatment.
	    const double w = 4.0*eps_self;
	    set_link( a, _n2d(p[0],p[1],p[2]), w );
	    cof += w;
	} else if( bnd ) {
	    // On the face: the neighbour outside is the mirror image of
	    // the one inside.
	    const bool lower = bnd & 0x01;
	    const double w = 2.0*eps_self;
	    set_link( a, lower ? _n2d(p[0],p[1],p[2]) : _n2d(m[0],m[1],m[2]), w );
	    (*_fd_vec)(a) += ( radial ? (1.0+0.5/j)*w : w )*_geom.h()*
		_geom.get_boundary( lower ? 1+2*c : 2+2*c ).value( x );
	    cof += w;
	} else {
	    double wm = vacuum_face_coefficient( ijk, c, -1, self_material, eps_self );
	    double wp = vacuum_face_coefficient( ijk, c, +1, self_material, eps_self );
	    if( radial ) {
		wm *= (1.0-0.5/j);
		wp *= (1.0+0.5/j);
	    }
	    set_link( a, _n2d(m[0],m[1],m[2]), wm );
	    set_link( a, _n2d(p[0],p[1],p[2]), wp );
	    cof += wm+wp;
	}
    }
    set_link( a, _n2d(i,j,k), -cof );

    // No space charge inside a dielectric, as in add_vacuum_node().
    if( self_material != 0 )
	_row_scharge[a] = 0;
}


void EpotMatrixSolver::reset_matrix( void )
{
    delete _fd_mat;
    delete _fd_vec;
    _fd_mat = 0;
    _fd_vec = 0;
    _dof = 0;
    _fd_mat_diag0.clear();
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


/* Has anything the linear system was built from changed since the last call?
 *
 * The matrix, the node numbering and the boundary part of the right hand
 * side depend on the geometry definition, on the node tags as left by
 * EpotSolver::preprocess() (which fixes nodes for forced potentials and
 * initial plasma) and on the potential of the fixed nodes. Solver settings
 * are not checked here: their setters call reset_problem().
 *
 * The current state is recorded as the reference for the next call.
 */
bool EpotMatrixSolver::linear_inputs_changed( void )
{
    const uint32_t n = _geom.nodecount();
    bool changed = ( _key_revision != _geom.revision() || _key_mesh.size() != n );
    _key_revision = _geom.revision();
    _key_mesh.resize( n );
    _key_epot.resize( n );

    // A Neumann value given as a functor can change without the geometry
    // knowing about it.
    for( uint32_t b = 1; b <= 6; b++ ) {
	const Bound &bound = _geom.boundary( b );
	if( bound.type() == BOUND_NEUMANN && !bound.is_constant() )
	    changed = true;
    }

    for( uint32_t a = 0; a < n; a++ ) {
	uint32_t mesh = _geom.mesh( a );
	// A Neumann mask node is fixed, but its potential is never read.
	double fixed_epot = 0.0;
	if( (mesh & SMESH_NODE_FIXED) && !SMESH_NODE_IS_NEUMANN_MASK( mesh ) )
	    fixed_epot = (*_epot)( a );
	if( mesh != _key_mesh[a] || fixed_epot != _key_epot[a] ) {
	    _key_mesh[a] = mesh;
	    _key_epot[a] = fixed_epot;
	    changed = true;
	}
    }

    return( changed );
}


bool EpotMatrixSolver::preprocess( MeshScalarField &epot, const MeshScalarField &scharge )
{
    _epot = &epot;
    _scharge = &scharge;
    EpotSolver::preprocess( epot );

    // The linear system is kept from one solve to the next, see the header.
    const bool rebuild = ( linear_inputs_changed() || !_fd_mat );
    if( rebuild )
	build_linear();
    ibsimu.message( 1 ) << "dof = " << _dof << "\n";

    // Constant right hand side of this solve: boundary part plus space
    // charge. With the shield model the charge enters through
    // update_nonlinear_node() instead.
    const bool add_scharge = ( _plasma != PLASMA_SHIELD );
    _fd_vec_base = _fd_vec_bc;
    for( uint32_t a = 0; a < _dof; a++ ) {
	if( add_scharge && _row_scharge[a] )
	    _fd_vec_base(a) += -(*_scharge)(_dof2node[a])*_geom.h()*_geom.h()/EPSILON0;
    }

    return( rebuild );
}


/* Build the part of the system that depends neither on the solution nor on
 * the space charge: the node numbering, the matrix J0 and the boundary part
 * of the right hand side.
 */
void EpotMatrixSolver::build_linear( void )
{
    auto time_t0 = Clock::now();
    int nthreads = (int)ibsimu.get_thread_count();

    reset_matrix();

    // Node numbering: _n2d holds the mesh node index for a fixed node and
    // the row number for a free one.
    _n2d.resize( _geom.size() );
    _dof2node.clear();
    _dof2node.reserve( _geom.nodecount() );
    for( int32_t a = 0; a < (int32_t)_geom.nodecount(); a++ ) {
	if( _geom.mesh(a) & SMESH_NODE_FIXED ) {
	    _n2d(a) = N2D_TYPE_FIXED | a;
	} else {
	    _n2d(a) = N2D_TYPE_FREE | _dof;
	    _dof2node.push_back( a );
	    _dof++;
	}
    }
    if( _dof == 0 )
        throw( Error( ERROR_LOCATION, "zero degrees of freedom" ) );

    _fd_mat = new CRowMatrix( _dof, _dof );
    _fd_vec = new Vector( _dof );
    _fd_vec->clear();
    _d_vec.resize( _dof );
    _row_scharge.assign( _dof, 1 );

    // CRowMatrix::construct_add() appends at one shared cursor, so the
    // parallel loop below can not write to the matrix directly. set_link()
    // collects each row's entries here instead. Capacity is reserved for
    // the largest stencil (3D: 6 neighbours + centre) so that no row
    // reallocates while other threads are working.
    _row_entries.assign( _dof, std::vector<std::pair<int32_t,double> >() );
    for( uint32_t a = 0; a < _dof; a++ )
	_row_entries[a].reserve( 7 );

    // Each node is handled by one thread, which writes only that node's
    // _row_entries[a], _fd_vec(a) and _row_scharge[a]. A user-supplied
    // boundary value functor is called from here without locking, so it
    // must be thread-safe.
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		Vec3D x( _geom.origo(0) + _geom.h()*i,
			 _geom.origo(1) + _geom.h()*j,
			 _geom.origo(2) + _geom.h()*k );

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;

		if( mesh & SMESH_NODE_FIXED )
		    continue;
		else if( node_id == SMESH_NODE_ID_PURE_VACUUM )
		    add_vacuum_node( i, j, k, x );
		else if( node_id == SMESH_NODE_ID_NEAR_SOLID )
		    add_near_solid_node( i, j, k, x );
		else if( node_id == SMESH_NODE_ID_NEUMANN )
		    add_neumann_node( i, j, k, x );
	    }
	}
    }

    // Copy the rows into the matrix in ascending row order, as
    // construct_add() requires.
    //
    // Entries with the same column are merged first: a link that
    // set_link() redirects from a Neumann mask lands on a column the row
    // has already (the mirror neighbour or the diagonal).
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
    _fd_mat->order_ascending();
    std::vector<std::vector<std::pair<int32_t,double> > >().swap( _row_entries );

    _fd_vec_bc = *_fd_vec;

    // Linear diagonal J0(a,a) and its position in the value array, for
    // build_mat_vec() and get_resjac().
    _fd_mat_diag0.resize( _dof );
    _fd_mat_diag_idx.resize( _dof );
    for( uint32_t a = 0; a < _dof; a++ ) {
	int idx = _fd_mat->value_index( a, a );
	if( idx < 0 ) {
	    // Every stencil builder emits a centre link, so this row was
	    // never built. Do not keep the broken system.
	    reset_matrix();
	    throw( Error( ERROR_LOCATION, "no diagonal element on matrix row "
			  + to_string(a) ) );
	}
	_fd_mat_diag_idx[a] = idx;
	_fd_mat_diag0[a] = _fd_mat->value_ptr()[idx];
    }

    _time_linbuild += std::chrono::duration<double>( Clock::now()-time_t0 ).count();
}


/* Set _fd_mat to the linear matrix J0 and _fd_vec to the right hand side
 * at the current solution _sol. The solution dependent part of the jacobian
 * diagonal is returned in _d_vec.
 */
void EpotMatrixSolver::build_mat_vec( void )
{
    if( !_fd_mat )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    int nthreads = (int)ibsimu.get_thread_count();

    // Start from the linear system. get_resjac() leaves the jacobian
    // diagonal in the matrix, so that is restored too.
    *_fd_vec = _fd_vec_base;
    double *mval = _fd_mat->value_ptr();
    for( uint32_t a = 0; a < _dof; a++ )
	mval[_fd_mat_diag_idx[a]] = _fd_mat_diag0[a];
    _fd_mat->values_changed();

    // Plasma and compensation terms, the only part that depends on the
    // solution.
    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP ||
	_plasma == PLASMA_SHIELD || _comp.enabled ) {

	auto time_t0 = Clock::now();

	// Confinement depends on the current iterate.
	if( _comp.enabled && _comp.ascent )
	    update_comp_mask();

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
    // R = J0*X - B(X) and J = J0 - I*D(X). build_mat_vec() gives J0 in
    // _fd_mat, B(X) in _fd_vec and D(X) in _d_vec; _d_vec is cleared first
    // because rows without a plasma term are not written.
    _d_vec.clear();
    _sol = &X;
    build_mat_vec();

    Vector w = (*_fd_mat) * X;

    // The diagonal is written straight into the value array.
    double *mval = _fd_mat->value_ptr();
    for( uint32_t a = 0; a < _dof; a++ ) {
	(*_fd_vec)(a) = w(a) - (*_fd_vec)(a);
	mval[_fd_mat_diag_idx[a]] = _fd_mat_diag0[a] - _d_vec(a);
    }
    _fd_mat->values_changed();

    *J = _fd_mat;
    *R = _fd_vec;
}


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
