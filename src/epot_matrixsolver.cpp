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
    : EpotSolver(geom), _dof(0), _fd_mat(0), _fd_vec(0), _d_vec(0),
      _linear_built(false), _fd_vec_base(0)
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
    if( _plasma != PLASMA_PEXP && _plasma != PLASMA_NSIMP && _plasma != PLASMA_SHIELD )
	return;

    bool inplasma = true;
    if( _plasma_calc_func )
	// Test if within plasma calculation region
	inplasma = (*_plasma_calc_func)(x);

    if( !inplasma ) {
	(*_d_vec)(a) = 0; // No plasma calculation
    } else if( _plasma == PLASMA_PEXP ) {
	double p = (*_sol)(a);
	double rhst, drhst;
	pexp_newton( rhst, drhst, p );
	(*_fd_vec)(a) += rhst;
	(*_d_vec)(a) = drhst;
    } else if( _plasma == PLASMA_NSIMP ) {
	double p = (*_sol)(a);
	double rhst, drhst;
	nsimp_newton( rhst, drhst, p );
	(*_fd_vec)(a) += rhst;
	(*_d_vec)(a) = drhst;
    } else if( _plasma == PLASMA_SHIELD ) {
	double p = (*_sol)(a);
	double rhst, drhst;
	shield_newton( rhst, drhst, p );
	double R = -(*_scharge)(i,j,k)*_geom.h()*_geom.h()/EPSILON0;
	(*_fd_vec)(a) += R*rhst;
	(*_d_vec)(a) = R*drhst;
    }
}


/* See the doc comment on node_epsilon_r() in the header. */
double EpotMatrixSolver::node_epsilon_r( uint32_t mesh_value ) const
{
    uint32_t node_id = mesh_value & SMESH_NODE_ID_MASK;
    if( node_id != SMESH_NODE_ID_PURE_VACUUM && node_id != SMESH_NODE_ID_PURE_VACUUM_FIX )
	return( 1.0 ); // near-solid/Neumann/Dirichlet/fine-boundary: never inside a dielectric here
    uint32_t solid_number = mesh_value & SMESH_NEAR_SOLID_INDEX_MASK;
    if( solid_number == 0 )
	return( 1.0 ); // plain vacuum
    return( _geom.get_boundary( solid_number ).value() );
}


double EpotMatrixSolver::face_epsilon( double eps_a, double eps_b )
{
    return( 2.0*eps_a*eps_b/(eps_a+eps_b) );
}


/* See the doc comment on neighbor_epsilon_r() in the header. */
double EpotMatrixSolver::neighbor_epsilon_r( uint32_t mesh_value, double eps_self ) const
{
    uint32_t node_id = mesh_value & SMESH_NODE_ID_MASK;
    if( node_id == SMESH_NODE_ID_PURE_VACUUM || node_id == SMESH_NODE_ID_PURE_VACUUM_FIX )
	return( node_epsilon_r( mesh_value ) );
    if( node_id == SMESH_NODE_ID_DIRICHLET )
	return( eps_self ); // no material info available at a fixed node; assume self's own medium
    return( 1.0 ); // near-solid/Neumann/fine-boundary: vacuum, by construction here
}


void EpotMatrixSolver::add_vacuum_node( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    (void)x; // no longer used here -- see update_nonlinear_node()
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;

    // Relative permittivity of this node and, per direction, of its
    // neighbour. Both are 1.0 (vacuum) unless a dielectric solid is
    // involved (see node_epsilon_r()), in which case face_epsilon()'s
    // harmonic mean gives the standard finite-volume flux-matching
    // coefficient for that face. When every eps involved is 1.0 (no
    // dielectric anywhere near this node -- the common case, and the
    // only case before dielectric support existed) every face_epsilon()
    // call below reduces to exactly 1.0, reproducing the original
    // fixed 1.0/-2.0/-4.0/-6.0 stencil unchanged. This is only correct
    // for grid-aligned dielectric surfaces: the two adjacent mesh
    // nodes are assumed to sit exactly one full cell width h apart on
    // either side of a flat material interface (no sub-cell/near-solid
    // fractional-distance handling yet -- see class documentation).
    double eps_self = node_epsilon_r( _geom.mesh(i,j,k) );
    double cof = 0.0;

    switch( _geom.geom_mode() ) {
    case MODE_1D: {
        double wm = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i-1,j,k), eps_self) );
        double wp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i+1,j,k), eps_self) );
        set_link( a, _n2d(i-1,j,k), wm );
        set_link( a, _n2d(i+1,j,k), wp );
        cof = wm+wp;
        break; }
    case MODE_2D: {
        double wxm = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i-1,j,k), eps_self) );
        double wxp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i+1,j,k), eps_self) );
        double wym = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j-1,k), eps_self) );
        double wyp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j+1,k), eps_self) );
        set_link( a, _n2d(i,j-1,k), wym );
        set_link( a, _n2d(i-1,j,k), wxm );
        set_link( a, _n2d(i+1,j,k), wxp );
        set_link( a, _n2d(i,j+1,k), wyp );
        cof = wxm+wxp+wym+wyp;
        break; }
    case MODE_CYL: {
        double wxm = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i-1,j,k), eps_self) );
        double wxp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i+1,j,k), eps_self) );
        double wym = (1.0-0.5/j)*face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j-1,k), eps_self) );
        double wyp = (1.0+0.5/j)*face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j+1,k), eps_self) );
        set_link( a, _n2d(i,j-1,k), wym );
        set_link( a, _n2d(i-1,j,k), wxm );
        set_link( a, _n2d(i+1,j,k), wxp );
        set_link( a, _n2d(i,j+1,k), wyp );
        cof = wxm+wxp+wym+wyp;
        break; }
    case MODE_3D: {
        double wxm = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i-1,j,k), eps_self) );
        double wxp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i+1,j,k), eps_self) );
        double wym = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j-1,k), eps_self) );
        double wyp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j+1,k), eps_self) );
        double wzm = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j,k-1), eps_self) );
        double wzp = face_epsilon( eps_self, neighbor_epsilon_r(_geom.mesh(i,j,k+1), eps_self) );
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
    // entirely on the coefficients above, not the rhs -- and rho
    // (scharge) should be exactly zero inside a solid dielectric's
    // bulk in any case, since particles are absorbed at solid surfaces
    // rather than depositing charge past them.
    if( _plasma != PLASMA_SHIELD )
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

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i), -2.0/(beta*beta) );
	set_link( a, _n2d(i+1), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a, _n2d(i), -2.0/(alpha*alpha) );
	set_link( a, _n2d(i-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	set_link( a, _n2d(i), -2.0/(alpha*beta) );
	set_link( a, _n2d(i-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1), 2.0/((alpha+beta)*beta) );
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

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
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
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j+1), 2.0/((alpha+beta)*beta) );
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

    // Factors for X axis
    if( bindex & EPOT_SOLVER_BXMIN ) {
	cof += 2.0/(beta*beta);
	set_link( a, _n2d(i+1,j), 2.0/(beta*beta) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value(x) / beta;
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i-1,j), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j), 2.0/((alpha+beta)*beta) );
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
	// On-axis
	cof += 4.0;
	set_link( a, _n2d(i,j+1), 4.0 );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	cof += 2.0/(alpha*alpha);
	set_link( a, _n2d(i,j-1), 2.0/(alpha*alpha) );
	(*_fd_vec)(a) += (1.0)/(2.0*alpha)*(2.0/alpha+1.0/j)*
	    2.0*_geom.h()*_geom.get_boundary(4).value(x) / alpha;
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1), 1.0/(alpha+beta)*(2.0/alpha-1.0/j) );
	set_link( a, _n2d(i,j+1), 1.0/(alpha+beta)*(2.0/beta+1.0/j) );
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
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i-1,j,k), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i+1,j,k), 2.0/((alpha+beta)*beta) );
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
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j-1,k), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j+1,k), 2.0/((alpha+beta)*beta) );
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
    } else {
	cof += 2.0/(alpha*beta);
	set_link( a, _n2d(i,j,k-1), 2.0/((alpha+beta)*alpha) );
	set_link( a, _n2d(i,j,k+1), 2.0/((alpha+beta)*beta) );
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

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i), -2.0 );
	set_link( a, _n2d(i+1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else {
	set_link( a, _n2d(i-1), 2.0 );
	set_link( a, _n2d(i), -2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    }
}


void EpotMatrixSolver::add_neumann_node_2d( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i+1,j), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a, _n2d(i-1,j), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    } else {
	set_link( a, _n2d(i-1,j), 1.0 );
	set_link( a, _n2d(i+1,j), 1.0 );
    }
    
    if( bindex & EPOT_SOLVER_BYMIN ) {
	set_link( a, _n2d(i,j+1), 2.0 );
	(*_fd_vec)(a) += -2.0*_geom.h()*_geom.get_boundary(3).value( x );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	set_link( a, _n2d(i,j-1), 2.0 );
	(*_fd_vec)(a) += -2.0*_geom.h()*_geom.get_boundary(4).value( x );
    } else {
	set_link( a, _n2d(i,j-1), 1.0 );
	set_link( a, _n2d(i,j+1), 1.0 );
    }

    set_link( a, _n2d(i,j), -4.0 );
}


void EpotMatrixSolver::add_neumann_node_cyl( uint32_t i, uint32_t j, const Vec3D &x )
{
    uint32_t a = _n2d(i,j) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j);

    if( bindex & EPOT_SOLVER_BYMIN ) {
	// On-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    set_link( a, _n2d(i+1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    set_link( a, _n2d(i-1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i-1,j), 1.0 );
	    set_link( a, _n2d(i+1,j), 1.0 );
	}
	set_link( a, _n2d(i,j+1), 4.0 );
	set_link( a, _n2d(i,j), -6.0 );
    } else {
	// Off-axis
	if( bindex & EPOT_SOLVER_BXMIN ) {
	    set_link( a, _n2d(i+1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
	} else if( bindex & EPOT_SOLVER_BXMAX ) {
	    set_link( a, _n2d(i-1,j), 2.0 );
	    (*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i-1,j), 1.0 );
	    set_link( a, _n2d(i+1,j), 1.0 );
	}
	
	if( bindex & EPOT_SOLVER_BYMAX ) {
	    set_link( a, _n2d(i,j-1), 2.0 );
	    (*_fd_vec)(a) += (1.0+0.5/j)*2.0*_geom.h()*_geom.get_boundary(2).value( x );
	} else {
	    set_link( a, _n2d(i,j-1), 1.0-0.5/j );
	    set_link( a, _n2d(i,j+1), 1.0+0.5/j );
	}
	
	set_link( a, _n2d(i,j), -4.0 );
    }
}


void EpotMatrixSolver::add_neumann_node_3d( uint32_t i, uint32_t j, uint32_t k, const Vec3D &x )
{
    uint32_t a = _n2d(i,j,k) & N2D_INDEX_MASK;
    uint8_t bindex = boundary_index(i,j,k);

    if( bindex & EPOT_SOLVER_BXMIN ) {
	set_link( a, _n2d(i+1,j,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(1).value( x );
    } else if( bindex & EPOT_SOLVER_BXMAX ) {
	set_link( a,_n2d(i-1,j,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(2).value( x );
    } else {
	set_link( a, _n2d(i-1,j,k), 1.0 );
	set_link( a, _n2d(i+1,j,k), 1.0 );
    }

    if( bindex & EPOT_SOLVER_BYMIN ) {
	set_link( a, _n2d(i,j+1,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(3).value( x );
    } else if( bindex & EPOT_SOLVER_BYMAX ) {
	set_link( a, _n2d(i,j-1,k), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(4).value( x );
    } else {
	set_link( a, _n2d(i,j-1,k), 1.0 );
	set_link( a, _n2d(i,j+1,k), 1.0 );
    }
    
    if( bindex & EPOT_SOLVER_BZMIN ) {
	set_link( a, _n2d(i,j,k+1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(5).value( x );
    } else if( bindex & EPOT_SOLVER_BZMAX ) {
	set_link( a, _n2d(i,j,k-1), 2.0 );
	(*_fd_vec)(a) += 2.0*_geom.h()*_geom.get_boundary(6).value( x );
    } else {
	set_link( a, _n2d(i,j,k-1), 1.0 );
	set_link( a, _n2d(i,j,k+1), 1.0 );
    }
    
    set_link( a, _n2d(i,j,k), -6.0 );
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
    _n2d.resize( _geom.size() );
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
	    _dof++;
	}
    }

    ibsimu.message( 1 ) << "dof = " << _dof << "\n";
    if( _dof == 0 )
        throw( Error( ERROR_LOCATION, "zero degrees of freedom" ) );
    
    // Allocate problem matrix and vector
    _fd_mat = new CRowMatrix( _dof, _dof );
    _fd_vec = new Vector( _dof );
}


void EpotMatrixSolver::build_mat_vec( void )
{
    if( _dof == 0  )
	throw( Error( ERROR_LOCATION, "preprocess not done" ) );

    int nthreads = (int)ibsimu.get_thread_count();

    if( !_linear_built ) {

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
	for( uint32_t a = 0; a < _dof; a++ ) {
	    const std::vector<std::pair<int32_t,double> > &row = _row_entries[a];
	    for( size_t n = 0; n < row.size(); n++ )
		_fd_mat->construct_add( a, row[n].first, row[n].second );
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
	_fd_mat_diag0.resize( _dof );
	for( uint32_t a = 0; a < _dof; a++ )
	    _fd_mat_diag0[a] = _fd_mat->get( a, a );

	_linear_built = true;

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
    if( _plasma == PLASMA_PEXP || _plasma == PLASMA_NSIMP || _plasma == PLASMA_SHIELD ) {
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
    _d_vec = new Vector( _dof );

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
    for( uint32_t a = 0; a < _dof; a++ )
	_fd_mat->set(a,a) = _fd_mat_diag0[a];

    // Linear part
    Vector w = (*_fd_mat) * X;

    for( uint32_t a = 0; a < _dof; a++ ) {

	(*_fd_vec)(a) = w(a) - (*_fd_vec)(a);
	// Set (not subtract from) the diagonal: it must be re-derived
	// fresh from the pure-linear J0(a,a) each call, not adjusted
	// relative to whatever a previous iteration left there.
	_fd_mat->set(a,a) = _fd_mat_diag0[a] - (*_d_vec)(a);
    }

    *J = _fd_mat;
    *R = _fd_vec;

    delete _d_vec;
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
