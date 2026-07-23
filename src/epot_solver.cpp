/*! \file epot_solver.cpp
 *  \brief Poisson equation problem for solving electric potential.
 */

/* Copyright (c) 2005-2013,2017,2021 Taneli Kalvas. All rights reserved.
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

#include "epot_solver.hpp"
#include "timer.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "compmath.hpp"
#include "constants.hpp"
#include <omp.h>
#include <atomic>



/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */


EpotSolver::EpotSolver( Geometry &geom ) 
    : _geom(geom), _plasma(PLASMA_NONE), 
      _rhoe(0.0), _Te(0.0), _Up(0.0), 
      _force_pot(0.0), _force_pot_func(0), _force_pot_func2(0), 
      _init_plasma_func(0), _plasma_calc_func(0)
{

}


EpotSolver::EpotSolver( const EpotSolver &epsolver, Geometry &geom )
    : _geom(geom), _plasma(epsolver._plasma),
      _rhoe(epsolver._rhoe), _Te(epsolver._Te), _Up(epsolver._Up),
      _rhoi(epsolver._rhoi), _Ei(epsolver._Ei), 
      _force_pot(epsolver._force_pot), _force_pot_func(epsolver._force_pot_func),
      _force_pot_func2(epsolver._force_pot_func2),
      _init_plasma_func(epsolver._init_plasma_func),
      _plasma_calc_func(epsolver._plasma_calc_func)
{

}


EpotSolver::EpotSolver( Geometry &geom, std::istream &s )
    : _geom(geom)
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


/* ************************************** *
 * EpotSolver constructing               *
 * ************************************** */


void EpotSolver::set_parameters( const EpotSolver &epsolver )
{
    _plasma = epsolver._plasma;
    _rhoe = epsolver._rhoe;
    _Te = epsolver._Te;
    _Up = epsolver._Up;
    _rhoi = epsolver._rhoi;
    _Ei = epsolver._Ei;
    _force_pot = epsolver._force_pot;
    _force_pot_func = epsolver._force_pot_func;
    _force_pot_func2 = epsolver._force_pot_func2;
    _init_plasma_func = epsolver._init_plasma_func;
    _plasma_calc_func = epsolver._plasma_calc_func;
}


void EpotSolver::set_forced_potential_volume( double force_pot, 
					      CallbackFunctorB_V *force_pot_func )
{
    _force_pot = force_pot;
    _force_pot_func = force_pot_func;
}


void EpotSolver::set_forced_potential_volume( CallbackFunctorD_V *force_pot_func )
{
    _force_pot_func2 = force_pot_func;
}


namespace {
    /* A dielectric solid's interior is tagged the same way as plain
     * vacuum (SMESH_NODE_ID_PURE_VACUUM[_FIX], material number in the
     * lower bits -- see Geometry::build_mesh_parallel_thread_3d()), so
     * the raw tag is enough to recognise it *unless* it reaches a
     * simulation box wall: Geometry::override_dielectric_box_boundary_3d()
     * resets such a node and lets it be reclassified by the ordinary
     * box-edge decision tree, which retags it SMESH_NODE_ID_NEUMANN or
     * SMESH_NODE_ID_DIRICHLET -- discarding the material number, even
     * though the node is still physically dielectric material. That's
     * the right call for the matrix discretisation itself (the box's
     * own condition, not the dielectric's bulk equation, applies right
     * at that surface -- see the doc comment on that function), but it
     * means the raw tag alone can no longer answer "is this dielectric"
     * for a Neumann-tagged node. This recovers the answer geometrically
     * when the fast path (raw tag) can't, mirroring
     * EpotMatrixSolver::self_material_at()'s slow path -- duplicated
     * rather than shared since this class doesn't depend on that one
     * (which is the derived class).
     */
    bool is_dielectric_at( const Geometry &geom, uint32_t mesh_value, const Vec3D &x )
    {
	uint32_t node_id = mesh_value & SMESH_NODE_ID_MASK;
	if( ( node_id == SMESH_NODE_ID_PURE_VACUUM || node_id == SMESH_NODE_ID_PURE_VACUUM_FIX ) &&
	    ( mesh_value & SMESH_NEAR_SOLID_INDEX_MASK ) >= 7 )
	    return( true ); // fast path -- tag already says dielectric-interior directly

	uint32_t solid_number = geom.inside( x );
	if( solid_number < 7 )
	    return( false ); // plain vacuum (or, shouldn't happen, outside the box)
	return( geom.get_boundary( solid_number ).type() == BOUND_DIELECTRIC );
    }
}


void EpotSolver::set_initial_plasma( double Up,
				     CallbackFunctorB_V *init_plasma_func )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = Up;
    if( !init_plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = init_plasma_func;
    reset_problem();
}


/* See the doc comment in epot_solver.hpp. Deliberately does not touch
 * _plasma/_init_plasma_func/reset_problem() at all -- unlike
 * set_initial_plasma(), this doesn't define a plasma model by itself,
 * it only fills in a starting guess for scharge that the caller then
 * passes to solve() alongside a real set_pexp_plasma()/
 * set_nsimp_plasma() call.
 *
 * Solid nodes are always skipped, regardless of what plasma_region_func
 * says: a Dirichlet (conductor) node is eliminated from the matrix
 * entirely, so scharge there is simply never read, but a dielectric
 * node is an ordinary *free* row (see is_dielectric_at()'s doc comment)
 * whose right-hand side directly includes -scharge*h*h/EPSILON0 (see
 * add_vacuum_node() and friends) -- setting rho there injects a large,
 * entirely spurious free-charge source term straight into the
 * dielectric's own Poisson equation, which has no mobile charge at all.
 * This is exactly the same category of mistake set_initial_plasma()'s
 * fixed-potential path had (see is_dielectric_at() and its call sites
 * in preprocess()), just landing on scharge here instead of a forced
 * potential. plasma_region_func only needs to describe the plasma's
 * geometric extent (e.g. "z < 0"); it does not need to know where the
 * solids are.
 */
void EpotSolver::set_initial_plasma_rho( MeshScalarField &scharge, double rho,
					 CallbackFunctorB_V *plasma_region_func ) const
{
    if( !plasma_region_func )
	throw( Error( ERROR_LOCATION, "NULL plasma region function" ) );

    int nthreads = (int)ibsimu.get_thread_count();
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		Vec3D x( _geom.origo(0) + _geom.h()*i,
			 _geom.origo(1) + _geom.h()*j,
			 _geom.origo(2) + _geom.h()*k );

		if( !(*plasma_region_func)( x ) )
		    continue;

		uint32_t mesh = _geom.mesh(i,j,k);
		if( (mesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET )
		    continue; // conductor -- eliminated row, never read
		if( is_dielectric_at( _geom, mesh, x ) )
		    continue; // real dielectric material -- no mobile charge, see above

		scharge(i,j,k) = rho;
	    }
	}
    }
}


void EpotSolver::set_plasma_calc_region( CallbackFunctorB_V *plasma_calc_func )
{
    _plasma_calc_func = plasma_calc_func;
    reset_problem();
}


void EpotSolver::set_pexp_plasma( double rhoe, double Te, double Up )
{
    _init_plasma_func = NULL;
    _plasma     = PLASMA_PEXP;
    _rhoe       = -fabs(rhoe); // Ensure correct sign of charge density
    _Te         = Te;
    _Up         = Up;
    reset_problem();
}



void EpotSolver::set_nsimp_initial_plasma( CallbackFunctorB_V *init_plasma_func )
{
    _plasma     = PLASMA_INITIAL;
    _Up         = 0.0;
    if( !init_plasma_func )
	throw( Error( ERROR_LOCATION, "NULL initial plasma function" ) );
    _init_plasma_func = init_plasma_func;
    reset_problem();
}


void EpotSolver::set_nsimp_plasma( double rhop, double Ep, 
				   std::vector<double> rhoi, std::vector<double> Ei )
{
    _init_plasma_func = NULL;
    _plasma = PLASMA_NSIMP;
    _rhoi.clear();
    _Ei.clear();

    if( rhop < 0 )
	ibsimu.message( MSG_WARNING, 1 ) << "Warning: negative fast species charge density\n";
    _rhoi.push_back( rhop );
    _Ei.push_back( Ep );

    if(  rhoi.size() != Ei.size() )
	throw( Error( ERROR_LOCATION, "different size rhoi and Ep vectors" ) );
    for( size_t a = 0; a < rhoi.size(); a++ ) {
	if( rhoi[a] < 0 )
	    ibsimu.message( MSG_WARNING, 1 ) << "Warning: negative thermal species " << a << " charge density\n";
	_rhoi.push_back( rhoi[a] );
	_Ei.push_back( Ei[a] );
    }

    reset_problem();
}


void EpotSolver::set_shield_plasma( double Tm, double Um )
{
    _init_plasma_func = NULL;
    _plasma     = PLASMA_SHIELD;
    _Te        = Tm;
    _Up        = Um;
    reset_problem();
}


void EpotSolver::shield_newton( double &rhs, double &drhs, double epot ) const
{
    double K = tanh(_plA*epot - _plB);
    rhs  = 0.5*(1.0+K);
    drhs = 0.5*_plA*(1.0-K*K);
}


void EpotSolver::pexp_newton( double &rhs, double &drhs, double epot ) const
{
    rhs  = _plA*exp( _plB*epot - _plC );
    drhs = _plB*rhs;
}


void EpotSolver::nsimp_newton( double &rhs, double &drhs, double epot ) const
{
    // fast
    double r = _plB*epot;
    if( r > 0.0 ) {
	rhs  = _plA*( 1.0 + 2.0/sqrt(M_PI)*r );
	drhs = _plC;
    } else {
	rhs  = _plA*( 1.0 + erf( r ) );
	drhs = _plC*exp( -r*r );
    }

    // thermal
    for( size_t i = 0; i < _plD.size(); i++ ) {
	double k = _plE[i]*epot;
	if( k > 10.0 ) {
	    // Prevent overruns with exponential function during iteration
	    // by changing to linear slope
	    double exp10 = exp( 10.0 );
	    double d = _plD[i]*exp10;
	    double w = _plD[i]*exp10 + d*(k - 10.0);
	    rhs  += w;
	    drhs += _plE[i]*d;
	} else {
	    double w = _plD[i]*exp( k );
	    rhs  += w;
	    drhs += _plE[i]*w;
	}
    }
}


void EpotSolver::preprocess( MeshScalarField &epot )
{
    if( _plasma == PLASMA_PEXP ) {
	// Calculate plasma parameters for positive ion extraction. 
	// Calculation done as rhs + A*exp(B*x-C)
	_plA = -_rhoe*_geom.h()*_geom.h()/EPSILON0;
	_plB = 1.0/_Te;
	_plC = _Up/_Te;
    } else if( _plasma == PLASMA_NSIMP ) {
	// Calculate plasma parameters for negative ion extraction. 
	// The total right hand side is rhs + r_fast + r_thermal, where
	// r_fast = A(1+erf(B*x)) and
	// r_thermal = D*exp(E*x), for i >= 0. The derivatives are
	// r_fast' = A*B*2/sqrt(pi)*exp(-B^2*x^2)
	// r_thermal' = D*E*exp(E*x)
	_plA = -_rhoi[0]*_geom.h()*_geom.h()/EPSILON0;
	_plB = -1.0/_Ei[0];
	_plC = _plA*_plB*2.0/sqrt(M_PI);
	_plD.clear();
	_plE.clear();
	for( uint32_t i = 1; i < _rhoi.size(); i++ ) {
	    _plD.push_back( -_rhoi[i]*_geom.h()*_geom.h()/EPSILON0 );
	    _plE.push_back( -1.0/_Ei[i] );
	}
    } else if( _plasma == PLASMA_SHIELD ) {
	_plA = 1.0/_Te;
	_plB = _Up/_Te;
    }


    // Set forced vacuum nodes and dirichlet nodes to correct
    // potential. Mark fixed vacuum nodes with a tag. Every node only
    // reads/writes its own (i,j,k) slot of _geom.mesh()/epot(), with no
    // dependency on any other node, so this is embarrassingly parallel.
    // The one requirement this imposes: _force_pot_func/_force_pot_func2/
    // _init_plasma_func (user-supplied callbacks) and
    // Boundary::value() must be safe to call concurrently from multiple
    // threads, since this loop calls them with no locking.
    // node_id should always be one of the three cases handled below --
    // the final "else" is a defensive, should-never-happen case. It used
    // to just throw directly, which is fine in a serial loop but is
    // undefined behavior (typically std::terminate) if it escapes an
    // OpenMP parallel region uncaught. So it's turned into a shared
    // flag instead, checked (and thrown from, serially) after the
    // parallel loop ends.
    std::atomic<bool> bad_node( false );
    int nthreads = (int)ibsimu.get_thread_count();
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		Vec3D x( _geom.origo(0) + _geom.h()*i,
			 _geom.origo(1) + _geom.h()*j,
			 _geom.origo(2) + _geom.h()*k );

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;

		if( node_id == SMESH_NODE_ID_NEAR_SOLID ||
		    node_id == SMESH_NODE_ID_PURE_VACUUM ) {

		    // Vacuum -- but this can *also* be genuinely dielectric
		    // material, in two different ways: directly (raw
		    // SMESH_NODE_ID_PURE_VACUUM[_FIX] tag with a solid
		    // number >=7 in the lower bits -- see
		    // Geometry::build_mesh_parallel_thread_3d() -- because
		    // it satisfies the same homogeneous Laplace stencil
		    // shape as vacuum, just with a position-dependent
		    // permittivity applied by add_vacuum_node()/
		    // vacuum_face_coefficient()), or more subtly, tagged
		    // SMESH_NODE_ID_NEAR_SOLID instead: a dielectric that
		    // reaches a simulation box edge gets reset and
		    // reclassified by the ordinary box-edge decision tree
		    // (Geometry::override_dielectric_box_boundary_3d()), and
		    // if a *real conductor* also happens to be adjacent to
		    // that same node along a different axis, that tree
		    // tags it NEAR_SOLID -- discarding the material number,
		    // so it looks exactly like an ordinary near-conductor
		    // vacuum node by raw tag alone. Either way it is
		    // emphatically not a region a plasma or a forced
		    // potential can occupy: unlike real vacuum, it is a
		    // fixed physical medium and must always stay a
		    // genuinely free node. is_dielectric_at() (see its doc
		    // comment above) catches both cases -- fast raw-tag
		    // path for the first, geometric fallback for the
		    // second -- and is deliberately checked *after* the
		    // cheaper force_pot_func/force_pot_func2/init_plasma_func
		    // predicates below, so the (potentially expensive,
		    // STL-backed) geometry query only runs on nodes where
		    // forcing would otherwise actually apply. Missing the
		    // NEAR_SOLID case here is exactly what let a dielectric
		    // squeezed between a box wall and a grounded electrode
		    // stay pinned to the plasma potential right at that
		    // corner, even after the box-edge Neumann case was
		    // fixed.
		    double val;
		    if( _force_pot_func2 &&
			comp_isfinite( (val = (*_force_pot_func2)( x ))) &&
			!is_dielectric_at( _geom, mesh, x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = val;

		    } else if( _force_pot_func && (*_force_pot_func)( x ) &&
			       !is_dielectric_at( _geom, mesh, x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = _force_pot;

		    } else if( (_plasma == PLASMA_PEXP_INITIAL ||
				_plasma == PLASMA_NSIMP_INITIAL) &&
			       _init_plasma_func && (*_init_plasma_func)( x ) &&
			       !is_dielectric_at( _geom, mesh, x ) ) {

			// Mark as fixed vacuum
			_geom.mesh(i,j,k) |= SMESH_NODE_FIXED;
			epot(i,j,k) = _Up;

		    }

		} else if( node_id == SMESH_NODE_ID_NEUMANN ) {

		    // A dielectric reaching this Neumann box wall had its
		    // material number discarded when its tag was
		    // reclassified (see the doc comment on
		    // is_dielectric_at() above) -- it must still be
		    // excluded from forcing, exactly like the dielectric-
		    // interior case above, even though the raw tag alone
		    // can't tell us that anymore. Without this check, a
		    // dielectric slab reaching the mirror/Neumann side
		    // walls would get pinned to the plasma's initial
		    // potential right at that edge (while the rest of the
		    // slab solved correctly), instead of behaving like an
		    // ordinary permittivity-weighted mirror node the way
		    // add_neumann_node() already treats it. is_dielectric_at()
		    // does a real geometric query (Geometry::inside()), so
		    // it's deliberately checked *after* the cheaper
		    // force_pot_func/init_plasma_func predicates, not on
		    // every Neumann node unconditionally -- only nodes where
		    // forcing would otherwise actually apply pay for it.
		    if( _force_pot_func && (*_force_pot_func)( x ) ) {

			if( !is_dielectric_at( _geom, mesh, x ) ) {
			    // Mark as fixed vacuum
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM_FIX;
			    epot(i,j,k) = _force_pot;
			}

		    } else if( (_plasma == PLASMA_PEXP_INITIAL ||
				_plasma == PLASMA_NSIMP_INITIAL) &&
			       _init_plasma_func && (*_init_plasma_func)( x ) ) {

			if( !is_dielectric_at( _geom, mesh, x ) ) {
			    // Mark as fixed vacuum
			    _geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM_FIX;
			    epot(i,j,k) = _Up;
			}
		    }

		} else if( node_id == SMESH_NODE_ID_DIRICHLET ) {
		    
		    // Dirichlet
		    uint32_t boundary = mesh & SMESH_BOUNDARY_NUMBER_MASK;
		    epot(i,j,k) = _geom.get_boundary( boundary ).value( x );

		} else {

		    bad_node = true;

		}
	    }
	}
    }

    if( bad_node )
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


void EpotSolver::postprocess( void )
{
    // Remove fixed vacuum tags. Every node only touches its own
    // _geom.mesh(i,j,k) slot (pure bitmask logic, no callbacks, no
    // cross-node reads), so this is a safe, no-caveat parallel loop --
    // the gain is smaller than preprocess()'s since there's no callback
    // work here, but it's essentially free.
    int nthreads = (int)ibsimu.get_thread_count();
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {

		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_NEAR_SOLID_FIX ) {
		    // Change to near solid node, keeping index pointer
		    uint32_t index = SMESH_NEAR_SOLID_INDEX_MASK & mesh;
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | index;
		} else if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX ) {
		    // Change to vacuum node if not on boundary
		    if( i == 0 || i == _geom.size(0)-1 )
			continue;
		    if( (_geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
			 _geom.geom_mode() == MODE_3D) && (j == 0 || j == _geom.size(1)-1) )
			continue;
		    if( _geom.geom_mode() == MODE_3D && (k == 0 || k == _geom.size(2)-1) )
			continue;
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM;
		}
	    }
	}
    }

    // Change PURE_VACUUM_FIX nodes on Neumann boundaries back to
    // Neumann nodes.

    // Xmin and Xmax
    for( uint32_t bound = 1; bound <= 2; bound++ ) {
	uint32_t i = 0;
	if( bound == 2 ) i = _geom.size(0)-1;
	for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		uint32_t mesh = _geom.mesh(i,j,k);
		uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
		    _geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
		    _geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL ||
	_geom.geom_mode() == MODE_3D ) {
	// Ymin and Ymax
	for( uint32_t bound = 3; bound <= 4; bound++ ) {
	    uint32_t j = 0;
	    if( bound == 4 ) j = _geom.size(1)-1;
	    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
			_geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		    }
		}
	    }
	}
    }
    if( _geom.geom_mode() == MODE_3D ) {
	// Zmin and Zmax
	for( uint32_t bound = 5; bound <= 6; bound++ ) {
	    uint32_t k = 0;
	    if( bound == 6 ) k = _geom.size(2)-1;
	    for( uint32_t j = 0; j < _geom.size(1); j++ ) {
		for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		    uint32_t mesh = _geom.mesh(i,j,k);
		    uint32_t node_id = mesh & SMESH_NODE_ID_MASK;
		    if( node_id == SMESH_NODE_ID_PURE_VACUUM_FIX &&
			_geom.get_boundary(bound).type() == BOUND_NEUMANN ) {
			_geom.mesh(i,j,k) = SMESH_NODE_ID_NEUMANN | bound;
		    }
		}
	    }
	}
    }
}


bool EpotSolver::linear( void ) const
{
    switch( _plasma ) {
    case PLASMA_NONE:
    case PLASMA_PEXP_INITIAL:
    case PLASMA_NSIMP_INITIAL:
	return( true );
	break;
    case PLASMA_PEXP:
    case PLASMA_NSIMP:
    case PLASMA_SHIELD:
	return( false );
	break;
    }

    throw( ErrorAssert( ERROR_LOCATION ) );
}


uint8_t EpotSolver::boundary_index( uint32_t i ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index( uint32_t i, uint32_t j ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    if( k == 0 )
	ret += EPOT_SOLVER_BZMIN;
    else if( k == _geom.size(2)-1 )
	ret += EPOT_SOLVER_BZMAX;
    return( ret );
}


uint8_t EpotSolver::boundary_index_general( uint32_t i, uint32_t j, uint32_t k ) const
{
    uint8_t ret = 0;
    if( i == 0 )
	ret += EPOT_SOLVER_BXMIN;
    else if( i == _geom.size(0)-1 )
	ret += EPOT_SOLVER_BXMAX;
    if( _geom.geom_mode() == MODE_1D )
	return( ret );

    if( j == 0 )
	ret += EPOT_SOLVER_BYMIN;
    else if( j == _geom.size(1)-1 )
	ret += EPOT_SOLVER_BYMAX;
    if( _geom.geom_mode() == MODE_2D || _geom.geom_mode() == MODE_CYL )
	return( ret );

    if( k == 0 )
	ret += EPOT_SOLVER_BZMIN;
    else if( k == _geom.size(2)-1 )
	ret += EPOT_SOLVER_BZMAX;
    return( ret );
}


/* ************************************** *
 * Solver                                 *
 * ************************************** */


MeshScalarField *EpotSolver::evaluate_scharge( const ScalarField &__scharge ) const
{
    // If scharge not defined in same points as geometry, evaluate scharge at nodes.
    // Every node is independent (writes only to its own (i,j,k) slot, no
    // read of any neighbour), so this is embarrassingly parallel -- the
    // only requirement is that __scharge's operator() is safe to call
    // concurrently from multiple threads, since it gets no locking here.
    MeshScalarField *scharge = new MeshScalarField( _geom );
    int nthreads = (int)ibsimu.get_thread_count();
#pragma omp parallel for num_threads(nthreads) collapse(3) schedule(dynamic,64)
    for( uint32_t k = 0; k < _geom.size(2); k++ ) {
	for( uint32_t j = 0; j < _geom.size(1); j++ ) {
	    for( uint32_t i = 0; i < _geom.size(0); i++ ) {
		double z = k*_geom.h()+_geom.origo(2);
		double y = j*_geom.h()+_geom.origo(1);
		double x = i*_geom.h()+_geom.origo(0);
		(*scharge)(i,j,k) = __scharge( Vec3D(x,y,z) );
	    }
	}
    }

    return( scharge );
}


void EpotSolver::solve( MeshScalarField &epot, const ScalarField &__scharge )
{
    Timer t;

    ibsimu.message( 1 ) << "Solving problem\n";
    ibsimu.inc_indent();

    // Check if geometry is built
    if( !_geom.built() )
	throw( Error( ERROR_LOCATION, "geomerty not built" ) );

    // Set scharge to be used by solver
    bool scharge_internal = false;
    const MeshScalarField *scharge = dynamic_cast<const MeshScalarField *>( &__scharge );
    if( scharge == 0 || *scharge != _geom ) {
	if( *scharge != _geom )
	    ibsimu.message( 1 ) << "Converting space charge density to match geometry.\n";
	scharge_internal = true;
	scharge = evaluate_scharge( __scharge );
    }

    // Resize epot if necessary, otherwise leave old solution as a starting point
    if( epot != _geom )
	epot.reset( _geom.geom_mode(), _geom.size(), _geom.origo(), _geom.h() );

    // Solve problem
    subsolve( epot, *scharge );

    // Free scharge if allocated
    if( scharge_internal )
	delete scharge;

    // End timer
    t.stop();

    ibsimu.message( 1 ) << "time used = " << t << "\n";
    ibsimu.dec_indent();
}


/* ************************************** *
 * Misc                                   *
 * ************************************** */


const Geometry &EpotSolver::geometry( void ) const
{
    return( _geom );
}


void EpotSolver::debug_print( std::ostream &os ) const 
{
    os << "**EpotSolver\n";
    os << "plasma = " << _plasma << "\n";
    os << "rhoe = " << _rhoe << "\n";
    os << "Tc = " << _Te << "\n";
    os << "Up = " << _Up << "\n";
}

