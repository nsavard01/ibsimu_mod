/*! \file particleiterator.hpp
 *  \brief %Particle iterator
 */

/* Copyright (c) 2005-2013,2018,2022 Taneli Kalvas, Tobin Jones. All rights reserved.
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

#ifndef PARTICLEITERATOR_HPP
#define PARTICLEITERATOR_HPP 1


#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include "geometry.hpp"
#include "mat3d.hpp"
#include "compmath.hpp"
#include "trajectory.hpp"
#include "particles.hpp"
#include "vectorfield.hpp"
#include "meshscalarfield.hpp"
#include "scharge.hpp"
#include "scheduler.hpp"
#include "polysolver.hpp"
#include "particledatabase.hpp"
#include "cfifo.hpp"
#include "ibsimu.hpp"


//#define DEBUG_PARTICLE_ITERATOR 1


#ifdef DEBUG_PARTICLE_ITERATOR 
#define DEBUG_MESSAGE(x) ibsimu.message(MSG_DEBUG_GENERAL,1) << x
#define DEBUG_INC_INDENT() ibsimu.inc_indent()
#define DEBUG_DEC_INDENT() ibsimu.dec_indent()
#else
#define DEBUG_MESSAGE(x) do {} while(0)
#define DEBUG_INC_INDENT() do {} while(0)
#define DEBUG_DEC_INDENT() do {} while(0)
#endif


#define COLLISION_EPS 1.0e-6


/*! \brief %Particle iterator type.
 */
enum particle_iterator_type_e {
    PARTICLE_ITERATOR_ADAPTIVE = 0,
    PARTICLE_ITERATOR_FIXED_STEP_LEN
};



/*! \brief %Mesh intersection (collision) coordinate data
 *
 * Contains one coordinate data and the direction of particle for one
 * intersection between particle trajectory and mesh plane. Templated
 * for particle point type (see ParticlePBase).
 *
 * Templated for particle point type (see ParticlePBase).
 */
template <class PP> class ColData {
public:
    PP                _x;         /*!< \brief %Mesh intersection coordinates. */
    int               _dir;       /*!< \brief Direction of particle at intersection.
				   *  i: -1/+1, j: -2/+2, k: -3:/+3. */
    
    /*! \brief Default constructor.
     */
    ColData() : _dir(0) {}

    /*! \brief Constructor for collision at \a x into direction \a dir.
     */
    ColData( PP x, int dir ) : _x(x), _dir(dir) {}
    
    /*! \brief Compare coldata entry times.
     *
     *  Used for sorting coldata entries.
     */
    bool operator<( const ColData &cd ) const {
	return( _x[0] < cd._x[0] );
    }

    /*! \brief Find mesh intersections of linearly interpolated
     *  particle trajectory segment.
     *
     *  Makes a linear interpolation between points \a x1 and \a x2
     *  and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_linear( std::vector<ColData> &coldata, const Mesh &mesh,
				      const PP &x1, const PP &x2 ) {
	
	DEBUG_MESSAGE( "Building coldata using linear interpolation\n" );
	DEBUG_INC_INDENT();

	coldata.clear();

	for( size_t a = 0; a < PP::dim(); a++ ) {
	    
            int a1 = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
            int a2 = (int)floor( (x2[2*a+1]-mesh.origo(a))/mesh.h() );
            if( a1 > a2 ) {
                int a = a2;
                a2 = a1;
                a1 = a;
            }

            for( int b = a1+1; b <= a2; b++ ) {

                // Save intersection coordinates
                double K = (b*mesh.h() + mesh.origo(a) - x1[2*a+1]) / 
                    (x2[2*a+1] - x1[2*a+1]);
                if( K < 0.0 ) K = 0.0;
                else if( K > 1.0 ) K = 1.0;

		DEBUG_MESSAGE( "Adding point " << x1 + (x2-x1)*K << "\n" );

                if( x2[2*a+1] > x1[2*a+1] )
                    coldata.push_back( ColData( x1 + (x2-x1)*K, a+1 ) );
                else
                    coldata.push_back( ColData( x1 + (x2-x1)*K, -a-1 ) );
            }
        }

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );

	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata built\n" );
    }

    /*! \brief Find mesh intersections of polynomially interpolated
        particle trajectory segment.
     *
     *  Makes a polynomial interpolation between points \ə x1 and \a
     *  x2 and searches intersection points of this line and \a
     *  mesh. Intersection points are saved to vector \a coldata in
     *  increasing time order.
     */
    static void build_coldata_poly( std::vector<ColData> &coldata, const Mesh &mesh,
				    const PP &x1, const PP &x2 ) {
	
	DEBUG_MESSAGE( "Building coldata using polynomial interpolation\n" );
	DEBUG_INC_INDENT();

	coldata.clear();

	// Construct trajectory representation.
	//
	// Stack array rather than new[]/delete[]: this function runs once
	// per accepted ODE step of every particle -- the innermost loop of
	// the whole mover -- so a heap allocation here is a malloc/free
	// pair per step, millions of them over a run, for what is at most
	// three tiny PODs. PP::dim() is 2 (ParticleP2D/ParticlePCyl) or 3
	// (ParticleP3D), fixed per particle type, so 3 always suffices;
	// the guard below keeps that honest if a higher-dimensional
	// particle type is ever added. Also removes the leak that the old
	// new[] ... delete[] pair had if traj[a].solve() below threw.
	const size_t ndim = PP::dim();
	if( ndim > 3 )
	    // Unreachable for every particle type that exists; PP::dim()
	    // is an inline function returning a literal, so this folds away
	    // entirely at -O2 rather than costing a branch per step.
	    throw( Error( ERROR_LOCATION, "particle dimension > 3 not supported here" ) );
	TrajectoryRep1D traj[3];
	for( size_t a = 0; a < ndim; a++ ) {
	    traj[a].construct( x2[0]-x1[0], 
			       x1[2*a+1], x1[2*a+2], 
			       x2[2*a+1], x2[2*a+2] );
	    DEBUG_MESSAGE( "Trajectory polynomial " << a << " order = " 
			   << traj[a].get_representation_order() << "\n" );
	}

	// Solve trajectory intersections
	for( size_t a = 0; a < ndim; a++ ) {

	    // Mesh number of x1 (start point)
	    int i = (int)floor( (x1[2*a+1]-mesh.origo(a))/mesh.h() );
	    
	    // Search to negative (dj = -1) and positive (dj = +1) mesh directions
	    for( int dj = -1; dj <= 1; dj += 2 ) {
		int j = i;
		if( dj == +1 )
		    j = i+1;
		int Kcount;  // Solution counter
		double K[3]; // Solution array
		while( 1 ) {

		    // Intersection point
		    double val = mesh.origo(a) + mesh.h() * j;
		    if( val < mesh.origo(a)-mesh.h() )
			break;
		    else if( val > mesh.max(a)+mesh.h() )
			break;

		    DEBUG_MESSAGE( "Searching intersections at coord(" << a << ") = " << val << "\n" );

		    Kcount = traj[a].solve( K, val );
		    if( Kcount == 0 )
			break; // No valid roots

		    // Save roots to coldata
		    DEBUG_MESSAGE( "Found " << Kcount << " valid roots\n" );
		    DEBUG_INC_INDENT();
		    for( int b = 0; b < Kcount; b++ ) {
			PP xcol;
			double x, v;
			xcol(0) = x1[0] + K[b]*(x2[0]-x1[0]);
			for( size_t c = 0; c < ndim; c++ ) {
			    traj[c].coord( x, v, K[b] );
			    if( a == c )
				xcol[2*c+1] = val; // limit numerical inaccuracy
			    else
				xcol[2*c+1] = x;
			    xcol[2*c+2] = v;
			}
			if( mesh.geom_mode() == MODE_CYL )
			    xcol[5] = x1[5] + K[b]*(x2[5]-x1[5]);

			DEBUG_MESSAGE( "K = " << K[b] << "\n" );
			DEBUG_MESSAGE( "Adding point " << xcol << "\n" );

			if( xcol[2*a+2] >= 0.0 )
			    coldata.push_back( ColData( xcol, a+1 ) );
			else
			    coldata.push_back( ColData( xcol, -a-1 ) );
		    }
		    DEBUG_DEC_INDENT();

		    j += dj;
		}
	    }
	}

	// Sort intersections in increasing time order
	sort( coldata.begin(), coldata.end() );

	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata built\n" );
    }

};


/*! \brief %Particle iterator class for continuous Vlasov-type iteration.
 *
 * Templated for particle point type (see ParticlePBase).
 *
 *  \todo Detailed documentation needed.
 *  \todo PIC style iterator needed.
 */
template <class PP> class ParticleIterator {

    typedef std::chrono::steady_clock Clock; /*!< \brief Lightweight clock for sub-timing accumulators. */

    gsl_odeiv2_system           _system;        /**< \brief GSL ODE integrator system. */
    gsl_odeiv2_step            *_step;          /**< \brief GSL ODE integrator stepper. */
    gsl_odeiv2_control         *_control;       /**< \brief GSL ODE integrator constrol. */
    gsl_odeiv2_evolve          *_evolve;        /**< \brief GSL ODE integrator integrator. */

    particle_iterator_type_e   _type;          /**< \brief Iteratory type. */

    trajectory_interpolation_e _intrp;         /*!< \brief Interpolation type. */
    scharge_deposition_e       _scharge_dep;   /*!< \brief Space charge deposition type. */
    double                     _epsabs;        /*!< \brief Absolute error limit. */
    double                     _epsrel;        /*!< \brief Relative error limit. */
    uint32_t                   _maxsteps;      /*!< \brief Maximum number of simulation steps for particle. */
    double                     _maxt;          /*!< \brief Maximum particle lifetime. */
    uint32_t                   _degenerate_bracket = 0; /*!< \brief Zero-length collision brackets seen,
							 *   limits the warning. */
    uint32_t                   _stuck_reported = 0;     /*!< \brief Step size underflows seen,
							 *   limits the warning. */
    bool                       _mask_reflect = true;    /*!< \brief Reflect particles at Neumann masks
							 *   (default). If false, a mask absorbs
							 *   like any other solid. */

    bool                       _save_points;   /*!< \brief Save all points? */
    uint32_t                   _trajdiv;       /*!< \brief Divisor for saved trajectories,
					        * if 3, every third trajectory is saved. */
    bool                       _mirror[6];     /*!< \brief Is particle mirrored on boundary? */

    std::vector<uint32_t>      _mask_solid;     /*!< \brief Solid numbers of the Neumann masks. */
    std::vector<uint32_t>      _analytic_solid; /*!< \brief Solid numbers of the solids with an
						 *   analytic surface. */
    bool                       _mask_reflected; /*!< \brief Set by check_collision_surface() if it
						 *   reflected the particle from a Neumann mask. */
    bool                       _surface_collision;

    ParticleIteratorData       _pidata;        /*!< \brief User data provided to PP::get_derivatives(). */
    TrajectoryHandlerCallback *_thand_cb;      /*!< \brief Trajectory handler callback. */
    TrajectoryEndCallback     *_tend_cb;       /*!< \brief Trajectory end callback. */
    TrajectorySurfaceCollisionCallback *_tsur_cb;   /*!< \brief Trajectory surface collision callback. */
    const TrajectoryEndCallback *_bsup_cb;     /*!< \brief B-field plasma suppression callback. */
    ParticleDataBase          *_pdb;           /*!< \brief Particle database pointer for adding secondary particles. */

    PP                         _xi;            /*!< \brief Previous mesh intersection coordinates 
					        *   or starting point. */
    std::vector<PP>            _traj;          /*!< \brief %Particle trajectory data for current trajectory. */
    std::vector<ColData<PP> >  _coldata;       /*!< \brief %Mesh intersection coordinate data. */
    CFiFo<PP,4>                _cdpast;        /*!< \brief Past three intersection coords. */

    ParticleStatistics         _stat;          /*!< \brief Particle statistics. */

    /* Timing accumulators. An iterator serves one thread only, so they
     * need no locking. The caller sums them over the iterators. */
    double                     _time_ode;        /*!< \brief Time spent in gsl_odeiv2_evolve_apply() (ODE stepping). */
    double                     _time_trajhandle;  /*!< \brief Time spent in handle_trajectory() (collision detection + space charge deposition). */

    /*! \brief Save trajectory point \a x.
     *
     *  Throw error if run out of memory.
     */
    void save_trajectory_point( PP x ) {

	try {
	    _traj.push_back( x );
	} catch( const std::bad_alloc & ) {
	    throw( ErrorNoMem( ERROR_LOCATION, "Out of memory saving trajectory" ) );
	}
    }

    /*! \brief Outward surface normal of analytic solid \a n at \a p, by
     *  differencing the indicator. On a face this is exactly the axis
     *  normal, which is what a box-shaped mask needs. */
    Vec3D analytic_normal( uint32_t n, const Vec3D &p ) {
	const double eps = 1.0e-3*_pidata._geom->h();
	Vec3D nv( 0.0, 0.0, 0.0 );
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    Vec3D pp = p, pm = p;
	    pp[a] += eps;
	    pm[a] -= eps;
	    nv[a] = (double)(_pidata._geom->inside( n, pp ) ? 1 : 0)
		  - (double)(_pidata._geom->inside( n, pm ) ? 1 : 0);
	}
	double l = nv.norm2();
	if( l > 0.0 )
	    nv /= l;
	return( nv );
    }

    /*! \brief First triangulated solid with a node on cube (i,j,k), 0 if
     *  none. A triangle hit in the cube is attributed to it. */
    uint32_t get_solid_mc( int i, int j, int k ) {
	for( int cz = 0; cz < 2; cz++ )
	    for( int cy = 0; cy < 2; cy++ )
		for( int cx = 0; cx < 2; cx++ ) {
		    uint32_t m = _pidata._geom->material( i+cx, j+cy, k+cz );
		    if( m != 0 && !_pidata._geom->solid_is_analytic( m ) )
			return( m );
		}
	return( 0 );
    }

    /*! \brief Check for particle collision with triangulated solid surface.
     *
     *  Checks for intersection between trajectory segment and surface
     *  triangles in grid cell i.
     */
    bool check_collision_surface( Particle<PP> &particle, const PP &x1, const PP &x2, 
				  PP &status_x, const int32_t i[3] ) {

	DEBUG_MESSAGE( "Checking collisions with solids - surface check\n" );
	DEBUG_INC_INDENT();

	if( i[0] < 0 )
	    return( true );
	else if( i[0] >= (int32_t)_pidata._geom->size(0)-1 )
	    return( true );
	else if( i[1] < 0 )
	    return( true );
	else if( i[1] >= (int32_t)_pidata._geom->size(1)-1 )
	    return( true );
	else if( i[2] < 0 )
	    return( true );
	else if( i[2] >= (int32_t)_pidata._geom->size(2)-1 )
	    return( true );

	Vec3D v1 = x1.location();
	Vec3D v2 = x2.location();

	DEBUG_MESSAGE( "v1 = " << v1 << "\n" );
	DEBUG_MESSAGE( "v2 = " << v2 << "\n" );

	/* Every solid that can be met in this cube is tested and the
	 * nearest hit wins: a cube can touch a mask and an electrode at
	 * once, and the surface reached first decides between reflection
	 * and absorption. Both tests below give the hit as the fraction t
	 * of the way from x1 to x2. */
	double tbest = 2.0;
	uint32_t sbest = 0;
	Vec3D nbest;
	bool nbest_valid = false;
	int32_t tribest = -1;
	double k1best = 0.0, k2best = 0.0;
	PP xbest;

	/* Analytic solids, all of them. They cannot be selected by the
	 * node materials of the cube: the surface of an analytic solid is
	 * on the plane of its last free node (set_link() mirrors about
	 * that node), so the face through which a particle enters it has
	 * no solid node. */
	for( size_t q = 0; q < _analytic_solid.size(); q++ ) {
	    const uint32_t m = _analytic_solid[q];

	    /* Only an entry is a hit. bracket_surface() needs x1 to be
	     * outside the solid. */
	    if( _pidata._geom->inside( m, v1 ) )
		continue;

	    /* The end of the segment is probed slightly past x2: x2 is on
	     * the node plane where the surface is, and the solid may not
	     * count that plane as inside. */
	    Vec3D dseg = v2-v1;
	    double dl = dseg.norm2();
	    Vec3D probe = v2;
	    if( dl > 0.0 )
		probe = v2 + (1.0e-3*_pidata._geom->h()/dl)*dseg;
	    if( !_pidata._geom->inside( m, probe ) )
		continue;

	    Vec3D vc;
	    if( v1 == probe )
		vc = probe;
	    else
		_pidata._geom->bracket_surface( m, probe, v1, vc );
	    double seg2 = dseg*dseg;
	    double t = ( seg2 > 0.0 ) ? ((vc-v1)*dseg)/seg2 : 0.0;
	    if( t < 0.0 ) t = 0.0;
	    if( t > 1.0 ) t = 1.0;
	    if( t < tbest ) {
		tbest = t;
		sbest = m;
		nbest_valid = false;
		tribest = -1;
		for( size_t b = 0; b < PP::size(); b++ )
		    xbest[b] = x1[b] + t*(x2[b]-x1[b]);
	    }
	}

	/* Triangulated solids: the surface triangles of this cube. */
	int32_t tric = _pidata._geom->surface_trianglec( i[0], i[1], i[2] );
	int32_t ptr = _pidata._geom->surface_triangle_ptr( i[0], i[1], i[2] );

	for( int32_t a = 0; a < tric; a++ ) {
	    const VTriangle &tri = _pidata._geom->surface_triangle( ptr+a );
	    const Vec3D &va = _pidata._geom->surface_vertex( tri[0] );
	    const Vec3D &vb = _pidata._geom->surface_vertex( tri[1] );
	    const Vec3D &vc = _pidata._geom->surface_vertex( tri[2] );

	    Mat3D m( v2[0]-v1[0], -vb[0]+va[0], -vc[0]+va[0],
		     v2[1]-v1[1], -vb[1]+va[1], -vc[1]+va[1],
		     v2[2]-v1[2], -vb[2]+va[2], -vc[2]+va[2] );
	    double mdet = m.determinant();
	    if( mdet == 0.0 )
		continue;
	    Mat3D minv = m.inverse( mdet );
	    Vec3D off( -v1[0]+va[0], -v1[1]+va[1], -v1[2]+va[2] );
	    Vec3D K = minv*off;
	    double K3 = K[1]+K[2];

	    if( K[0] > -COLLISION_EPS && K[0] < 1.0+COLLISION_EPS &&
		K[1] > -COLLISION_EPS && K[1] < 1.0+COLLISION_EPS &&
		K[2] > -COLLISION_EPS && K[2] < 1.0+COLLISION_EPS &&
		K3 > -COLLISION_EPS && K3 < 1.0+COLLISION_EPS &&
		K[0] < tbest ) {
		tbest = K[0];
		sbest = get_solid_mc( i[0], i[1], i[2] );
		nbest = cross( vb-va, vc-va );
		nbest_valid = true;
		tribest = ptr+a;
		k1best = K[1];
		k2best = K[2];
		for( size_t b = 0; b < PP::size(); b++ )
		    xbest[b] = x1[b] + K[0]*(x2[b]-x1[b]);
	    }
	}

	if( sbest == 0 ) {
	    DEBUG_MESSAGE( "No collisions\n" );
	    DEBUG_DEC_INDENT();
	    return( true );
	}

	status_x = xbest;

	for( int32_t b = _traj.size()-1; b > 0; b-- ) {
	    if( _traj[b][0] > status_x[0] )
		_traj.pop_back();
	    else
		break;
	}

	/* A Neumann mask is a symmetry surface and reflects the particle.
	 * Any other solid absorbs it. */
	if( _mask_reflect && is_mask_solid( sbest ) ) {
	    Vec3D n = nbest_valid ? nbest
				  : analytic_normal( sbest, status_x.location() );
	    double nn = n.norm2();
	    if( nn > 0.0 ) {
		n /= nn;
		double vn = 0.0;
		for( size_t b = 0; b < PP::dim(); b++ )
		    vn += status_x[2*b+2]*n[b];
		const double sgn = ( vn > 0.0 ? -1.0 : 1.0 );
		for( size_t b = 0; b < PP::dim(); b++ ) {
		    status_x[2*b+2] -= 2.0*vn*n[b];
		    status_x[2*b+1] += sgn*1.0e-3*_pidata._geom->h()*n[b];
		}
	    }
	    _mask_reflected = true;
	    DEBUG_MESSAGE( "Neumann surface reflection\n" );
	    DEBUG_DEC_INDENT();
	    return( false );
	}

	particle.set_status( PARTICLE_COLL );
	_stat.add_bound_collision( sbest, particle.IQ() );

	if( _tsur_cb && tribest >= 0 )
	    (*_tsur_cb)( &particle, &status_x, tribest, k1best, k2best );

	DEBUG_MESSAGE( "Solid collision detected\n" );
	DEBUG_DEC_INDENT();
	return( false );
    }

    /*! \brief Check for particle collision with solid surface with inside().
     *
     *  Particle propagates from x1 to x2, where x1 is in
     *  vacuum. Checks if x2 is inside solid and if it is, the
     *  collision point is bracketed between x1 and x2. Particle
     *  coordinates at status_x are set to collision coordinates.
     *  Returns false if particle collided.
     */
    bool check_collision_solid( Particle<PP> &particle, const PP &x1, const PP &x2, 
				PP &status_x ) {

	DEBUG_MESSAGE( "Checking collisions with solids - inside check\n" );
	DEBUG_INC_INDENT();

	// If inside solid, bracket for collision point
	Vec3D v2 = x2.location();
	int32_t bound = _pidata._geom->inside( v2 );
	if( bound < 7 ) {
	    DEBUG_MESSAGE( "No collisions\n" );
	    DEBUG_DEC_INDENT();
	    return( true );
	}
	Vec3D vc;
	Vec3D v1 = x1.location();

	// x1 and x2 can be the same point, which bracket_surface() cannot
	// bracket. The point is inside the solid, so the collision is
	// taken to be there (K = 0). It also means that x1 was not in
	// vacuum as it should be, hence the warning.
	if( v1 == v2 && _degenerate_bracket < 10 ) {
	    _degenerate_bracket++;
	    ibsimu.message( MSG_WARNING, 1 )
		<< "Warning: zero-length collision bracket, solid " << bound
		<< " at " << v2 << "; treating as a collision there.\n"
		<< "         x1 = " << x1 << "\n"
		<< "         x2 = " << x2 << "\n"
		<< ( _degenerate_bracket == 10
		     ? "         Further occurrences on this thread suppressed.\n"
		     : "" );
	} else if( v1 == v2 ) {
	    _degenerate_bracket++;
	}

	double K = ( v1 == v2 )
	    ? 0.0
	    : _pidata._geom->bracket_surface( bound, v2, v1, vc );

	// Calculate new PP
	for( size_t a = 0; a < PP::size(); a++ )
	    status_x[a] = x2[a] + K*(x1[a]-x2[a]);

	// Remove all points from trajectory after time status_x[0].
	for( size_t a = _traj.size()-1; a > 0; a-- ) {
	    if( _traj[a][0] > status_x[0] )
		_traj.pop_back();
	    else
		break;
	}

	// Save last trajectory point and update status
	//save_trajectory_point( status_x );
	particle.set_status( PARTICLE_COLL );

	// Update collision statistics for boundary
	_stat.add_bound_collision( bound, particle.IQ() );

	DEBUG_MESSAGE( "Solid collision detected\n" );
	DEBUG_DEC_INDENT();

	return( false );
    }

    /*! \brief Check for particle collision with solid
     *
     *  Two different algorithms.
     */
    bool check_collision( Particle<PP> &particle, const PP &x1, const PP &x2, 
			  PP &status_x, const int i[3] ) {

	if( _surface_collision )
	    return( check_collision_surface( particle, x1, x2, status_x, i ) );
	return( check_collision_solid( particle, x1, x2, status_x ) );
    }


    /*! \brief End the current step at mesh crossing \a c.
     *
     *  For use when the particle coordinates are changed discontinuously
     *  at the crossing. The rest of the crossings describe a path that
     *  is not taken, and the history of the adaptive stepper is not
     *  valid any more. The crossing point has its own time, so the next
     *  step continues from the instant of the crossing.
     */
    void truncate_step( size_t c ) {
	save_trajectory_point( _coldata[c]._x );
	_coldata.resize( c+1 );      // ends the loop of handle_trajectory()
	gsl_odeiv2_step_reset( _step );
	gsl_odeiv2_evolve_reset( _evolve );
    }


    /*! \brief Reflect the particle at mesh crossing \a c from the plane
     *  x[a] = \a plane, which it met moving in direction \a side (-1 or
     *  +1) along axis \a a.
     *
     *  Serves the mirror boundaries of the simulation box and the
     *  Neumann masks alike. The particle is stopped at the plane,
     *  slightly on the side it came from, with its normal velocity
     *  reversed, and the step ends there. The part of the step beyond
     *  the plane is not mirrored back: beyond a mask there is no field
     *  solution to have followed, and in cylindrical symmetry the image
     *  can land on r <= 0, where get_derivatives() fails.
     */
    void reflect_at_plane( size_t c, size_t a, int side, double plane, PP &x2 ) {

	DEBUG_MESSAGE( "Reflect trajectory at plane " << plane << "\n" );

	const double eps = 1.0e-3 * _pidata._geom->h();
	x2 = _coldata[c]._x;
	x2[2*a+1] = plane + ( side < 0 ? eps : -eps );
	x2[2*a+2] *= -1.0;
	truncate_step( c );
    }


    void handle_collision( Particle<PP> &particle, uint32_t bound, size_t c, PP &status_x ) {

	DEBUG_MESSAGE( "Handle collision\n" );

	//save_trajectory_point( _coldata[c]._x );
	status_x = _coldata[c]._x;
	particle.set_status( PARTICLE_OUT );
	_stat.add_bound_collision( bound, particle.IQ() );
    }


    /*! \brief Is any node of a face of mesh cell \a i a solid node?
     *
     *  The face is the one across axis \a a, on the low side of the
     *  cell for \a side = -1 and on the high side for +1.
     */
    bool face_has_solid( const int i[3], size_t a, int side ) {

	const Geometry &geom = *_pidata._geom;
	int n[3] = { i[0], i[1], i[2] };
	if( side > 0 )
	    n[a]++;

	if( PP::dim() == 2 ) {
	    if( geom.material_check( n[0], n[1] ) != 0 )
		return( true );
	    n[1-a]++;
	    return( geom.material_check( n[0], n[1] ) != 0 );
	}

	const size_t b = (a+1)%3;
	const size_t d = (a+2)%3;
	for( int kd = 0; kd < 2; kd++ ) {
	    for( int kb = 0; kb < 2; kb++ ) {
		int m[3] = { n[0], n[1], n[2] };
		m[b] += kb;
		m[d] += kd;
		if( geom.material_check( m[0], m[1], m[2] ) != 0 )
		    return( true );
	    }
	}
	return( false );
    }


    /*! \brief Is the particle entering a Neumann mask at mesh crossing
     *  \a c, which it makes in direction \a side (-1 or +1) along axis
     *  \a a?
     *
     *  Asked from the mask solids at a point slightly past the crossing,
     *  not from the node materials: the two disagree between the surface
     *  of a mask and its first solid node, and check_collision_solid()
     *  absorbs by the solid. Only the masks are asked, so the cost does
     *  not depend on the other solids.
     */
    bool entering_mask( size_t c, size_t a, int side ) {
	Vec3D v = _coldata[c]._x.location();
	v[a] += side * 1.0e-3 * _pidata._geom->h();
	for( size_t s = 0; s < _mask_solid.size(); s++ ) {
	    if( _pidata._geom->inside( _mask_solid[s], v ) )
		return( true );
	}
	return( false );
    }


    /*! \brief Is solid number \a n a Neumann mask?
     */
    bool is_mask_solid( uint32_t n ) const {
	return( std::find( _mask_solid.begin(), _mask_solid.end(), n ) != _mask_solid.end() );
    }


    /*! \brief Handle particle mesh intersection.
     *
     *  Particle mesh coordinates \a i are advanced through
     *  intersection described in \a _coldata[c]. Collisions with
     *  solids and particle mirroring on boundaries are handled. Final
     *  particle coordinates \a x2 are updated accordingly. Returns
     *  false if particle collided.
     */
    bool handle_trajectory_advance( Particle<PP> &particle, size_t c, int i[3], PP &x2 ) {

	DEBUG_MESSAGE( "Handle trajectory advance\n" );
	DEBUG_INC_INDENT();

	// The particle leaves cell i through its face across axis a, on
	// the low (side = -1) or high (side = +1) side of the cell.
	const int dir = _coldata[c]._dir;
	const size_t a = (size_t)((dir < 0 ? -dir : dir) - 1);
	const int side = ( dir < 0 ? -1 : +1 );

	// Solid collision model: a mask reflects here, before the
	// collision test, which would absorb on it. The surface collision
	// model finds the masks in check_collision_surface() instead.
	if( !_mask_solid.empty() && _mask_reflect && !_surface_collision &&
	    entering_mask( c, a, side ) ) {
	    reflect_at_plane( c, a, side, _coldata[c]._x[2*a+1], x2 );
	    DEBUG_DEC_INDENT();
	    return( true );
	}

	// Test for a collision if the face has a solid node. With the
	// surface collision model also if there is any analytic solid,
	// see check_collision_surface().
	_mask_reflected = false;
	const bool collided =
	    ( face_has_solid( i, a, side ) || (_surface_collision && !_analytic_solid.empty()) ) &&
	    !check_collision( particle, _xi, _coldata[c]._x, x2, i );

	if( _mask_reflected ) {
	    // check_collision_surface() reflected the particle from a
	    // mask and set x2. The step ends there, before the crossing.
	    _coldata[c]._x = x2;
	    truncate_step( c );
	    DEBUG_DEC_INDENT();
	    return( true );
	}

	i[a] += side;

	if( collided ) {
	    DEBUG_MESSAGE( "Surface collision!\n" );
	    DEBUG_DEC_INDENT();
	    return( false );
	}

	// Simulation box boundary: i is now outside the cell range
	// 0...size-2 on the axis where the particle left the box.
	for( size_t b = 0; b < PP::dim(); b++ ) {

	    int border;
	    if( i[b] < 0 )
		border = -1;
	    else if( i[b] >= (int32_t)(_pidata._geom->size(b)-1) )
		border = +1;
	    else
		continue;

	    const uint32_t bound = ( border < 0 ? 1+2*b : 2+2*b );
	    DEBUG_MESSAGE( "Boundary collision at boundary " << bound << "\n" );
	    DEBUG_DEC_INDENT();
	    if( !_mirror[bound-1] ) {
		handle_collision( particle, bound, c, x2 );
		return( false );
	    }

	    // A reflected particle stays in the first or last cell.
	    i[b] = ( border < 0 ? 0 : (int32_t)_pidata._geom->size(b)-2 );
	    reflect_at_plane( c, b, border,
			      border < 0 ? _pidata._geom->origo(b) : _pidata._geom->max(b), x2 );
	    return( true );
	}

	DEBUG_MESSAGE( "No collision.\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }

    /*! \brief Limit trajectory advance from \a x1 to \a x2 to double
     *  the simulation box.
     *
     *  Return true if limitation is done, false if data is left untouched.
     */
    bool limit_trajectory_advance( const PP &x1, PP &x2 ) {

	bool touched = false;
	for( size_t a = 0; a < PP::dim(); a++ ) {

	    double lim1 = _pidata._geom->origo(a) - 
		(_pidata._geom->size(a)-1)*_pidata._geom->h();
	    double lim2 = _pidata._geom->origo(a) + 
		2*(_pidata._geom->size(a)-1)*_pidata._geom->h();

	    if( x2[2*a+1] < lim1 ) {
		double K = (lim1 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
		DEBUG_MESSAGE( "Limiting step to " << x2 << "\n" );
	    } else if(x2[2*a+1] > lim2 ) {
		double K = (lim2 - x1[2*a+1]) / (x2[2*a+1] - x1[2*a+1]);
		x2 = x1 + K*(x2-x1);
		touched = true;
		DEBUG_MESSAGE( "Limiting step to " << x2 << "\n" );
	    }
	}

	return( touched );
    }

    /*! \brief Build coldata.
     *
     *  Throw error if run out of memory.
     */
    void build_coldata( bool force_linear, const PP &x1, const PP &x2 ) {

	try {
	    if( _intrp == TRAJECTORY_INTERPOLATION_POLYNOMIAL && !force_linear )
		ColData<PP>::build_coldata_poly( _coldata, *_pidata._geom, x1, x2 );
	    else
		ColData<PP>::build_coldata_linear( _coldata, *_pidata._geom, x1, x2 );
	} catch( const std::bad_alloc & ) {
	    throw( ErrorNoMem( ERROR_LOCATION, "out of memory building collision data" ) );
	}
    }


    /*! \brief Handle particle iteration step from coordinates \a x1 to \a x2.
     *
     *  Searches mesh intersections between points \a x1 and \a x2 and
     *  builds ColData. Checks for collisions with solids and
     *  boundaries and sets space charge. If particle collides with
     *  mirroring boundary, the \a x2 coordinates are changed and GSL
     *  ODE integrator is resetted. Space charge is deposited at each
     *  mesh line crossing.
     *
     *  If \a force_linear is true, linear interpolation of trajectory
     *  is used regardless of interpolation settings.
     *
     *  If \a first_step is true, this is the first step and the
     *  particle is allowed to get into the simulation box if the
     *  definition point (x1) was outside the geometry.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_trajectory( Particle<PP> &particle, const PP &x1, PP &x2, 
			    bool force_linear, bool first_step ) {

	DEBUG_MESSAGE( "Handle trajectory from x1 to x2:\n" <<
		       "  x1 = " << x1 << "\n" << 
		       "  x2 = " << x2 << "\n" );
	DEBUG_INC_INDENT();

	// Limit trajectory advance to double the simulation box
	// If limitation done, force to linear interpolation
	if( limit_trajectory_advance( x1, x2 ) )
	    force_linear = true;

	build_coldata( force_linear, x1, x2 );

	// TODO
	// Remove entrance to geometry if coming from outside or make
	// code to skip the collision detection for these particles

	// No intersections, nothing to do
	if( _coldata.size() == 0 ) {
	    DEBUG_MESSAGE( "No coldata\n" );
	    DEBUG_DEC_INDENT();
	    return( true );
	}

	// Starting mesh index
	int i[3] = {0, 0, 0};
	for( size_t cdir = 0; cdir < PP::dim(); cdir++ )
	    i[cdir] = (int)floor( (x1[2*cdir+1]-_pidata._geom->origo(cdir))*_pidata._geom->div_h() );

	// Process intersection points
	DEBUG_MESSAGE( "Process coldata points:\n" );
	DEBUG_INC_INDENT();
	for( size_t a = 0; a < _coldata.size(); a++ ) {

	    if( _save_points )
		save_trajectory_point( _coldata[a]._x );

	    DEBUG_MESSAGE( "Coldata " << a << "\n" <<
			   "  x = " << _coldata[a]._x << "\n" <<
			   "  i = (" << std::setw(3) << i[0] << ", "
			   << std::setw(3) << i[1] << ", "
			   << std::setw(3) << i[2] << ")\n" << 
			   "  dir = " << _coldata[a]._dir << "\n" );

	    // Update space charge for mesh volume i
	    if( _scharge_dep == SCHARGE_DEPOSITION_LINEAR && _pidata._scharge ) {
		_cdpast.push( _coldata[a]._x );
		scharge_add_from_trajectory_linear( *_pidata._scharge, particle.IQ(),
						    _coldata[a]._dir, _cdpast, i );
	    }

	    // Advance particle in mesh, update i, check for possible
	    // collisions and do mirroring for trajectory and coldata.
	    handle_trajectory_advance( particle, a, i, x2 );

	    if( _scharge_dep == SCHARGE_DEPOSITION_PIC && _pidata._scharge ) {
		scharge_add_from_trajectory_pic( *_pidata._scharge, particle.IQ(),
						 _xi, _coldata[a]._x );
	    }

	    // Call trajectory handler callback
	    if( _thand_cb )
		(*_thand_cb)( &particle, &_coldata[a]._x, &x2 );

	    // Clear coldata and exit if particle collided.
	    if( particle.get_status() != PARTICLE_OK ) {
		save_trajectory_point( x2 );
		_coldata.clear();
		DEBUG_DEC_INDENT();
		DEBUG_MESSAGE( "Interrupting coldata processing\n" );
		DEBUG_DEC_INDENT();
		return( false );
	    }

	    // Update last accepted intersection point xi.
	    _xi = _coldata[a]._x;
	}

	_coldata.clear();
	DEBUG_DEC_INDENT();
	DEBUG_MESSAGE( "Coldata done\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }


    /*! \brief Is particle mirroring required at axis in cylindrical symmetry?
     */
     bool axis_mirror_required( const PP &x2 ) {
	 return( _pidata._geom->geom_mode() == MODE_CYL && 
		 x2[4] < 0.0 && 
		 x2[3] <= 0.01*_pidata._geom->h() &&
		 x2[3]*fabs(x2[5]) <= 1.0e-9*fabs(x2[4]) );
		 
     }


    /*! \brief Handle particle mirrored at axis in cylindrical symmetry.
     *
     *  Return true if particle status is PARTICLE_OK after trajectory
     *  step, false otherwise.
     */
    bool handle_axis_mirror_step( Particle<PP> &particle, const PP &x1, PP &x2 ) {

	// Get acceleration at x2
	double dxdt[5];
	PP::get_derivatives( x2[0], &x2[1], dxdt, (void *)&_pidata );

	// Calculate crossover point assuming zero acceleration in
	// r-direction and constant acceleration in x-direction
	double dt = -x2[3]/x2[4];
	PP xc;
	xc.clear();
	xc[0] = x2[0]+dt;
	xc[1] = x2[1]+(x2[2]+0.5*dxdt[1]*dt)*dt;
	xc[2] = x2[2];
	xc[3] = x2[3]+x2[4]*dt;
	xc[4] = x2[4];
	xc[5] = x2[5];

	// Mirror x2 to x3
	PP x3 = 2*xc - x2;
	x3[3] *= -1.0;
	x3[4] *= -1.0;
	x3[5] *= -1.0;

	DEBUG_MESSAGE( "Particle mirror:\n" <<
		       "  x1: " << x1 << "\n" <<
		       "  x2: " << x2 << "\n" << 
		       "  xc: " << xc << "\n" << 
		       "  x3: " << x3 << "\n" );
	
	// Handle step with linear interpolation to avoid going to r<=0
	{
	    Clock::time_point time_t0 = Clock::now();
	    bool ok = handle_trajectory( particle, x2, x3, true, false );
	    _time_trajhandle += std::chrono::duration<double>( Clock::now()-time_t0 ).count();
	    if( !ok )
		return( false ); // Particle done
	}

	// Save trajectory calculation points
	save_trajectory_point( x2 );
	save_trajectory_point( xc );
	xc[4] *= -1.0;
	xc[5] *= -1.0;
	save_trajectory_point( xc );
	
	// Next step not a continuation of previous one, reset
	// integrator
	gsl_odeiv2_step_reset( _step );
	gsl_odeiv2_evolve_reset( _evolve );

	// Continue iteration at mirrored point
	x2 = x3;
	return( true );
    }
    
    /*! \brief Check particle definition.
     *
     *  False is returned if particle is defined out of simulation
     *  area or inside a solid. False is also returned if illegal
     *  input (NaN) is received.
     *
     *  Removed feature: False is also returned if a particle
     *  is defined on the simulation border going out of simulation
     *  with the corresponding border not having mirroring enabled. If
     *  the mirroring is enabled, the particle is mirrored and true is
     *  returned. Otherwise the particle definition is ok and true is
     *  returned.
     */
    bool check_particle_definition( PP &x ) {

	DEBUG_MESSAGE( "Particle defined at x = " << x << "\n" );
	DEBUG_INC_INDENT();

	// Check for NaN
	for( size_t a = 0; a < PP::size(); a++ ) {
	    if( comp_isnan( x[a] ) ) {
		DEBUG_MESSAGE( "Particle coordinate NaN\n" );
		DEBUG_DEC_INDENT();
		return( false );
	    }
	}

	// Check if inside solids or outside simulation geometry box.
	if( _surface_collision ) {
	    if( _pidata._geom->surface_inside( x.location() ) ) {
		DEBUG_MESSAGE( "Particle inside solid or outside simulation\n" );
		DEBUG_DEC_INDENT();
		return( false );
	    }
	} else { 
	    if( _pidata._geom->inside( x.location() ) ) {
		DEBUG_MESSAGE( "Particle inside solid or outside simulation\n" );
		DEBUG_DEC_INDENT();
		return( false );
	    }
	}

	// Check if particle on simulation geometry border and directed outwards
	/*
	for( size_t a = 0; a < PP::dim(); a++ ) {
	    if( x[2*a+1] == _pidata._geom->origo(a) && x[2*a+2] < 0.0 ) {
		if( _mirror[2*a] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Mirroring to:\n";
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}

	    } else if( x[2*a+1] == _pidata._geom->max(a) & x[2*a+2] > 0.0 ) {
		if( _mirror[2*a+1] ) {
		    x[2*a+2] *= -1.0;
#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Mirroring to:\n";
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "  x = " << x << "\n";
#endif
		} else {
		    return( false );
		}
	    }
	}

#ifdef DEBUG_PARTICLE_ITERATOR
	ibsimu.message(MSG_DEBUG_GENERAL,1) << "Definition ok\n";
#endif

	*/
	DEBUG_MESSAGE( "Ok\n" );
	DEBUG_DEC_INDENT();
	return( true );
    }
    
    double calculate_dt( const PP &x, const double *dxdt ) {

	double spd = 0.0, acc = 0.0;

	for( size_t a = 0; a < (PP::size()-1)/2; a++ ) {
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd += " << dxdt[2*a]*dxdt[2*a] << "\n";
	    spd += dxdt[2*a]*dxdt[2*a];
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc += " << dxdt[2*a+1]*dxdt[2*a+1] << "\n";
	    acc += dxdt[2*a+1]*dxdt[2*a+1];
	}
	if( _pidata._geom->geom_mode() == MODE_CYL ) {
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "MODE_CYL\n";
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd += " << x[3]*x[3]*x[5]*x[5] << "\n";
	    spd += x[3]*x[3]*x[5]*x[5];
	    //ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc += " << x[3]*x[3]*dxdt[4]*dxdt[4] << "\n";
	    acc += x[3]*x[3]*dxdt[4]*dxdt[4];
	}
	//ibsimu.message(MSG_DEBUG_GENERAL,1) << "spd = " << sqrt(spd) << "\n";
	//ibsimu.message(MSG_DEBUG_GENERAL,1) << "acc = " << sqrt(acc) << "\n";
	spd = _pidata._geom->h() / sqrt(spd);
	acc = sqrt( 2.0*_pidata._geom->h() / sqrt(acc) );

	return( spd < acc ? spd : acc );
    }

public:

    /*! \brief Constructor for new particle iterator.
     *
     *  New particle iterator is initialized with given settings.
     *
     *  \param type Particle iterator type used
     *  \param epsabs Absolute error limit in iteration
     *  \param epsrel Relative error limit in iteration
     *  \param intrp Interpolation type.
     *  \param scharge_dep Space charge deposition type.
     *  \param maxsteps Maximum number of steps to take before particle is killed
     *  \param maxt Maximum flight time for a particle
     *  \param save_points Flag for saving all intersection points of trajectories
     *  \param trajdiv Trajectory saving divisor. Only every trajdiv:th particle 
     *  trajectory saved.
     *  \param mirror %Particle mirroring on surfaces
     *  \param scharge Space charge field to save to
     *  \param efield Electric field in the geometry
     *  \param bfield Magnetic field in the geometry
     *  \param geom %Geometry definition
     *
     *  The particle iterator is given the settings for calculation
     *  and geometry, electric field and space charge map to
     *  build. Pointer to first particle in the particle database
     *  vector is used to calculate the particle number from the
     *  particle memory location.
     */
    ParticleIterator( particle_iterator_type_e type, double epsabs, double epsrel,
		      trajectory_interpolation_e intrp, scharge_deposition_e scharge_dep,
		      uint32_t maxsteps, double maxt, bool save_points,
		      uint32_t trajdiv, bool mirror[6], MeshScalarField *scharge,
		      const VectorField *efield, const VectorField *bfield,
		      const Geometry *geom )
	: _type(type), _intrp(intrp), _scharge_dep(scharge_dep), _epsabs(epsabs), _epsrel(epsrel),
	  _maxsteps(maxsteps), _maxt(maxt), _save_points(save_points), _trajdiv(trajdiv),
	  _mask_reflected(false), _surface_collision(false), _pidata(scharge,efield,bfield,geom),
	  _thand_cb(0), _tend_cb(0), _tsur_cb(0), _bsup_cb(0), _pdb(0),
	  _stat(geom->number_of_boundaries()),
	  _time_ode(0.0), _time_trajhandle(0.0) {

	for( uint32_t n = 7; n <= geom->number_of_boundaries(); n++ ) {
	    if( geom->get_boundary(n).type() == BOUND_NEUMANN )
		_mask_solid.push_back( n );
	    if( geom->solid_is_analytic( n ) )
		_analytic_solid.push_back( n );
	}

	// Initialize mirroring
	_mirror[0] = mirror[0];
	_mirror[1] = mirror[1];
	_mirror[2] = mirror[2];
	_mirror[3] = mirror[3];
	_mirror[4] = mirror[4];
	_mirror[5] = mirror[5];

	// Initialize system of ordinary differential equations (ODE)
	_system.jacobian  = NULL;
	_system.params    = (void *)&_pidata;
	_system.function  = PP::get_derivatives;
	_system.dimension = PP::size()-1; // Time is not part of differential equation dimensions

	// Make scale
	// 2D:  x vx y vy
	// Cyl: x vx r vr omega
	// 3D:  x vx y vy z vz
	double scale_abs[PP::size()-1];
	for( uint32_t a = 0; a < (uint32_t)PP::size()-2; a+=2 ) {
	    scale_abs[a+0] = 1.0;
	    scale_abs[a+1] = 1.0e6;
	}
	if( _pidata._geom->geom_mode() == MODE_CYL )
	    scale_abs[4] = 1.0;

	// Initialize ODE solver
	_step    = gsl_odeiv2_step_alloc( gsl_odeiv2_step_rkck, _system.dimension );
	//_control = gsl_odeiv2_control_standard_new( _epsabs, _epsrel, 1.0, 1.0 );
	_control = gsl_odeiv2_control_scaled_new( _epsabs, _epsrel, 1.0, 1.0, scale_abs, PP::size()-1 );
	_evolve  = gsl_odeiv2_evolve_alloc( _system.dimension );
    }


    /*! \brief Destructor.
     */
    ~ParticleIterator() {
	gsl_odeiv2_evolve_free( _evolve );
	gsl_odeiv2_control_free( _control );
	gsl_odeiv2_step_free( _step );
    }


    /*! \brief Enable/disable surface collision model.
     */
    void set_surface_collision( bool surface_collision ) {
	if( surface_collision && _pidata._geom->geom_mode() == MODE_2D )
	    throw( Error( ERROR_LOCATION, "2D surface collision not supported" ) );
	if( surface_collision && !_pidata._geom->surface_built() )
	    throw( Error( ERROR_LOCATION, "surface model not built" ) );
	_surface_collision = surface_collision;
    }


    /*! \brief Reflect particles at Neumann masks (true, default) or
     *  absorb them (false).
     */
    void set_mask_reflection( bool reflect ) { _mask_reflect = reflect; }

    /*! \brief Set trajectory handler callback. 
     */
    void set_trajectory_handler_callback( TrajectoryHandlerCallback *thand_cb ) {
	_thand_cb = thand_cb;
    }


    /*! \brief Set trajectory end callback. 
     */
    void set_trajectory_end_callback( TrajectoryEndCallback *tend_cb, ParticleDataBase *pdb ) {
	_tend_cb = tend_cb;
	_pdb = pdb;
    }


    /*! \brief Set trajectory surface collision callback. 
     */
    void set_trajectory_surface_collision_callback( TrajectorySurfaceCollisionCallback *tsur_cb ) {
	_tsur_cb = tsur_cb;
    }


    /*! \brief Set B-field potential dependent suppression callback.
     */
    void set_bfield_suppression_callback( const CallbackFunctorD_V *bsup_cb ) {
	_pidata.set_bfield_suppression_callback( bsup_cb );
    }

    /*! \brief Set relativistic particle iteration.
     */
    void set_relativistic( bool enable ) {
	_pidata.set_relativistic( enable );
    }

    /*! \brief Get particle iterator statistics.
     */
    const ParticleStatistics &get_statistics( void ) const {
	return( _stat );
    }

    /*! \brief Get accumulated time spent in ODE stepping (gsl_odeiv2_evolve_apply()).
     *
     *  Accumulated over all trajectories handled by this iterator
     *  instance since construction. Caller sums across all per-thread
     *  instances after Scheduler::finish().
     */
    double get_time_ode( void ) const {
	return( _time_ode );
    }

    /*! \brief Get accumulated time spent in handle_trajectory() (collision
     *  detection and space charge deposition combined).
     *
     *  Accumulated over all trajectories handled by this iterator
     *  instance since construction. Caller sums across all per-thread
     *  instances after Scheduler::finish().
     */
    double get_time_trajhandle( void ) const {
	return( _time_trajhandle );
    }


    /*! \brief Iterate a particle from start to end.
     *
     *  Iterate particle \a particle from start to end. This function
     *  is called by the Scheduler \a scheduler, which provides
     *  particles to be solved. Reference to \a scheduler is provided
     *  for the possibility to add secondary particles to particle
     *  database.
     */
    void operator()( Particle<PP> *particle, uint32_t pi ) {

	DEBUG_MESSAGE( "Calculating particle " << pi << "\n" );
	DEBUG_INC_INDENT();

	// Check particle status
	if( particle->get_status() != PARTICLE_OK ) {
	    DEBUG_DEC_INDENT();
	    return;
	}

	// Copy starting point to x and 
	PP x = particle->x();

	// Check particle definition
	if( !check_particle_definition( x ) ) {
	    particle->set_status( PARTICLE_BADDEF );
	    _stat.inc_end_baddef();
	    DEBUG_DEC_INDENT();
	    return;
	}
	particle->x() = x;

	// Reset trajectory and save first trajectory point.
	_traj.clear();
	save_trajectory_point( x );
	_pidata._qm = particle->qm();
	_xi = x;
	if( _scharge_dep == SCHARGE_DEPOSITION_LINEAR ) {
	    _cdpast.clear();
	    _cdpast.push( x );
	}

	// Reset integrator
	gsl_odeiv2_step_reset( _step );
	gsl_odeiv2_evolve_reset( _evolve );
	
	// Make initial guess for step size
	double dxdt[PP::size()-1];
	PP::get_derivatives( 0.0, &x[1], dxdt, (void *)&_pidata );
	double dt = calculate_dt( x, dxdt );

	DEBUG_MESSAGE( "Starting values for iteration:\n" << 
		       "  dxdt = " );
	for( size_t a = 0; a < PP::size()-1; a++ )
	    DEBUG_MESSAGE( dxdt[a] << " " );
	DEBUG_MESSAGE( "\n" <<
		       "  dt = " << dt << "\n" );
	DEBUG_MESSAGE( "Starting iteration\n" );
	
	// Iterate ODEs until maximum steps are done, time is used 
	// or particle collides.
	PP x2;
	size_t nstp = 0; // Steps taken
	while( nstp < _maxsteps && x[0] < _maxt ) {

	    // Take a step.
	    x2 = x;

	    DEBUG_MESSAGE( "\nStep *** *** *** ****\n" <<
			   "  x  = " << x2 << "\n" <<
			   "  dt = " << dt << " (proposed)\n" );

	    {
		Clock::time_point time_t0 = Clock::now();
		while( true ) {
		    int retval = gsl_odeiv2_evolve_apply( _evolve, _control, _step, &_system,
							 &x2[0], _maxt, &dt, &x2[1] );
		    if( retval == IBSIMU_DERIV_ERROR ) {
			DEBUG_MESSAGE( "Step rejected\n" <<
				       "  x2 = " << x2 << "\n" <<
				       "  dt = " << dt << "\n" );
			/* Restore the FULL state, not just the time.
			 *
			 * x2 = x is done once, above, OUTSIDE this retry loop.
			 * gsl_odeiv2_evolve_apply() is handed &x2[0] as the
			 * time and &x2[1] as the state, and on a rejected step
			 * it can leave the state partially advanced -- the
			 * original code here restored only x2[0], on the
			 * suspicion (see the note it carried about GSL-1.12)
			 * that the time was not being rolled back.
			 *
			 * The state has exactly the same problem and was not
			 * being rolled back at all. So a step that overshoots
			 * into an invalid region -- r <= 0 at the axis, say --
			 * leaves x2 sitting THERE, and every retry then starts
			 * from the invalid point instead of from x. The first
			 * RK stage is evaluated at the start of the step, so
			 * they all fail regardless of dt, and dt halves to zero:
			 * "too small step size", thrown for a particle that was
			 * never actually stuck.
			 *
			 * Restoring the whole state makes the retry mean what
			 * it says -- take the same step again, smaller.
			 */
			x2 = x;

			/* GSL zeroes *h on some failure paths, and 0*0.5 is
			 * still 0 -- which would look like an underflow on the
			 * FIRST rejection rather than after a genuine sequence
			 * of halvings. Re-seed from the geometric estimate in
			 * that case so the retry has something to work with,
			 * and let the counter below decide when to give up. */
			if( dt == 0.0 ) {
			    double dxdt_rs[PP::size()-1];
			    if( PP::get_derivatives( 0.0, &x[1], dxdt_rs,
						     (void *)&_pidata ) == GSL_SUCCESS )
				dt = 0.5*calculate_dt( x, dxdt_rs );
			} else {
			    dt *= 0.5;
			}
			if( dt == 0.0 ) {

			    /* STUCK PARTICLE -- kill it, do not abort the run.
			     *
			     * Halving only helps if the START of the step is a
			     * valid sampling point: the first RK stage is
			     * evaluated there. Once the current state itself
			     * fails get_derivatives()'s range test -- r <= 0,
			     * or outside the box by more than h -- every
			     * candidate step is rejected however small, and dt
			     * grinds to zero. It is a stuck state, not a step
			     * size problem, and no dt can rescue it.
			     *
			     * With the state restored properly above, this
			     * should now be unreachable in practice: a retry
			     * from a valid x with a smaller dt has somewhere to
			     * go. It is kept as a backstop, because throwing
			     * was a poor trade -- one pathological particle out
			     * of 192000 destroyed a 9000-second run at cycle 30
			     * with nothing kept -- and because a genuinely
			     * stuck state is still conceivable (a particle
			     * placed outside the domain by a mirror, for
			     * instance). If this fires, the state was already
			     * bad on ENTRY to the step, which is a different
			     * bug and the printed location says where.
			     * Killing it as PARTICLE_BADDEF loses that
			     * particle's contribution to the space charge --
			     * one part in 1e5 -- and the run continues.
			     *
			     * Reported with its coordinates, rate-limited,
			     * because WHERE it gets stuck is the diagnosis: r
			     * at or below zero means the axis handling, x
			     * beyond the box means a mirror or reflection put
			     * it there. A steady trickle is tolerable; a burst
			     * means something upstream is placing particles
			     * outside the domain and should be fixed rather
			     * than absorbed here.
			     */
			    if( _stuck_reported < 10 ) {
				_stuck_reported++;
				ibsimu.message( MSG_WARNING, 1 )
				    << "Warning: particle stuck (step size underflow) at "
				    << x2.location() << ", killed as BADDEF."
				    << ( _stuck_reported == 10
					 ? " Further occurrences on this thread suppressed.\n"
					 : "\n" );
			    }
			    particle->set_status( PARTICLE_BADDEF );
			    _stat.inc_end_baddef();
			    DEBUG_DEC_INDENT();
			    return;
			}
			//nstp++;
			continue;
		    } else if( retval == GSL_SUCCESS ) {
			break;
		    } else {

			/* A THIRD return value: neither success nor our own
			 * IBSIMU_DERIV_ERROR (201). GSL has its own failure
			 * modes here -- notably it can decide the step size
			 * cannot be reduced further and return GSL_ENOPROG or
			 * GSL_FAILURE, and on some paths it ZEROES *h before
			 * returning, which is why a run can reach this branch
			 * rather than the underflow one above.
			 *
			 * This used to throw, discarding the code without ever
			 * saying which one it was -- so the message named the
			 * symptom and hid the diagnosis. It now prints the code
			 * and gsl_strerror(), and kills the particle instead of
			 * the run, consistent with the other two exits.
			 *
			 * If these appear in numbers, the code is the thing to
			 * report: ENOPROG means GSL gave up on the step size
			 * (same pathology as the underflow path, detected
			 * inside GSL), while EBADFUNC would mean the derivative
			 * function is returning something GSL cannot interpret,
			 * which would be a different bug entirely.
			 */
			if( _stuck_reported < 10 ) {
			    _stuck_reported++;
			    ibsimu.message( MSG_WARNING, 1 )
				<< "Warning: gsl_odeiv2_evolve_apply returned "
				<< retval << " (" << gsl_strerror( retval )
				<< ") for particle " << pi << " at "
				<< x2.location() << ", killed as BADDEF."
				<< ( _stuck_reported == 10
				     ? " Further occurrences on this thread suppressed.\n"
				     : "\n" );
			}
			particle->set_status( PARTICLE_BADDEF );
			_stat.inc_end_baddef();
			DEBUG_DEC_INDENT();
			return;
		    }
		}
		_time_ode += std::chrono::duration<double>( Clock::now()-time_t0 ).count();
	    }
	    
	    // Check step count number and step size validity
	    if( nstp >= _maxsteps )
	        break;
	    if( x2[0] == x[0] ) {

		/* The step was ACCEPTED but advanced time by zero, so the
		 * particle cannot progress. Kill it rather than aborting the
		 * run -- same reasoning as the underflow backstop above, and
		 * the same trade: losing one particle's space charge beats
		 * losing every completed cycle.
		 *
		 * The original also dumped the entire trajectory to the
		 * message stream. That is worse than useless here: particle
		 * iteration is threaded, so the dump interleaves with 31
		 * other threads, and a long trajectory is thousands of lines.
		 * One line with the location is the part that identifies
		 * where it happened.
		 */
		if( _stuck_reported < 10 ) {
		    _stuck_reported++;
		    ibsimu.message( MSG_WARNING, 1 )
			<< "Warning: particle " << pi << " made no progress (dt "
			<< "accepted as zero) at " << x2.location()
			<< ", killed as BADDEF."
			<< ( _stuck_reported == 10
			     ? " Further occurrences on this thread suppressed.\n"
			     : "\n" );
		}
		particle->set_status( PARTICLE_BADDEF );
		_stat.inc_end_baddef();
		DEBUG_DEC_INDENT();
		return;
	    }
	    
	    // Increase step count.
	    nstp++;

	    DEBUG_MESSAGE( "Step accepted from x1 to x2:\n" <<
			   "  x1 = " << x << "\n" <<
			   "  x2 = " << x2 << "\n" );

	    // Handle collisions and space charge of step.
	    {
		Clock::time_point time_t0 = Clock::now();
		bool ok = handle_trajectory( *particle, x, x2, false, nstp == 0 );
		_time_trajhandle += std::chrono::duration<double>( Clock::now()-time_t0 ).count();
		if( !ok ) {
		    x = x2;
		    break; // Particle done
		}
	    }

	    // Check if particle mirroring is required to avoid 
	    // singularity at symmetry axis.
	    if( axis_mirror_required( x2 ) ) {
		if( !handle_axis_mirror_step( *particle, x, x2 ) )
		    break; // Particle done
	    }

	    // Propagate coordinates
	    x = x2;

	    // Save trajectory point
	    save_trajectory_point( x2 );
	}

	DEBUG_MESSAGE( "Done iterating\n" );

	// Check if step count or time limited 
	if( nstp == _maxsteps ) {
	    particle->set_status( PARTICLE_NSTP );
	    _stat.inc_end_step();
	} else if( x[0] >= _maxt ) {
	    particle->set_status( PARTICLE_TIME );
	    _stat.inc_end_time();
	}

	// Save step count
	_stat.inc_sum_steps( nstp );

	// Save trajectory of current particle
	if( _trajdiv != 0 && pi % _trajdiv == 0 )
	    particle->copy_trajectory( _traj );

	// Save last particle location
	particle->x() = x;

	// Call trajectory end callback
	if( _tend_cb )
	    (*_tend_cb)( particle, _pdb );

	DEBUG_DEC_INDENT();
    }

};


#endif

