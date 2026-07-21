/*! \file solid.hpp
 *  \brief Base for solid definition
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include "solid.hpp"


const double Solid::BBOX_HUGE = 1.0e30;


Solid::Solid()
{

}


Solid::~Solid()
{

}


bool Solid::bbox_from_local_box( const Vec3D &lo, const Vec3D &hi, Vec3D &wmin, Vec3D &wmax ) const
{
    Transformation Tinv = _T.inverse();

    Vec3D corner[8] = {
	Vec3D( lo[0], lo[1], lo[2] ),
	Vec3D( hi[0], lo[1], lo[2] ),
	Vec3D( lo[0], hi[1], lo[2] ),
	Vec3D( hi[0], hi[1], lo[2] ),
	Vec3D( lo[0], lo[1], hi[2] ),
	Vec3D( hi[0], lo[1], hi[2] ),
	Vec3D( lo[0], hi[1], hi[2] ),
	Vec3D( hi[0], hi[1], hi[2] )
    };

    wmin = Vec3D(  BBOX_HUGE,  BBOX_HUGE,  BBOX_HUGE );
    wmax = Vec3D( -BBOX_HUGE, -BBOX_HUGE, -BBOX_HUGE );

    for( int c = 0; c < 8; c++ ) {
	Vec3D w = Tinv.transform_point( corner[c] );
	for( int a = 0; a < 3; a++ ) {
	    if( w[a] < wmin[a] )
		wmin[a] = w[a];
	    if( w[a] > wmax[a] )
		wmax[a] = w[a];
	}
    }

    // A local axis flagged with +-BBOX_HUGE stays huge (or becomes
    // huge scaled by whatever _T does to it) after an affine
    // transform unless it happens to land on a matrix row that zeros
    // it out exactly -- and a genuinely bounded local axis can never
    // produce a huge world coordinate for any physically sensible
    // mesh/solid. So any world axis found here at a "mesh-implausible"
    // magnitude is not a real bound: flag it as unrestricted instead
    // of handing the caller a meaningless huge number.
    bool any_finite = false;
    const double sanity_limit = 1.0e10;
    for( int a = 0; a < 3; a++ ) {
	if( fabs(wmin[a]) >= sanity_limit || fabs(wmax[a]) >= sanity_limit ) {
	    wmin[a] = -BBOX_HUGE;
	    wmax[a] =  BBOX_HUGE;
	} else {
	    any_finite = true;
	}
    }

    return( any_finite );
}


void Solid::reset_transformation( void ) 
{
    _T.reset();
}


void Solid::set_transformation( const Transformation &T ) 
{
    // T is inverse of the matrix presented to user
    _T = T.inverse();
}


void Solid::translate( const Vec3D &dx ) 
{
    _T = _T * Transformation::translation( -dx );
}


void Solid::scale( double sx ) 
{
    _T = _T * Transformation::scaling( Vec3D(1.0/sx, 1.0/sx, 1.0/sx) );
}


void Solid::scale( const Vec3D &sx ) 
{
    _T = _T * Transformation::scaling( Vec3D(1.0/sx[0], 1.0/sx[1], 1.0/sx[2]) );
}


void Solid::rotate_x( double a ) 
{
    _T = _T * Transformation::rotation_x( -a );
}


void Solid::rotate_y( double a ) 
{
    _T = _T * Transformation::rotation_y( -a );
}


void Solid::rotate_z( double a ) 
{
    _T = _T * Transformation::rotation_z( -a );
}


