/*! \file dxf_solid.cpp
 *  \brief %Solid definition using MyDXF
 */

/* Copyright (c) 2010-2012,2014 Taneli Kalvas. All rights reserved.
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


#include <iostream>
#include <sstream>
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "compmath.hpp"
#include "ibsimu.hpp"


DXFSolid::DXFSolid( MyDXFFile *dxffile, const std::string &layername )
    : _func(&unity), _has_bbox(false)
{
    ibsimu.message( 1 ) << "Defining electrode \'" << layername << "\'\n";

    MyDXFEntities *ent = dxffile->get_entities();
    MyDXFEntitySelection *layer = ent->selection_layer( layername );
    if( layer->size() == 0 ) {
	delete layer;
	std::stringstream se;
	se << "No entities in layer\n. Layers in dxf-file:\n";
	
	// Add list of layer names
	MyDXFTables *tables = dxffile->get_tables();
	MyDXFTable *layers = tables->get_layers_table();
	for( size_t a = 0; a < layers->size(); a++ ) {
	    MyDXFTableEntry *entry = layers->get_entry( a );
	    MyDXFTableEntryLayer *layer = dynamic_cast<MyDXFTableEntryLayer *>( entry );
	    se << a << ": \'" << layer->name() << "\'\n";
	}

	throw( Error( ERROR_LOCATION, se.str() ) );
    }
    MyDXFEntitySelection *loop = ent->selection_path_loop( layer );
    if( loop->size() == 0 ) {
	delete layer;
	delete loop;
	throw( Error( ERROR_LOCATION, "No loops defined in layer" ) );
    }

    _entities = new MyDXFEntities( dxffile, ent, loop );
    _selection = _entities->selection_all();

    if( (int)layer->size()-(int)loop->size() > 0 )
	ibsimu.message( 1 ) << "  removed "
			    << (int)layer->size()-(int)loop->size()
			    << " entities\n";
    ibsimu.message( 1 ) << "  solid defined using "
			<< _entities->size() << " entities\n";

    // Cache the 2D extent of the loop now, while dxffile is still
    // guaranteed to be alive -- the class intentionally keeps no
    // dependency on dxffile after construction.
    Transformation ident;
    _entities->get_bbox( _selection, _bbox_min, _bbox_max, dxffile, &ident );
    _has_bbox = true;

    delete layer;
    delete loop;
}


DXFSolid::DXFSolid( MyDXFFile *dxffile, MyDXFEntities *ent )
    : _func(&unity), _has_bbox(false)
{
    if( ent->size() == 0 ) {
	throw( Error( ERROR_LOCATION, "No entities" ) );
    }

    MyDXFEntitySelection *all = _entities->selection_all();
    MyDXFEntitySelection *loop = ent->selection_path_loop( all );
    if( loop->size() == 0 ) {
	delete all;
	delete loop;
	throw( Error( ERROR_LOCATION, "No loops defined by entities" ) );
    }

    _entities = new MyDXFEntities( dxffile, ent, loop );
    _selection = _entities->selection_all();

    if( (int)all->size()-(int)loop->size() > 0 )
	ibsimu.message( 1 ) << "  removed "
			    << (int)all->size()-(int)loop->size()
			    << " entities\n";
    ibsimu.message( 1 ) << "  solid defined using "
			<< _entities->size() << " entities\n";

    // Cache the 2D extent of the loop now, while dxffile is still
    // guaranteed to be alive -- the class intentionally keeps no
    // dependency on dxffile after construction.
    Transformation ident;
    _entities->get_bbox( _selection, _bbox_min, _bbox_max, dxffile, &ident );
    _has_bbox = true;

    delete all;
    delete loop;
}


DXFSolid::DXFSolid( std::istream &is )
    : _func(&unity), _entities(NULL), _selection(NULL), _has_bbox(false)
{
    ibsimu.message( MSG_WARNING, 1 ) << "Warning: loading of DXFSolid not implemented\n";
    ibsimu.flush();
}


DXFSolid::~DXFSolid()
{
    if( _entities )
	delete _entities;
    if( _selection )
	delete _selection;
}


Vec3D DXFSolid::unity( const Vec3D &x )
{
    return( x );
}


Vec3D DXFSolid::rotx( const Vec3D &x )
{
    return( Vec3D( x[0], sqrt( x[1]*x[1] + x[2]*x[2] ) ) );
}


Vec3D DXFSolid::roty( const Vec3D &x )
{
    return( Vec3D( x[1], sqrt( x[0]*x[0] + x[2]*x[2] ) ) );
}


Vec3D DXFSolid::rotz( const Vec3D &x )
{
    return( Vec3D( x[2], sqrt( x[0]*x[0] + x[1]*x[1] ) ) );
}


void DXFSolid::define_2x3_mapping( Vec3D (*func)(const Vec3D &) )
{
    if( func )
	_func = func;
    else
	_func = unity;
}


bool DXFSolid::inside( const Vec3D &x ) const
{
    if( !_entities )
	return( false );

    // Transform 3D -> 3D
    // T is direct transformation
    Vec3D y = _T.transform_point( x );

    // Transform 3D -> 2D
    Vec3D z = _func( y );

    if( isnan(z[0]) || isnan(z[1]) )
	return( true );
    else if( isinf(z[0]) || isinf(z[1]) )
	return( false );

    return( _entities->inside_loop( _selection, z[0], z[1] ) );
}


bool DXFSolid::get_bbox( Vec3D &min, Vec3D &max ) const
{
    if( !_entities || !_has_bbox )
	return( false );

    double xlo = _bbox_min[0], xhi = _bbox_max[0];
    double ylo = _bbox_min[1], yhi = _bbox_max[1];

    // Local (pre-_T) box implied by the current 2D<-3D mapping. Only
    // the built-in mappings are handled -- a user-supplied _func can
    // return NaN to force "inside" for arbitrary input, which makes
    // it impossible to infer any bound from the 2D loop extent alone.
    Vec3D lo, hi;
    if( _func == &unity ) {
	// (x,y) come straight from the dxf loop; z is the extrusion
	// direction and is not constrained at all.
	lo = Vec3D( xlo, ylo, -BBOX_HUGE );
	hi = Vec3D( xhi, yhi,  BBOX_HUGE );
    } else if( _func == &rotx ) {
	// x comes straight from the loop; (y,z) lie on a disc of
	// radius up to max(|ylo|,|yhi|) around the x-axis.
	double rmax = fabs(ylo) > fabs(yhi) ? fabs(ylo) : fabs(yhi);
	lo = Vec3D( xlo, -rmax, -rmax );
	hi = Vec3D( xhi,  rmax,  rmax );
    } else if( _func == &roty ) {
	double rmax = fabs(ylo) > fabs(yhi) ? fabs(ylo) : fabs(yhi);
	lo = Vec3D( -rmax, xlo, -rmax );
	hi = Vec3D(  rmax, xhi,  rmax );
    } else if( _func == &rotz ) {
	double rmax = fabs(ylo) > fabs(yhi) ? fabs(ylo) : fabs(yhi);
	lo = Vec3D( -rmax, -rmax, xlo );
	hi = Vec3D(  rmax,  rmax, xhi );
    } else {
	// Custom user-defined mapping -- no safe bound available.
	return( false );
    }

    return( bbox_from_local_box( lo, hi, min, max ) );
}


void DXFSolid::debug_print( std::ostream &os ) const
{
    os << "**DXFSolid\n";

    if( !_entities )
	return;

    for( size_t a = 0; a < _entities->size(); a++ ) {
	os << _entities->get_entity( a );
    }
}


void DXFSolid::save( std::ostream &os ) const
{
    write_int32( os, FILEID_DXFSOLID );
    ibsimu.message( MSG_WARNING, 1 ) << "Warning: saving of DXFSolid not implemented\n";
    ibsimu.flush();
}

