/*! \file vtriangle.hpp
 *  \brief Vertex-based triangle representation
 */

/* Copyright (c) 2011-2013 Taneli Kalvas. All rights reserved.
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


#include <iomanip>
#include <limits>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "vtriangle.hpp"
#include "ibsimu.hpp"


VTriangle::VTriangle( uint32_t v1, uint32_t v2, uint32_t v3 )
{
    _v[0] = v1;
    _v[1] = v2;
    _v[2] = v3;
}


VTriangle::VTriangle( const uint32_t v[3] )
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
}


VTriangle::~VTriangle()
{
    
}


void VTriangle::debug_print( std::ostream &os ) const
{
    os << "**VTriangle\n";    
    os << "  v = "  
       << std::setw(6) << _v[0] << " "
       << std::setw(6) << _v[1] << " "
       << std::setw(6) << _v[2] << "\n";
}


/*
 * VTriangleSurface
 * ****************************************************** */


VTriangleSurface::VTriangleSurface( double vertex_matching_eps )
    : _vertex_matching_eps(vertex_matching_eps*vertex_matching_eps)
{

}


VTriangleSurface::~VTriangleSurface()
{

}


void VTriangleSurface::set_vertex_matching_eps( double vertex_matching_eps )
{
    _vertex_matching_eps = vertex_matching_eps*vertex_matching_eps;
}


void VTriangleSurface::clear( void )
{
    _vertex.clear();
    _triangle.clear();
}


uint32_t VTriangleSurface::add_vertex( const Vec3D &x )
{
    /* Search for vertex matching to x, allow equal to matching
     * epsilon to allow the use of zero as vertex_matching_eps.
    */
    uint32_t v;
    for( v = 0; v < _vertex.size(); v++ )
	if( ssqr(x-_vertex[v]) <= _vertex_matching_eps )
	    return( v ); // Vertex found

    // New vertex needed
    _vertex.push_back( x ); 
    return( _vertex.size()-1 );    
}


uint32_t VTriangleSurface::add_triangle( const Vec3D x[3] )
{
    uint32_t v[3];

    for( int a = 0; a < 3; a++ )
	v[a] = add_vertex( x[a] );

    _triangle.push_back( VTriangle( v ) );    
    return( _triangle.size()-1 );
}


uint32_t VTriangleSurface::add_triangle( const Vec3D &x1, const Vec3D &x2, const Vec3D &x3 )
{
    uint32_t v[3];

    v[0] = add_vertex( x1 );
    v[1] = add_vertex( x2 );
    v[2] = add_vertex( x3 );
    
    _triangle.push_back( VTriangle( v ) );
    return( _triangle.size()-1 );
}


void VTriangleSurface::debug_print( std::ostream &os ) const
{
    os << "**VTriangleSurface\n";

    os << "vertex_matching_eps = " << _vertex_matching_eps << "\n";
    os << "trianglec = " << _triangle.size() << "\n";
    os << "vertexc = " << _vertex.size() << "\n";

    for( uint32_t a = 0; a < _triangle.size(); a++ ) 
	os << "  triangle[" << a << "] = (" 
	   << _triangle[a][0] << ", "
	   << _triangle[a][1] << ", "
	   << _triangle[a][2] << ")\n";

    for( uint32_t a = 0; a < _vertex.size(); a++ ) 
	os << "  vertex[" << a << "] = " << _vertex[a] << "\n";
}


/*
 * VTriangleSurfaceSolid
 * ****************************************************** */


VTriangleSurfaceSolid::VTriangleSurfaceSolid( double vertex_matching_eps, 
					      double signed_volume_eps )
    : VTriangleSurface(vertex_matching_eps), 
      _signed_volume_eps(signed_volume_eps)
{

}


VTriangleSurfaceSolid::~VTriangleSurfaceSolid()
{

}


#define VTRI_EDGE1   0
#define VTRI_EDGE2   1
#define VTRI_EDGE3   2
#define VTRI_FACE    3
#define VTRI_OUTSIDE 4
#define VTRI_INSIDE  5


int VTriangleSurfaceSolid::signvol4( const Vec3D &q0, const Vec3D &q1, 
				     const Vec3D &q2, const Vec3D &q3 ) const
{
    double x = q0[0]*(-q1[1]*q2[2] + q1[1]*q3[2] + q2[1]*q1[2] - 
		       q2[1]*q3[2] - q3[1]*q1[2] + q3[1]*q2[2] ) +
	       q1[0]*( q0[1]*q2[2] - q0[1]*q3[2] - q2[1]*q0[2] + 
		       q2[1]*q3[2] + q3[1]*q0[2] - q3[1]*q2[2] ) +
	       q2[0]*(-q0[1]*q1[2] + q0[1]*q3[2] + q1[1]*q0[2] - 
		       q1[1]*q3[2] - q3[1]*q0[2] + q3[1]*q1[2] ) +
	       q3[0]*( q0[1]*q1[2] - q0[1]*q2[2] - q1[1]*q0[2] + 
		       q1[1]*q2[2] + q2[1]*q0[2] - q2[1]*q1[2] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


int VTriangleSurfaceSolid::signvol3( const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const
{
    double x = q1[0]*( q2[1]*q3[2] - q3[1]*q2[2] ) -
	       q1[1]*( q2[0]*q3[2] - q3[0]*q2[2] ) +
	       q1[2]*( q2[0]*q3[1] - q3[0]*q2[1] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


/* Classify inclusion of point p in tetrahedron (o,q1,q2,q3), when
 * The sense of tetrahedron (o,q1,q2,q3) is ss.
 */
int VTriangleSurfaceSolid::classify_original_tetrahedron( int ss, const Vec3D &p, 
							  const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 ) const
{
    int s0, s1, s2, s3;

    if( ss > 0 ) {
	// Positive sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) < 0 )
	    return( VTRI_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) < 0 )
	    return( VTRI_OUTSIDE );
    } else {
	// Negative sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) > 0 )
	    return( VTRI_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) > 0 )
	    return( VTRI_OUTSIDE );
    }

    // Edges and faces
    if( s1 == 0 ) {
	if( s2 == 0 )
	    return( VTRI_EDGE3 );
	else if( s3 == 0 )
	    return( VTRI_EDGE2 );
	return( VTRI_FACE );
    } else if( s2 == 0 ) {
	if( s3 == 0 )
	    return( VTRI_EDGE1 );
	return( VTRI_FACE );
    } else if( s3 == 0 ) {
	return( VTRI_FACE );
    }

    return( VTRI_INSIDE );
}


/* Exhaustive form: the signed tetrahedron decomposition over EVERY
 * triangle. Kept as written, both as the fallback when no grid has been
 * built and as the reference the accelerated path must reproduce exactly.
 */
bool VTriangleSurfaceSolid::inside_exhaustive( const Vec3D &x ) const
{
    // Cleared positive and negative vertex arrays
    std::vector<bool> vpos( _vertex.size(), false );
    std::vector<bool> vneg( _vertex.size(), false );

    int incl = 0;
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {

	int ss = signvol3( _vertex[_triangle[a][0]],
			   _vertex[_triangle[a][1]],
			   _vertex[_triangle[a][2]] );
	int stat = classify_original_tetrahedron( ss, x,
						  _vertex[_triangle[a][0]],
						  _vertex[_triangle[a][1]],
						  _vertex[_triangle[a][2]] );
	if( stat == VTRI_INSIDE ) {
	    incl += 2*ss;
	} else if( stat == VTRI_FACE ) {
	    incl += ss;
	} else if( stat != VTRI_OUTSIDE ) {
	    if( ss > 0 && !vpos[_triangle[a][stat]] ) {
		vpos[_triangle[a][stat]] = true;
		incl += 2*ss;
	    } else if( ss < 0 && !vneg[_triangle[a][stat]] ) {
		vneg[_triangle[a][stat]] = true;
		incl += 2*ss;
	    }
	}
    }

    if( incl > 0 )
	return( true ); 
    return( false );
}


void VTriangleSurfaceSolid::update_bbox( Vec3D &min, Vec3D &max, const Vec3D x ) const
{
    for( int a = 0; a < 3; a++ ) {
	if( x[a] < min[a] )
	    min[a] = x[a];
	if( x[a] > max[a] )
	    max[a] = x[a];
    }
}


/* Build a uniform grid holding, for each cell, the triangles whose
 * bounding box overlaps it. Cells are sized so the average occupancy is
 * about one triangle.
 */
/* Exact triangle / axis-aligned-box overlap (Akenine-Moller separating
 * axis test). Used to decide which grid cells a triangle really touches.
 *
 * Inserting by BOUNDING BOX instead is far simpler and was tried first,
 * but these surfaces are CAD tessellations of annular faces: a flat disc
 * from r = 6.5 to r = 48 mm comes out as long radial slivers whose bbox
 * spans most of the grid while the triangle itself touches a line of
 * cells. Bbox insertion put such a triangle in ~1600 cells instead of
 * ~60, which inflated both the table and every candidate list drawn from
 * it, and left the query only 2.6x faster than the exhaustive sweep.
 */
static bool plane_box_overlap( const Vec3D &normal, const Vec3D &vert,
			       const Vec3D &maxbox )
{
    Vec3D vmin, vmax;
    for( int q = 0; q < 3; q++ ) {
	if( normal[q] > 0.0 ) {
	    vmin[q] = -maxbox[q] - vert[q];
	    vmax[q] =  maxbox[q] - vert[q];
	} else {
	    vmin[q] =  maxbox[q] - vert[q];
	    vmax[q] = -maxbox[q] - vert[q];
	}
    }
    if( normal*vmin > 0.0 ) return( false );
    if( normal*vmax >= 0.0 ) return( true );
    return( false );
}

static bool tri_box_overlap( const Vec3D &boxcenter, const Vec3D &boxhalf,
			     const Vec3D &t0, const Vec3D &t1, const Vec3D &t2 )
{
    Vec3D v0 = t0-boxcenter, v1 = t1-boxcenter, v2 = t2-boxcenter;
    Vec3D e0 = v1-v0, e1 = v2-v1, e2 = v0-v2;

#define AXISTEST(a0,a1,b0,b1,c0,c1,ea,eb)				    do {									double p0 = (a0)*(b0) - (a1)*(b1);					double p1 = (a0)*(c0) - (a1)*(c1);					double mn = p0 < p1 ? p0 : p1, mx = p0 < p1 ? p1 : p0;			double rad = fabs(a0)*boxhalf[ea] + fabs(a1)*boxhalf[eb];		if( mn > rad || mx < -rad ) return( false );			    } while( 0 )

    AXISTEST( e0[2], e0[1], v0[1], v0[2], v2[1], v2[2], 1, 2 );
    AXISTEST( e0[2], e0[0], v0[0], v0[2], v2[0], v2[2], 0, 2 );
    AXISTEST( e0[1], e0[0], v1[0], v1[1], v2[0], v2[1], 0, 1 );
    AXISTEST( e1[2], e1[1], v0[1], v0[2], v2[1], v2[2], 1, 2 );
    AXISTEST( e1[2], e1[0], v0[0], v0[2], v2[0], v2[2], 0, 2 );
    AXISTEST( e1[1], e1[0], v0[0], v0[1], v1[0], v1[1], 0, 1 );
    AXISTEST( e2[2], e2[1], v0[1], v0[2], v1[1], v1[2], 1, 2 );
    AXISTEST( e2[2], e2[0], v0[0], v0[2], v1[0], v1[2], 0, 2 );
    AXISTEST( e2[1], e2[0], v1[0], v1[1], v2[0], v2[1], 0, 1 );
#undef AXISTEST

    for( int q = 0; q < 3; q++ ) {
	double mn = v0[q], mx = v0[q];
	if( v1[q] < mn ) mn = v1[q];
	if( v1[q] > mx ) mx = v1[q];
	if( v2[q] < mn ) mn = v2[q];
	if( v2[q] > mx ) mx = v2[q];
	if( mn > boxhalf[q] || mx < -boxhalf[q] ) return( false );
    }

    Vec3D normal = cross( e0, e1 );
    return( plane_box_overlap( normal, v0, boxhalf ) );
}


void VTriangleSurfaceSolid::build_grid( void )
{
    _gstart.clear();
    _gtri.clear();
    _cand.clear();

    const size_t ntri = _triangle.size();
    if( ntri == 0 )
	return;

    // Grid spans the SHIFTED vertices, i.e. the frame inside() works in.
    Vec3D gmin(  std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    Vec3D gmax( -std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity() );
    for( size_t a = 0; a < _vertex.size(); a++ )
	update_bbox( gmin, gmax, _vertex[a] );

    Vec3D ext = gmax-gmin;
    for( int a = 0; a < 3; a++ )
	if( ext[a] <= 0.0 )
	    ext[a] = 1.0;
    double cell = pow( ext[0]*ext[1]*ext[2]/(double)ntri, 1.0/3.0 );
    if( cell <= 0.0 )
	return;
    int64_t ncell = 1;
    for( int a = 0; a < 3; a++ ) {
	_gn[a] = (int32_t)floor( ext[a]/cell );
	if( _gn[a] < 1 )   _gn[a] = 1;
	if( _gn[a] > 128 ) _gn[a] = 128;
	_gh[a] = ext[a]/_gn[a];
	ncell *= _gn[a];
    }
    _gmin = gmin;

    // Two passes: count per cell, prefix sum, then scatter.
    std::vector<uint32_t> count( (size_t)ncell+1, 0 );
    std::vector<int32_t> lo( 3*ntri ), hi( 3*ntri );
    for( size_t t = 0; t < ntri; t++ ) {
	Vec3D tmin(  std::numeric_limits<double>::infinity(),
		     std::numeric_limits<double>::infinity(),
		     std::numeric_limits<double>::infinity() );
	Vec3D tmax( -std::numeric_limits<double>::infinity(),
		    -std::numeric_limits<double>::infinity(),
		    -std::numeric_limits<double>::infinity() );
	for( int v = 0; v < 3; v++ )
	    update_bbox( tmin, tmax, _vertex[_triangle[t][v]] );
	for( int a = 0; a < 3; a++ ) {
	    int32_t l = (int32_t)floor( (tmin[a]-_gmin[a])/_gh[a] );
	    int32_t h = (int32_t)floor( (tmax[a]-_gmin[a])/_gh[a] );
	    /* Clamp BOTH ends. A triangle touching the grid's upper edge
	     * gives floor()==_gn[a] for the low index too, and the
	     * "if( h < l ) h = l" below would then carry that out-of-range
	     * value straight into the cell index. */
	    if( l < 0 ) l = 0;
	    if( l > _gn[a]-1 ) l = _gn[a]-1;
	    if( h < 0 ) h = 0;
	    if( h > _gn[a]-1 ) h = _gn[a]-1;
	    if( h < l ) h = l;
	    lo[3*t+a] = l;
	    hi[3*t+a] = h;
	}
	const Vec3D &q0 = _vertex[_triangle[t][0]];
	const Vec3D &q1 = _vertex[_triangle[t][1]];
	const Vec3D &q2 = _vertex[_triangle[t][2]];
	/* Inflate the test box by a hair. Vertices of an axis-aligned face
	 * can land exactly on a cell boundary -- the flat faces of these
	 * electrodes sit exactly on the grid's outer plane -- and the
	 * separating-axis test then rejects or accepts on the last bit of
	 * the box-centre arithmetic. Growing the box can only ever add a
	 * triangle to a cell, which costs a little work and cannot lose a
	 * contribution; shrinking it silently drops surface. */
	const double sat_eps = 1.0e-9;
	const Vec3D half( 0.5*_gh[0]*(1.0+sat_eps) + sat_eps*_gh[0],
			  0.5*_gh[1]*(1.0+sat_eps) + sat_eps*_gh[1],
			  0.5*_gh[2]*(1.0+sat_eps) + sat_eps*_gh[2] );
	for( int32_t k = lo[3*t+2]; k <= hi[3*t+2]; k++ )
	    for( int32_t j = lo[3*t+1]; j <= hi[3*t+1]; j++ )
		for( int32_t i = lo[3*t+0]; i <= hi[3*t+0]; i++ ) {
		    Vec3D ctr( _gmin[0]+(i+0.5)*_gh[0],
			       _gmin[1]+(j+0.5)*_gh[1],
			       _gmin[2]+(k+0.5)*_gh[2] );
		    if( tri_box_overlap( ctr, half, q0, q1, q2 ) )
			count[(size_t)((k*_gn[1]+j)*_gn[0]+i)+1]++;
		}
    }
    for( size_t c = 1; c <= (size_t)ncell; c++ )
	count[c] += count[c-1];
    _gstart = count;
    _gtri.resize( count[ncell] );
    std::vector<uint32_t> fill( _gstart.begin(), _gstart.end()-1 );
    for( size_t t = 0; t < ntri; t++ ) {
	const Vec3D &q0 = _vertex[_triangle[t][0]];
	const Vec3D &q1 = _vertex[_triangle[t][1]];
	const Vec3D &q2 = _vertex[_triangle[t][2]];
	/* Inflate the test box by a hair. Vertices of an axis-aligned face
	 * can land exactly on a cell boundary -- the flat faces of these
	 * electrodes sit exactly on the grid's outer plane -- and the
	 * separating-axis test then rejects or accepts on the last bit of
	 * the box-centre arithmetic. Growing the box can only ever add a
	 * triangle to a cell, which costs a little work and cannot lose a
	 * contribution; shrinking it silently drops surface. */
	const double sat_eps = 1.0e-9;
	const Vec3D half( 0.5*_gh[0]*(1.0+sat_eps) + sat_eps*_gh[0],
			  0.5*_gh[1]*(1.0+sat_eps) + sat_eps*_gh[1],
			  0.5*_gh[2]*(1.0+sat_eps) + sat_eps*_gh[2] );
	for( int32_t k = lo[3*t+2]; k <= hi[3*t+2]; k++ )
	    for( int32_t j = lo[3*t+1]; j <= hi[3*t+1]; j++ )
		for( int32_t i = lo[3*t+0]; i <= hi[3*t+0]; i++ ) {
		    Vec3D ctr( _gmin[0]+(i+0.5)*_gh[0],
			       _gmin[1]+(j+0.5)*_gh[1],
			       _gmin[2]+(k+0.5)*_gh[2] );
		    if( tri_box_overlap( ctr, half, q0, q1, q2 ) )
			_gtri[fill[(size_t)((k*_gn[1]+j)*_gn[0]+i)]++] = (uint32_t)t;
		}
    }
    ibsimu.message( 1 ) << "  inside() grid: " << _gn[0] << "x" << _gn[1]
			<< "x" << _gn[2] << " cells, " << _gtri.size()
			<< " triangle refs for " << ntri << " triangles\n";

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    _cand.resize( nthreads > 0 ? nthreads : 1 );
}


/* Point-in-solid, accelerated but ARITHMETICALLY IDENTICAL to
 * inside_exhaustive().
 *
 * The exhaustive form sums a signed contribution over every triangle,
 * where triangle t contributes only if the query point lies inside the
 * tetrahedron (O, t) -- O being the origin of the shifted frame, which
 * prepare_for_inside() places just outside the solid's minimum corner.
 * Every other triangle returns VTRI_OUTSIDE and contributes exactly zero.
 *
 * The point lies inside tet(O, t) precisely when the ray leaving x along
 * the direction x -- i.e. straight away from O -- strikes triangle t. So
 * the triangles that can contribute anything at all are exactly those the
 * ray meets, and they can be gathered by walking the grid along that ray
 * instead of by scanning the whole list. Everything skipped would have
 * added zero.
 *
 * Two details keep it bit-exact rather than merely equivalent:
 *
 *   - the gathered candidates are SORTED back into triangle order before
 *     they are accumulated. The vpos/vneg de-duplication is order
 *     dependent -- for a vertex shared by several triangles the first one
 *     seen claims it -- so visiting in grid order rather than index order
 *     would be a different, equally defensible answer, and we want the
 *     same one.
 *   - a triangle spanning several cells appears in each, so duplicates are
 *     removed after sorting; counting it twice would double its term.
 */
bool VTriangleSurfaceSolid::inside( const Vec3D &p ) const
{
    // Fast test if outside bbox
    for( uint32_t a = 0; a < 3; a++ ) {
	if( p[a] < _bbox[0][a] )
	    return( false );
	if( p[a] > _bbox[1][a] )
	    return( false );
    }

    // Offset point
    Vec3D x = p+_offset;

    if( _gstart.empty() || _cand.empty() )
	return( inside_exhaustive( x ) );

    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if( tid < 0 || tid >= (int)_cand.size() )
	return( inside_exhaustive( x ) );
    std::vector<uint32_t> &cand = _cand[tid];
    cand.clear();

    /* Walk the grid from x along direction d = x (away from O) with a 3D
     * DDA, collecting every triangle in every cell entered. Conservative:
     * a cell is visited whenever the ray passes through it, so no triangle
     * the ray actually hits can be missed. */
    Vec3D d = x;
    int32_t ijk[3], step[3];
    double tmax[3], tdelta[3];
    for( int a = 0; a < 3; a++ ) {
	ijk[a] = (int32_t)floor( (x[a]-_gmin[a])/_gh[a] );
	if( ijk[a] < 0 ) ijk[a] = 0;
	if( ijk[a] > _gn[a]-1 ) ijk[a] = _gn[a]-1;
	if( d[a] > 0.0 ) {
	    step[a] = 1;
	    tmax[a] = (_gmin[a] + (ijk[a]+1)*_gh[a] - x[a])/d[a];
	    tdelta[a] = _gh[a]/d[a];
	} else if( d[a] < 0.0 ) {
	    step[a] = -1;
	    tmax[a] = (_gmin[a] + ijk[a]*_gh[a] - x[a])/d[a];
	    tdelta[a] = -_gh[a]/d[a];
	} else {
	    step[a] = 0;
	    tmax[a] = std::numeric_limits<double>::infinity();
	    tdelta[a] = std::numeric_limits<double>::infinity();
	}
    }

    while( true ) {
	const size_t c = (size_t)((ijk[2]*_gn[1]+ijk[1])*_gn[0]+ijk[0]);
	for( uint32_t q = _gstart[c]; q < _gstart[c+1]; q++ )
	    cand.push_back( _gtri[q] );

	int a = 0;
	if( tmax[1] < tmax[0] ) a = 1;
	if( tmax[2] < tmax[a] ) a = 2;
	if( step[a] == 0 )
	    break;
	ijk[a] += step[a];
	if( ijk[a] < 0 || ijk[a] > _gn[a]-1 )
	    break;
	tmax[a] += tdelta[a];
    }

    if( cand.empty() )
	return( false );
    std::sort( cand.begin(), cand.end() );
    cand.erase( std::unique( cand.begin(), cand.end() ), cand.end() );

    /* De-duplication of the degenerate vertex case.
     *
     * The exhaustive form allocates two std::vector<bool> the size of the
     * whole VERTEX list on every call and clears them, which costs the
     * same whether one triangle is examined or ten thousand. Here only the
     * candidates can contribute, so the same bookkeeping is done with two
     * small lists scanned linearly -- typically a handful of entries, and
     * touched at all only in the degenerate branch. Same rule, same
     * result: for a vertex shared by several triangles the first one seen
     * in triangle order claims it. */
    uint32_t vpos[8], vneg[8];
    size_t nvpos = 0, nvneg = 0;
    bool vpos_of = false, vneg_of = false;
    std::vector<uint32_t> vpos_x, vneg_x;

    int incl = 0;
    for( size_t q = 0; q < cand.size(); q++ ) {
	const uint32_t a = cand[q];
	int ss = signvol3( _vertex[_triangle[a][0]],
			   _vertex[_triangle[a][1]],
			   _vertex[_triangle[a][2]] );
	int stat = classify_original_tetrahedron( ss, x,
						  _vertex[_triangle[a][0]],
						  _vertex[_triangle[a][1]],
						  _vertex[_triangle[a][2]] );
	if( stat == VTRI_INSIDE ) {
	    incl += 2*ss;
	} else if( stat == VTRI_FACE ) {
	    incl += ss;
	} else if( stat != VTRI_OUTSIDE ) {
	    const uint32_t v = _triangle[a][stat];
	    if( ss > 0 ) {
		bool seen = false;
		for( size_t z = 0; z < nvpos && !seen; z++ )
		    if( vpos[z] == v ) seen = true;
		if( vpos_of )
		    for( size_t z = 0; z < vpos_x.size() && !seen; z++ )
			if( vpos_x[z] == v ) seen = true;
		if( !seen ) {
		    if( nvpos < 8 ) vpos[nvpos++] = v;
		    else { vpos_x.push_back( v ); vpos_of = true; }
		    incl += 2*ss;
		}
	    } else if( ss < 0 ) {
		bool seen = false;
		for( size_t z = 0; z < nvneg && !seen; z++ )
		    if( vneg[z] == v ) seen = true;
		if( vneg_of )
		    for( size_t z = 0; z < vneg_x.size() && !seen; z++ )
			if( vneg_x[z] == v ) seen = true;
		if( !seen ) {
		    if( nvneg < 8 ) vneg[nvneg++] = v;
		    else { vneg_x.push_back( v ); vneg_of = true; }
		    incl += 2*ss;
		}
	    }
	}
    }

    if( incl > 0 )
	return( true );
    return( false );
}


void VTriangleSurfaceSolid::prepare_for_inside()
{
    // Calculate bbox
    Vec3D min = Vec3D( std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity() );
    Vec3D max = Vec3D( -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _vertex.size(); a++ )
	update_bbox( min, max, _vertex[a] );

    _bbox[0] = min;
    _bbox[1] = max;

    // Set offset
    _offset = -_bbox[0] + 1.0e-3*(_bbox[1]-_bbox[0]);

    // Apply offset to vertex data
    for( size_t a = 0; a < _vertex.size(); a++ )
	_vertex[a] += _offset;

    build_grid();
}


void VTriangleSurfaceSolid::set_signed_volume_eps( double signed_volume_eps )
{
    _signed_volume_eps = signed_volume_eps;
}


void VTriangleSurfaceSolid::get_bbox( Vec3D &min, Vec3D &max ) const
{
    min = _bbox[0];
    max = _bbox[1];
}


void VTriangleSurfaceSolid::remove_duplicate_triangles( void )
{
    int removed = 0;

    /* Check all pairs
     */
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {
	const VTriangle &tria = _triangle[a];
	int a1 = tria[0];
	int a2 = tria[1];
	int a3 = tria[2];
	for( uint32_t b = a+1; b < _triangle.size(); b++ ) {
	    const VTriangle &trib = _triangle[b];
	    int b1 = trib[0];
	    int b2 = trib[1];
	    int b3 = trib[2];
	    if( a1 == b1 && a2 == b2 && a3 == b3 ) {
		// Duplicate found, remove triangle b
		removed++;
		_triangle.erase(_triangle.begin()+b);
	    }
	}
    }

    if( removed ) {
	ibsimu.message(1) << "Removed " << removed << " duplicate triangles\n";
    }
}


void VTriangleSurfaceSolid::check_data( void ) const
{
    for( uint32_t a = 0; a < _triangle.size(); a++ ) {

	const VTriangle &tri = _triangle[a];
	Vec3D V1 = _vertex[tri[1]] - _vertex[tri[0]];
	Vec3D V2 = _vertex[tri[2]] - _vertex[tri[0]];
	double area = 0.5*norm2( cross( V1, V2 ) );
	if( fabs(area) == 0.0 )
	    throw( Error( ERROR_LOCATION, "Zero area triangle " + to_string(a) ) );
	for( uint32_t b = 0; b < 3; b++ ) {
	    
	    int e1 = tri[b];
	    int e2 = tri[(b+1)%3];

	    bool found = false;
	    uint32_t neighbour;
	    uint32_t c;
	    for( c = 0; c < _triangle.size(); c++ ) {

		if( c == a ) continue; // No self-checking
		const VTriangle &tric = _triangle[c];
		for( uint32_t d = 0; d < 3; d++ ) {
		    
		    int f1 = tric[d];
		    int f2 = tric[(d+1)%3];
		    if( f1 == e1 && f2 == e2 ) {
			// Incorrect orientation
			throw( Error( ERROR_LOCATION, "Incorrect orientation between neighbouring triangles " +
				      to_string(a) + " and " + to_string(c) ) );
		    } else if( f1 == e2 && f2 == e1 ) {
			// Correct orientation
			if( found == true ) {
			    const VTriangle &trin = _triangle[neighbour];
			    throw( Error( ERROR_LOCATION, "Double neighbours (" + to_string(neighbour) +
					  " and " + to_string(c) + ") found for triangle " +
					  to_string(a) + ".\n" + 
					  "Triangle " + to_string(a) + ": " + to_string(tri[0]) + ", " + to_string(tri[1]) + ", " + to_string(tri[2]) + ":\n" +
					  "  " + to_string(_vertex[tri[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tri[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tri[2]]-_offset) + "\n" +
					  "Triangle " + to_string(neighbour) + ": " + to_string(trin[0]) + ", " + to_string(trin[1]) + ", " + to_string(trin[2]) + ":\n" +
					  "  " + to_string(_vertex[trin[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[trin[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[trin[2]]-_offset) + "\n" +
					  "Triangle " + to_string(c) + ": " + to_string(tric[0]) + ", " + to_string(tric[1]) + ", " + to_string(tric[2]) + ":\n" +
					  "  " + to_string(_vertex[tric[0]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tric[1]]-_offset) + "\n" +
					  "  " + to_string(_vertex[tric[2]]-_offset) ) );
			}
			found = true;
			neighbour = c;
		    }
		}
	    }
	    if( !found ) {
		throw( Error( ERROR_LOCATION, "Triangle " + to_string(a) + " neighbour not found."
			      "Triangle " + to_string(a) + ": " + to_string(tri[0]) + ", " + to_string(tri[1]) + ", " + to_string(tri[2]) + ":\n" +
			      "  " + to_string(_vertex[tri[0]]-_offset) + "\n" +
			      "  " + to_string(_vertex[tri[1]]-_offset) + "\n" +
			      "  " + to_string(_vertex[tri[2]]-_offset) + "\n" ) );
	    }
	}
    }
}


void VTriangleSurfaceSolid::clear( void )
{
    VTriangleSurface::clear();
    //_vpos.clear();
    //_vneg.clear();
}


void VTriangleSurfaceSolid::debug_print( std::ostream &os ) const
{
    VTriangleSurface::debug_print( os );

    os << "**VTriangleSurfaceSolid\n";

    os << "signed_volume_eps = " << _signed_volume_eps << "\n";

    os << "offset = " << _offset << "\n";
    os << "bbox[0] = " << _bbox[0] << "\n";
    os << "bbox[1] = " << _bbox[1] << "\n";
}

