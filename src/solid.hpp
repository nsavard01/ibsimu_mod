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

#ifndef SOLID_HPP
#define SOLID_HPP 1


#include <iostream>
#include "vec3d.hpp"
#include "transformation.hpp"


/*! \brief %Solid base class.
 *
 *  %Solid class holds the definition for one solid. %Solid class is
 *  a base class, different implementation exist.
 */
class Solid {

protected:

    Transformation         _T;

    /*! \brief Constructor.
     */
    Solid();

    /*! \brief Sentinel magnitude used by bbox_from_local_box() (and by
     *  get_bbox() implementations) to flag a world-space axis along
     *  which no finite bound could be established. Callers of
     *  get_bbox() must treat any returned axis whose magnitude is
     *  \f$ \ge 10^{10} \f$ as unrestricted -- not as a literal
     *  coordinate limit.
     */
    static const double BBOX_HUGE;

    /*! \brief Turn a solid-local-space (pre-\a _T) axis-aligned box
     *  into a conservative world-space axis-aligned box.
     *
     *  \a lo and \a hi are corners of a box in the solid's native
     *  coordinate frame, i.e. the frame in which inside() tests
     *  actually happen after \a _T has been applied to the query
     *  point (world = _T.inverse() * local). Components of \a lo/\a
     *  hi may be +-BBOX_HUGE to mean "unbounded along this local
     *  axis" (do not pass +-infinity: multiplying an infinite
     *  coordinate by an exactly-zero matrix entry, which happens
     *  routinely for axis-aligned rotations, produces NaN; a large
     *  finite sentinel does not).
     *
     *  All 8 corners of the local box are mapped to world space with
     *  _T.inverse() and the axis-aligned bounding box of the result
     *  is returned in \a wmin/\a wmax. Because arbitrary rotations in
     *  \a _T can mix a genuinely unbounded local axis into any world
     *  axis, each world axis found to exceed a sanity threshold is
     *  reset to +-BBOX_HUGE to flag it as unrestricted rather than
     *  returning a meaningless huge number as if it were real.
     *
     *  Returns \a false if every world axis ended up unrestricted
     *  (i.e. the box is useless and the caller should fall back to
     *  scanning the whole mesh), \a true if at least one axis carries
     *  a usable finite bound.
     */
    bool bbox_from_local_box( const Vec3D &lo, const Vec3D &hi, Vec3D &wmin, Vec3D &wmax ) const;

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Solid();

    /*! \brief Return if point x is inside solid.
     */
    virtual bool inside( const Vec3D &x ) const = 0;

    /*! \brief Return a conservative world-space axis-aligned bounding
     *  box for the solid in \a min and \a max.
     *
     *  Returns \a false if no useful bound is known (the default),
     *  in which case \a min/\a max are left untouched and the solid
     *  must be assumed to possibly occupy the entire simulation
     *  domain. A \a true return guarantees every point actually
     *  inside the solid satisfies \a min <= x <= \a max
     *  component-wise, but does not guarantee the reverse -- points
     *  inside the box need not be inside the solid. Some axes of a
     *  box returned as \a true may still be individually unrestricted;
     *  see bbox_from_local_box().
     */
    virtual bool get_bbox( Vec3D &min, Vec3D &max ) const {
	(void)min; (void)max;
	return( false );
    }

    /*! \brief Set transformation to unity.
     *
     *  Resets the primary 3D to 3D transformation to unity.
     */
    void reset_transformation( void );

    /*! \brief Set transformation.
     *
     *  Sets the primary 3D to 3D transformation as a copy of
     *  transformation of \a T.
     */
    void set_transformation( const Transformation &T );

    /*! \brief Translate solid.
     */
    void translate( const Vec3D &dx );

    /*! \brief Scale solid.
     */
    void scale( double sx );

    /*! \brief Scale solid.
     */
    void scale( const Vec3D &sx );

    /*! \brief Rotate solid around x-axis.
     *
     *  Rotate around x-axis for \a a radians.
     */
    void rotate_x( double a );

    /*! \brief Rotate solid around y-axis.
     *
     *  Rotate around y-axis for \a a radians.
     */
    void rotate_y( double a );

    /*! \brief Rotate solid around z-axis.
     *
     *  Rotate around z-axis for \a a radians.
     */
    void rotate_z( double a );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const = 0;

    /*! \brief Saves solid data to stream.
     */
    virtual void save( std::ostream &s ) const = 0;
};


#endif
