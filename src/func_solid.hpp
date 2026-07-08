/*! \file func_solid.hpp
 *  \brief %Solid definition based on C functions.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef FUNC_SOLID_HPP
#define FUNC_SOLID_HPP 1


#include <iostream>
#include "solid.hpp"


/*! \brief Function solid class.
 *
 *  FuncSolid class holds the definition for one solid defining
 *  C-function. This solid implementation suffers from the inability
 *  of saving to file. If FuncSolid is constructed from stream the
 *  function pointer inside is set to NULL and error is thrown if
 *  function is evaluated using inside().
 */
class FuncSolid : public Solid {

    bool (*_func)(double,double,double);

public:

    /*! \brief Constructor.
     */
    FuncSolid( bool (*func)(double,double,double) ) : _func(func) {}

    /*! \brief Constructor for loading solid data from stream \a is.
     */
    FuncSolid( std::istream &is );

    /*! \brief Destructor.
     */
    ~FuncSolid() {}

    /*! \brief Return if point x is inside funcsolid.
     */
    bool inside( const Vec3D &x ) const;

    /*! \brief Print debugging information to \a os.
     */
    void debug_print( std::ostream &os ) const;

    /*! \brief Saves solid data to stream \a os
     */
    void save( std::ostream &os ) const;
};


#endif
