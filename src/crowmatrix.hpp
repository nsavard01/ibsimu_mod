/*! \file crowmatrix.hpp
 *  \brief Compressed row sparse matrix algebra
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

#ifndef CROWMATRIX_HPP
#define CROWMATRIX_HPP 1

#include "config.h"
#ifdef USE_MKL
#include <mkl.h>
#endif
#include <cstdlib>
#include <iostream>

#include "matrix.hpp"
#include "error.hpp"


/*! \brief Compressed row sparse matrix class.
 *
 *  The matrix is stored in the standard compressed row sparse matrix
 *  storage mode.  In compressed row storage method all non-zero
 *  matrix elements are stored in array \a val in row-by-row
 *  order. The corresponding column indices are stored in another
 *  array \a col in the same order. The third array \a ptr contains
 *  "pointer" indices indicating start and end of each row in the
 *  first two arrays. The format itself does not require a certain
 *  ordering of elements, but some implementation might need/be faster
 *  on some ordering. The elements can be ordered into ascending index 
 *  order with order_ascending(). Our example matrix \code
 *      | 1  2  0  0  3|
 *      | 4  5  6  0  0|
 *  A = | 0  7  8  0  9|
 *      | 0  0  0 10  0|
 *      |11  0  0  0 12|
 *  \endcode is represented in compressed row sparse matrix class as: \code
 *  ptr[] = {0, 3, 6, 9, 10, 12}
 *  val[] = {1, 2, 3, 4,  5,  6, 7, 8, 9, 10, 11, 12}
 *  col[] = {0, 1, 4, 0,  1,  2, 1, 2, 4,  3,  0,  4}
 *  \endcode
 */
class CRowMatrix : public Matrix {
    int       _n;      //!< Number of rows.
    int       _m;      //!< Number of columns.
    int       _nz;     //!< Number of nonzero elements.
    int       _asize;  //!< Allocation size of \a _col and \a _val.
    int      *_ptr;    //!< Row pointers, size \a _n+1.
    int      *_col;    //!< Column indices(j), \a _nz elements in use, \a _asize elements allocated.
    double   *_val;    //!< Element values, \a _nz elements in use, \a _asize elements allocated.

#ifdef USE_MKL
    /* 1-based copies of _ptr and _col for mkl_dcsrgemv/mkl_dcsrtrsv.
     * Plain malloc/free pointers, same style as _ptr/_col/_val above.
     * Rebuilt lazily when _mkl_dirty is true (any structural change). */
    
    mutable sparse_matrix_t _mklA;

    mutable bool _mkl_valid;
    mutable bool _mkl_structure_dirty;
    mutable bool _mkl_values_dirty;


    mutable bool _mkl_dirty;

    void mkl_build_handle() const;
#endif

    /*! \brief Mark any cached MKL sparse handle stale.
     *
     *  Must be called from EVERY mutation of _ptr/_col/_val -- values,
     *  structure and reallocation alike -- not just the structural ones.
     *  The handle created by mkl_sparse_d_create_csr() aliases these
     *  arrays, and mkl_sparse_optimize() is additionally permitted to
     *  keep its own reordered copy, so:
     *
     *    - a realloc (reserve()/reallocate()) or free (resize()/clear())
     *      leaves the handle holding dangling pointers -- a use-after-free
     *      on the next mkl_sparse_d_mv();
     *    - an in-place value overwrite (set_no_check()'s existing-element
     *      path, used by EpotMatrixSolver::get_resjac() on every Newton
     *      iteration to refresh the Jacobian diagonal) may be invisible to
     *      an optimized internal copy, silently multiplying by a stale
     *      matrix;
     *    - a permutation (order_ascending()) changes which value belongs
     *      to which column without changing any pointer at all.
     *
     *  Whether the last two actually bite depends on MKL's internal
     *  choices for a given version/matrix, which is exactly why this is
     *  handled defensively rather than relying on today's behaviour.
     *  Compiles to nothing without USE_MKL.
     */
    inline void mkl_invalidate( void ) const
    {
#ifdef USE_MKL
	_mkl_dirty = true;
#endif
    }

    void reallocate( void );
    void allocate( void );

    double get_check( int i, int j ) const;
    double &set_check( int i, int j );
    double get_no_check( int i, int j ) const;
    double &set_no_check( int i, int j );

    void clear_check( int i, int j );
    void clear_no_check( int i, int j );

    void build( const class CColMatrix &mat );
    void build( const class CRowMatrix &mat );
    void build( const class CoordMatrix &mat );

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Default constructor.
     */
    CRowMatrix();

    /*! \brief Constructor to make empty \a n x \a m matrix.
     */
    CRowMatrix( int n, int m );

    /*! \brief Constructor to make \a n x \a m matrix from compressed
     *  row matrix data.
     *
     *  The constructed matrix uses \a ptr, \a col and \a val as its
     *  internal data. The arrays are not copied! For memory
     *  allocation compatibility reasons, arrays \a ptr, \a col and \a
     *  val should be allocated using \a malloc and/or \a realloc.  \a
     *  nz is the number of non-zeros and \a asize is the number of
     *  allocated elements.
     */
    CRowMatrix( int n, int m, int nz, int asize,
		int *ptr, int *col, double *val );

    /*! \brief Copy constructor.
     */
    CRowMatrix( const CRowMatrix &mat );

    /*! \brief Constructor for conversion from compressed column matrix.
     *
     *  The compressed row matrix will be built in ascending column
     *  index order.
     */
    CRowMatrix( const class CColMatrix &mat );

    /*! \brief Constructor for conversion from coordinate matrix.
     */
    CRowMatrix( const class CoordMatrix &mat );

    /*! \brief Constructor for conversion from unknown matrix type.
     */
    CRowMatrix( const class Matrix &mat );

    /*! \brief Destructor.
     */
    ~CRowMatrix();

/* ************************************** *
 * Access and information                 *
 * ************************************** */

    /*! \brief Returns the number of columns in the matrix.
     */
    int columns( void ) const { return( _m ); }

    /*! \brief Returns the number of rows in the matrix.
     */
    int rows( void ) const { return( _n ); }

    /*! \brief Returns the number of columns and number of columns in
     *  \a nn and \a mm.
     */
    void size( int &n, int &m ) const { n = _n; m = _m; }

    /*! \brief Returns the number of non-zero elements in the matrix.
     */
    int nz_elements( void ) const { return( _nz ); }

    /*! \brief Returns the number of elements allocated for matrix.
     */
    int capacity( void ) const { return( _asize ); }

/* ************************************** *
 * User level control                     *
 * ************************************** */

    /*! \brief Resizes the matrix to size \a n x \a m.
     *
     *  All existing non-zero elements are cleared.
     */
    void resize( int n, int m );

    /*! \brief Merges matrix \a mat into the matrix leaving \a mat empty.
     *
     *  Copies contents of matrix \a mat into the matrix and sets
     *  contents of matrix \a mat to \a n = 0 and \a m = 0.
     *  \param mat Matrix to copy from.
     */
    void merge( CRowMatrix &mat );

    /*! \brief Clear non-zero matrix elements, set all elements to zero.
     */
    void clear( void );

    /*! \brief Clear matrix element (i,j).
     *
     *  Removes element (i,j) from the list of non-zero matrix elements.
     */
    void clear( int i, int j );

    /*! \brief Reserve memory for \a size matrix elements.
     */
    void reserve( int size );

    /*! \brief Order (sort) matrix data in ascending column index
     *  order within each row.
     */
    void order_ascending( void );

    /*! \brief Check if matrix data is in ascending column index
     *  order within each row.
     */
    bool check_ascending( void ) const;

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;

/* ************************************** *
 * User level matrix element access       *
 * ************************************** */

    /*! \brief Function to get a matrix element value at (i,j).
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double get( int i, int j ) const;

    /*! \brief Function to get a reference to matrix element value at (i,j).
     *
     *  This function can be used to set or modify matrix element
     *  value. See following examples: \code
     *  A.set(0,0) = 1.2
     *  A.set(0,1) *= 2
     *  A.set(0,1) += 0.1 
     *  \endcode Note that a
     *  reference is actually a pointer to the memory location of the
     *  element and therefore unexpected things can happen if matrix
     *  is modified while using set, for example \code 
     *  A.set(0,0) = A.set(0,1) = A.set(0,2) = 5.0 
     *  \endcode does not do what you
     *  would expect it to do.
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double &set( int i, int j );

    /*! \brief Return the flat index into the value array of the element
     *  at (\a i, \a j), or -1 if that element is not stored.
     *
     *  Intended to be resolved ONCE for a set of elements whose
     *  positions are then reused many times -- the canonical case being
     *  EpotMatrixSolver's Jacobian diagonal, rewritten on every Newton
     *  iteration. set() has to linear-scan row \a i on every single call
     *  to find the element, so refreshing d diagonal entries costs
     *  O(nnz) scanning per sweep; caching these indices once turns each
     *  later sweep into d direct array writes.
     *
     *  Only valid while the matrix STRUCTURE is unchanged. Any call that
     *  inserts/removes an element or reorders a row (set() on a
     *  not-yet-present element, clear(), order_ascending(), resize(),
     *  reserve(), set_row(), construct_add(), ...) invalidates every
     *  previously returned index. Values may be freely overwritten
     *  through value_ptr() in the meantime.
     */
    int value_index( int i, int j ) const;

    /*! \brief Raw, writable pointer to the value array, for use with
     *  value_index().
     *
     *  Callers writing through this pointer MUST call values_changed()
     *  afterwards, since bypassing set() also bypasses its automatic
     *  invalidation of any cached MKL sparse handle. Doing it this way
     *  round is deliberate: it lets a bulk update mark the matrix stale
     *  exactly once instead of once per element.
     */
    double *value_ptr( void ) { return( _val ); }

    /*! \brief Const overload of value_ptr().
     */
    const double *value_ptr( void ) const { return( _val ); }

    /*! \brief Declare that element values were modified directly through
     *  value_ptr(), so any cached MKL sparse handle must be rebuilt.
     *
     *  Not needed after set()/construct_add()/set_row(), which already
     *  do this themselves.
     */
    void values_changed( void ) { mkl_invalidate(); }

    /*! \brief Function to set matrix row elements.
     *
     *  The existing elements of the row \a i are replaced by \a N
     *  elements at column indices described in array \a col and with
     *  values described in array \a val. Range checking is always
     *  done for indexes \a i and \a j. Throws ErrorRange exception on
     *  range checking errors.
     */
    void set_row( int i, int N, const int *col, const double *val );

    /*! \brief Adds an element to matrix while constructing the whole
        matrix.
     *
     *  This is a special function for constructing the whole matrix
     *  row-by-row in ascending row order. The elements within a row
     *  can be defined in any order. With this function, every row of
     *  the matrix has to be defined. No checking of the definitions
     *  are made.
     *
     *  Using this function for defining a large matrix gains
     *  drasticly in speed. The function leaves all but the next row
     *  pointers unmodified. Therefore the matrix is invalid and
     *  should not be accessed with other functions before all rows
     *  have been defined.
     */
    void construct_add( int i, int j, double val );

    const int *ptr_data() const { return _ptr; }
    const int *col_data() const { return _col; }
    const double *val_data() const { return _val; }

/* ************************************** *
 * Low level access                       *
 * ************************************** */

    /*! \brief Returns a reference to the to the internal pointer index
     *  data \a ptr of the matrix.
     */
    int &ptr( int i ) { return( _ptr[i] ); }

    /*! \brief Returns a reference to the to the internal column data
     *  of the matrix.
     */
    int &col( int i ) { return( _col[i] ); }

    /*! \brief Returns a reference to the to the internal value data
     *  of the matrix.
     */
    double &val( int i ) { return( _val[i] ); }

    /*! \brief Returns a const reference to the to the internal
     *  pointer index data \a ptr of the matrix.
     */
    const int &ptr( int i ) const { return( _ptr[i] ); }

    /*! \brief Returns a const reference to the to the internal column
     *  data of the matrix.
     */
    const int &col( int i ) const { return( _col[i] ); }

    /*! \brief Returns a const reference to the to the internal value
     *  data of the matrix.
     */
    const double &val( int i ) const { return( _val[i] ); }

    /*! \brief Set number of non-zero elements in the matrix.
     *
     *  This function is to be used with low level access
     *  functions. The number of non-zero elements should be set to
     *  the same value as \a ptr[n]. Internal arrays are resized if \a
     *  nz is larger than the allocated size.
     */
    void set_nz( int nz );

/* ************************************** *
 * Assignent operators                    *
 * ************************************** */

    /*! \brief Assignment operator.
     */
    CRowMatrix &operator=( const CRowMatrix &mat );

    /*! \brief Assignment operator for conversion from compressed column
     *  matrix.
     *
     *  The compressed row matrix will be built in ascending column
     *  index order.
     */
    CRowMatrix &operator=( const class CColMatrix &mat );

    /*! \brief Assignment operator for conversion from coordinate matrix.
     */
    CRowMatrix &operator=( const class CoordMatrix &mat );

    /*! \brief Assignment operator for conversion from unknown matrix
     *  type.
     */
    CRowMatrix &operator=( const class Matrix &mat );

/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */

    /*  \brief Calculates \a x = \a A*b.
     */
    void multiply_by_vector( Vector &x, const Vector &b ) const;

    /*! \brief Solves \a A*x = \a b for lower unit diagonal matrix.
     *
     *  Matrix has to have elements only in the lower triangle. Unit
     *  diagonal is implied, it is not to be saved to matrix.
     */
    void lower_unit_solve( Vector &x, const Vector &b ) const;

    /*! \brief Solves \a A*x = \a b for upper diagonal matrix.
     *
     *  Matrix has elements only in the upper triangle.
     *
     *  Assumes matrix is in sorted ascending order.
     */
    void upper_diag_solve( Vector &x, const Vector &b ) const;

    /*! \brief Solves \a A*x = \a b for packed LU matrix.
     *
     *  The \a L and \a U matrices are packed in a single matrix \a A.
     *  The \a L matrix is a lower triangular matrix with unit
     *  diagonal (unit diagonal not saved in \a A) and the \a U matrix
     *  is upper diagonal matrix.
     *
     *  Assumes matrix is in sorted ascending order.
     */
    void LU_solve( Vector &x, const Vector &b ) const;

    friend class CColMatrix;
    friend class CoordMatrix;
    friend class Vector;
};


inline double CRowMatrix::get( int i, int j ) const
{
#ifdef SPM_RANGE_CHECK
    return( get_check( i, j ) );
#else
    return( get_no_check( i, j ) );
#endif
}    


inline double &CRowMatrix::set( int i, int j )
{
#ifdef SPM_RANGE_CHECK
    return( set_check( i, j ) );
#else
    return( set_no_check( i, j ) );
#endif
}    


inline void CRowMatrix::clear( int i, int j )
{
#ifdef SPM_RANGE_CHECK
    clear_check( i, j );
#else
    clear_no_check( i, j );
#endif
}


#endif