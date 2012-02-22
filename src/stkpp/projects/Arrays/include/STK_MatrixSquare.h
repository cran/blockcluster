/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  stkpp::Arrays
 * Purpose:  Define the MatrixSquare class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_MatrixSquare.h
  * @brief In this file, we define MatrixSquare class.
 **/

#ifndef STK_MATRIXSQUARE_H
#define STK_MATRIXSQUARE_H

#include "STK_Matrix.h"

namespace STK
{

class MatrixSquare;

/** @ingroup Arrays
 *  Specialization of the traits class ContainerTraits for Matrix
 *  class.
 */
template<>
struct ContainerTraits<Real, MatrixSquare >
{
   typedef Point TContainerHo;
   typedef Vector TContainerVe;
};


/** @ingroup Arrays
  * @brief Derivation of the MatrixSquare class for square arrays of
  * Real.
  *
  * A MatrixSquare is a column oriented two dimensional
  * container of Real with the same number of rows and columns.
  *
  * The range of the rows and the columns is the same.
  **/
class MatrixSquare : public Matrix
{
  public:
    /** Default constructor with rangeHo_=(1:0) and rangeVe_=(1:0).
     *  @param I range of the Rows and Cols
     **/
    MatrixSquare( Range const& I = Range());

    /** constructor with rangeHo_and rangeVe_ givens,
     *  initialization with a constant.
     *  @param I range of the Rows and Cols
     *  @param v initial value of the container
     **/
    MatrixSquare( Range const& I, Real const& v);

    /** Copy constructor.
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    MatrixSquare( MatrixSquare const&T, bool ref=false);

    /** constructor by reference, ref_=1 in the range given by I.
     *  @param T the Container to wrap
     *  @param I range of the container to wrap
     **/
    MatrixSquare( Matrix const& T, Range const& I);

    /** destructor.                                                         */
    virtual ~MatrixSquare();

    /** get range of the Matrix.
     * @return the range of the matrix
     **/
    inline Range const& range() const
    { return rangeVe();}

    /** get the first index of the Matrix.
     * @return the first index of the matrix
     **/
    inline Integer const& first() const
    { return firstCol();}

    /** get the last index of the Matrix.
     * @return the last index of the matrix
     **/
    inline Integer const& last() const
    { return lastCol();}

    /** access to a square sub-array.
     *  @param I range of the data to wrap
     **/
    inline MatrixSquare sub(Range const& I) const
    { return MatrixSquare(*this, I);}

    /** New beginning index for the object.
     *  @param beg first index of the container
     **/
    inline void shift(Integer beg)
    { Matrix::shift(beg, beg);}

    /** New size for the container.
     *  @param I range of the columns and rows of the container
     **/
    inline void resize( Range const& I =Range())
    { Matrix::resize(I, I);}

    /** Insert n rows and column at the given position to the container.
     *  @param pos position to insert the Rows and Cols
     *  @param n number of Rows and Cols insert
     **/
    inline void insert( Integer const& pos, Integer const& n =1)
    {
      Matrix::insertRows(pos, n);
      Matrix::insertCols(pos, n);
    }

    /** Delete n rows and columns at the specified position to
     *  the container.
     *  @param pos position to erase the Rows and Cols
     *  @param n number of Rows and Cols erase
     **/
    inline void erase( Integer const& pos, Integer const& n=1)
    {
      Matrix::eraseCols(pos, n);
      Matrix::eraseRows(pos, n);
    }

    /** Add n rows and columns to the container.
     *  @param n number of Rows and Cols to add
     **/
    inline void pushBack(Integer n=1)
    {
      Matrix::pushBackRows(n);
      Matrix::pushBackCols(n);
    }

    /** Delete n rows and columns at the end of the container.
     *  @param n number of Rows and Cols to delete
     **/
    inline void popBack(Integer const& n=1)
    {
      Matrix::popBackCols(n);
      Matrix::popBackRows(n);
    }

    /** operator = : overwrite the MatrixSquare with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline MatrixSquare& operator=( MatrixSquare const& T)
    {
            Matrix *p1(this); // convert this to Array2D
      const Matrix *p2(&T);   // convert T    to Array2D
      (*p1) = (*p2);                   // use = of Array2D class

      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline MatrixSquare& operator=(Real const& v)
    {
      Matrix *p1(this); // convert this to Array2D
      (*p1) = v;        // use = of Array2D class

      return *this;
    }

    /** operator = overwriting a MatrixSquare with a container.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    template < class TContainer2D>
    inline MatrixSquare& operator=( ITArray2D< Real, TContainer2D> const& T)
    {
      if (T.rangeVe() != T.rangeHo())
      { throw runtime_error("MatrixSquare::operator=(T) "
                                 "dimensions of T are not square.");
      }
      // use Matrix = operator
      Matrix *p1(this); // convert this to Matrix
      (*p1) = T;        // use = of Matrix class
      return *this;
    }

  private:
    /** Resize the container. The range of the rows and of the
     * columns have to be the same.
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    inline void resize( Range const& I, Range const& J)
    {
      if (I != J)
      { throw runtime_error("MatrixSquare::resize(I, J) "
                                 "I != J");
      }
      Matrix::resize(I, J);
    }

    /** overloading of the shift method. The first index of the row
     * must be the first index of the columns.
     *  @param rbeg first index of the rows
     *  @param cbeg first index of the columns
     **/
    inline void shift(Integer rbeg, Integer cbeg)
    {
      if (rbeg != cbeg)
      { throw runtime_error("MatrixSquare::shift(rbeg, cbeg) "
                                 "rbeg != cbeg");
      }
      Matrix::shift(rbeg, cbeg);
    }



};

} // namespace STK


#endif
// STK_MATRIXSQUARE_H
