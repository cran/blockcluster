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
 * Purpose:  Define the MatrixLowerTriangular class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_MatrixLowerTriangular.h
  * @brief In this file we define the MatrixLowerTriangular class
 **/

#ifndef STK_MATRIXLOWERTRIANGULAR_H
#define STK_MATRIXLOWERTRIANGULAR_H

#include "STK_Vector.h"
#include "STK_Point.h"

#include "STK_ITArray2D.h"
#include "STK_Display2D.h"

#include "../../Sdk/include/STK_ContainerTraits.h"

namespace STK
{
class MatrixLowerTriangular;
/** @ingroup Arrays
 *  Specialization of the traits class ContainerTraits for MatrixLowerTriangular
 *  class.
 */
template<>
struct ContainerTraits<Real, MatrixLowerTriangular>
{
   typedef Point TContainerHo;
   typedef Vector TContainerVe;
};


/** @ingroup Arrays
  * @brief Declaration of the lower triangular matrix class.
  *
  * A MatrixLowerTriangular is a column oriented 2D lower triangular
  * container of Real. It is possible to add/remove rows and columns
  * but in this case the container will no more be triangular.
  * The container can be set lower triangular again using the method
  * ITArray2D::update().
 **/
class MatrixLowerTriangular : public ITArray2D< Real, MatrixLowerTriangular>
{
  private:
    /** Default constant value (Real(0)).**/
    const Real defaultConst_;

    /** Default non constant value.**/
    Real default_;

  public:
    /** Real for the Interface Class.                                 */
    typedef ITArray2D<Real, MatrixLowerTriangular> _ITArray2DType;

    /** Default constructor Default is I=(1:0) and J=(1:0)
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    MatrixLowerTriangular( Range const& I = Range(), Range const& J = Range());

    /** constructor with rangeVe_ and rageHo_ specified, initialization with a
     *  specified value.
     *  @param I range of the Rows
     *  @param J range of the Cols
     *  @param v initial value in the container
     **/
    MatrixLowerTriangular( Range const& I, Range const& J, Real const& v);

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    MatrixLowerTriangular( const MatrixLowerTriangular &T, bool ref=false);

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    MatrixLowerTriangular( const _ITArray2DType& T, Range const& I, Range const& J);

    /** Wrapper constructor Contruct a reference container.
     *  @param q pointer on the data
     *  @param I range of the  Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    MatrixLowerTriangular( Real** q, Range const& I, Range const& J);

    /** virtual destructor : use destructor of Array2D.                           */
    virtual ~MatrixLowerTriangular();

    /** Compute the first index of the column icol.
     *  @param icol the index of the column we want to compute the
     *  first index
     **/
    inline Integer compFirstVe( Integer const& icol) const
    {
      return max( this->firstRow()+ icol - this->firstCol(), this->firstRow());
    }

    /** Compute the last index of the column icol.
     *  For a lower triangular matrix, this is the index of the last row.
     *  @param icol the column we want to know the last index
     **/
    inline Integer compLastVe( Integer const& icol) const
    {
      return this->lastRow();
    }

    /** Compute the size of the column icol.
     *  @param icol the column we want to know the size
     **/
    inline Integer compSizeVe( Integer const& icol) const
    {
      return max(this->lastRow()-compFirstVe(icol)+1, Integer(0));
    }

    /** compute the range of the effectively stored elements in the col
     *  icol.
     *  @param icol the number of the column to compute the range
     **/
    inline Range compRangeVe( Integer const& icol) const
    {
      return Range(compFirstVe(icol), this->lastRow());
    }

    /** compute the beginning of the row irow.
     *  @param irow the column to compute the beginning
     **/
    inline Integer compFirstHo( Integer const& irow) const
    {
      return this->firstCol();
    }

    /** compute the end of the row irow.
     *  @param irow the column to compute the end
     **/
    inline Integer compLastHo( Integer const& irow) const
    {
      return min( this->firstCol()+ irow - this->firstRow(),
                 this->lastCol());
    }

    /** compute the true number of element of the row irow.
     *  @param irow the number of the column to compute the size
     **/
    inline Integer compSizeHo( Integer const& irow) const
    {
      return max(compLastHo(irow)-this->firstCol()+1, Integer(0));
    }

    /** compute the range of the effectively stored elements in the row
     *  irow.
     *  @param irow the number of the row to compute the range
     **/
    inline Range compRangeHo( Integer const& irow) const
    {
      return Range(this->firstCol(), this->compLastHo(irow));
    }

    /** function for determining if the row i of the col
     *  j is in the lower triangular part.
     * @param i the number of the row
     * @param j the number of the col
     **/
    inline bool isInside(Integer const& i, Integer const& j) const
    {
      return (i>=compFirstVe(j));
    }

    /** access to one element.
     *  @param i index of the row
     *  @param j index of the col
     **/
    inline Real& elt(Integer i, Integer j)
    {
      return isInside(i, j) ? this->data(j)[i] : default_;
    }

    /** access to one element const.
     *  @param i index of the row
     *  @param j index of the col
     **/
    inline Real const& elt(Integer i, Integer j) const
    {
      // we dont use ? : operator in order to avoid a message from the compiler:
      // return reference to temporary (???)
      if (isInside(i, j))
      {
        return this->data(j)[i];
      }
      else
      {
        return defaultConst_;
      }
    }

    /** access to a sub-array.
     *  @param I range of the rows
     *  @param J range of the cols
     **/
    inline MatrixLowerTriangular sub(Range const& I, Range const& J) const
    {
      return MatrixLowerTriangular(*this, I, J);
    }

    /** access to a part of a column.
     *  @param I range of the rows
     *  @param j index of the col
     *  @return a reference in the range I of the column j of this
     **/
    inline Vector col(Range const& I, Integer j) const
    {
      return Vector( this->data(j), Range::inf(I, compRangeVe(j)), j);
    }

    /** access to a part of a column.
     *  @param j index of the column
     *  @return a reference in the range I of the column j of this
     **/
    inline Vector col( Integer j) const
    {
      return Vector( this->data(j), compRangeVe(j), j);
    }

    /** access to a part of a row.
     *  @param i index of the row
     *  @param J range of the columns
     *  @return a reference of the row i.
     **/
    inline Point row(Integer i, Range const& J) const
    {
      return Point(*this, Range::inf(J, compRangeHo(i)), i);
    }

    /** access to a part of a row.
     *  @param i index of the row
     *  @return a reference of the row i.
     **/
    inline Point row( Integer i) const
    {
      return Point(*this, compRangeHo(i), i);
    }

    /** operator = : overwrite the MatrixLowerTriangular with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    MatrixLowerTriangular& operator=( MatrixLowerTriangular const& T);

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline MatrixLowerTriangular& operator=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      {
        Real* p(this->data(j));
        Integer beg(this->rangeCols_[j].first());
        Integer end(this->rangeCols_[j].last());

        for (Integer i=beg; i<=end; i++) p[i]= v;
      }
      return *this;
    }
};

/** @ingroup Arrays
 *  ostream for MatrixLowerTriangular.
 *  @param s the output stream
 *  @param V the MatrixLowerTriangular to write
 **/
ostream& operator<<(ostream& s, MatrixLowerTriangular const& V);

} // namespace STK

#endif
// STK_MATRIXLOWERTRIANGULAR_H
