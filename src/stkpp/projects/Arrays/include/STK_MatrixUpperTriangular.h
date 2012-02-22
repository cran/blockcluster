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
 * Purpose:  Define the Upper Triangular Matrix class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_MatrixUpperTriangular.h
  * @brief In this file we define the MatrixTriangular class
 **/

#ifndef STK_MATRIXUPPERTRIANGULAR_H
#define STK_MATRIXUPPERTRIANGULAR_H

#include "STK_Vector.h"
#include "STK_Point.h"

#include "STK_ITArray2D.h"
#include "STK_Display2D.h"

#include "../../Sdk/include/STK_ContainerTraits.h"

namespace STK
{
class MatrixUpperTriangular;
/** @ingroup Arrays
 *  Specialization of the traits class ContainerTraits for MatrixLowerTriangular
 *  class.
 */
template<>
struct ContainerTraits<Real, MatrixUpperTriangular>
{
   typedef Point TContainerHo;
   typedef Vector TContainerVe;
};


/** @ingroup Arrays
  * @brief Declaration of the upper triangular matrix class.
  *
  * A MatrixUpperTriangular is a column oriented 2D upper triangular
  * container of Real. It is possible to add/remove rows and columns
  * but in this case the container will no more be triangular.
  * The container can be set upper triangular again using the method
  * ITArray2D::update().
  **/
class MatrixUpperTriangular : public ITArray2D< Real, MatrixUpperTriangular>
{
  private:
    /** Default constant value (Real(0)).**/
    const Real defaultConst_;

    /** Default non constant value.**/
    Real default_;

  public:
    /** Column Real.                                                  */
    typedef Array1D<Real> Vector;

    /** Row Real.                                                     */
    typedef ArrayHo<Real> Point;

    /** Real for the Base reference Class.                            */
    typedef AllocatorBase<Real*> _AllocatorBaseType_;

    /** Real for the Interface Class.                                 */
    typedef ITArray2D< Real, MatrixUpperTriangular >
           _ITArray2DType;

    /** Default constructor Default is I=(1:0) and J=(1:0).
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    MatrixUpperTriangular( Range const& I = Range(), Range const& J = Range());

    /** constructor with rangeVe_ and rageHo_ specified, initialization with a
     *  specified value.
     *  @param I range of the Rows
     *  @param J range of the Cols
     *  @param v initial value in the container
     **/
    MatrixUpperTriangular( Range const& I, Range const& J, Real const& v);

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    MatrixUpperTriangular( const MatrixUpperTriangular& T
                         , bool ref=false
                         );

    /** constructor by reference in a given range, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    MatrixUpperTriangular( const _ITArray2DType& T
                         , Range const& I
                         , Range const& J
                         );

    /** Wrapper constructor Contruct a reference container.
     *  @param q pointer on the data
     *  @param I range of the  Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    MatrixUpperTriangular( Real** q, Range const& I, Range const& J);

    /** virtual destructor.                                                 */
    virtual ~MatrixUpperTriangular();

    /** compute the beginning of the row irow.
     *  @param irow the row to compute the begining
     **/
    inline Integer compFirstHo( Integer const& irow) const
    {
      return max( this->firstCol() + irow - this->firstRow()
                , this->firstCol()
                );
    }

    /** compute the end of the row irow.
     *  @param irow the row to compute the beginning
     **/
    inline Integer compLastHo( Integer const& irow) const
    {
      return this->lastCol();
    }

    /** compute the true number of element of the row irow.
     *  @param irow the index of the row we want to compute the size
     **/
    inline Integer compSizeHo( Integer const& irow) const
    {
      return max(this->lastCol() - compFirstHo(irow) + 1, (Integer )0);
    }

    /** Compute the range of the effectively stored elements in the row
     *  irow.
     *  @param irow the number of the row we want to compute the range
     **/
    inline Range compRangeHo( Integer const& irow) const
    {
      return Range(compFirstHo(irow), this->lastCol());
    }

    /** compute the begin of the column icol.
     *  @param icol the number of the column we want to compute the first index
     **/
    inline Integer compFirstVe( Integer const& icol) const
    {
      return this->firstRow();
    }

    /** compute the end of the column icol.
     *  @param icol the row to compute the end
     **/
    inline Integer compLastVe( Integer const& icol) const
    {
      return min( this->firstRow() + icol - this->firstCol()
                , this->lastRow()
                );
    }

    /** compute the true number of element of the column icol.
     *  @param icol the number of the column to compute the size
     **/
    inline Integer compSizeVe( Integer const& icol) const
    {
      return max(compLastVe(icol) - this->firstRow() + 1, (Integer )0);
    }

    /** compute the range of the effectively stored elements in the col
     *  icol.
     *  @param icol the number of the column to compute the range
     **/
    inline Range compRangeVe( Integer const& icol) const
    {
      return Range(this->firstRow(), this->compLastVe(icol));
    }

    /** private function for determining if the element i of the col
     *  j is in the upper triangular part.
     **/
   inline bool isInside(Integer const& i, Integer const& j) const
   { return (i<=this->rangeCols_[j].last());}

    /** access to one element.                                        */
    inline Real& elt(Integer i, Integer j)
    { return isInside(i,j) ? this->data(j)[i] : default_;}

    /** access to one element const.                                  */
    inline Real const& elt(Integer i, Integer j) const
    {
      if (isInside(i,j)) { return this->data(j)[i];}
      else               { return defaultConst_;}
    }
    /** access to a sub-array.                                        */
    inline MatrixUpperTriangular sub(Range const& I, Range const& J) const
    { return MatrixUpperTriangular(*this, I, J);}

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

    /** Operator = : overwrite the MatrixUpperTriangular with T.      */
    MatrixUpperTriangular& operator=(const MatrixUpperTriangular &T);

    /** Operator = : overwrite with a constant value.
     **/
    inline MatrixUpperTriangular& operator=(Real const& v)
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
 *  ostream for MatrixUpperTriangular.
 *  @param s the output stream
 *  @param V the MatrixUpperTriangular to write
 **/
ostream& operator<<(ostream& s, MatrixUpperTriangular const& V);

} // namespace STK

#endif
// STK_MATRIXUPPERTRIANGULAR_H
