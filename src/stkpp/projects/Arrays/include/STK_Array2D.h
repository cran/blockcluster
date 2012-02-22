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
 * Project: stkpp::Arrays
 * created on: 07 jul. 2007
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Array2D.h
  * @brief In this file, we define the final class @c Array2D.
 **/

#ifndef STK_ARRAY2D_H
#define STK_ARRAY2D_H

#include "STK_ITArray2D.h"
#include "STK_Display2D.h"

#include "../../Sdk/include/STK_ContainerTraits.h"


namespace STK
{
template<typename> class Array2D;
/** @ingroup Arrays
 *  Specialization of the traits class ContainerTraits for Matrix
 *  class.
 */
template<class TYPE>
struct ContainerTraits<TYPE, Array2D<TYPE> >
{
   typedef ArrayHo<TYPE> TContainerHo;
   typedef Array1D<TYPE> TContainerVe;
};

/** @ingroup Arrays
  * @brief Templated two dimensional column (vertically) oriented Array.
  *
  * A Array2D is a two-dimensional implementation of an ITArray2D.
  *
  * A column of an Array2D is (almost like) a Array1D
  * and a Row of an Array2D is (almost like) a ArrayHo.
  *
  * When accessing to a row or a column, the methods return
  * the ArrayHo or Array1D as a reference-like and without copying the data.
 *
 * The elements of the Array2D can be acceded with the operators () and [].
 * - T(2,3) allow to access the third element of the second row of T
 * - T(3,Range(1,2)) allow to access to the first two members of the third
 *  row of T
 * - T(Range(2,3),2) allow to access to the second and third element of the
 *   second column of  T.
 * - T(Range(2,3),Range(1,2)) allow to access to the specified sub-matrix in T.
 * - T[2] allow to access to the second column of T.
 * - T(3) allow to access to the third row of T.
 **/
template<class TYPE >
class Array2D : public ITArray2D< TYPE, Array2D<TYPE> >
{
  /** Type for the Interface base Class. */
  typedef ITArray2D<TYPE, Array2D<TYPE> > _ITArray2DType;

  public:
//    /** Column type.  */
//    typedef Array1D<TYPE> TContainerVe;
//
//    /** Row type. */
//    typedef ArrayHo<TYPE> TContainerHo;

    /** Default constructor
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    Array2D( Range const& I = Range(), Range const& J = Range())
           : _ITArray2DType(I, J)
    {
      // initialize vertically the container
      this->initializeCols(J);
    }

    /** constructor with rbeg, rend, cbeg and cend specified,
     *  initialization with a constant.
     *  @param I range of the Rows
     *  @param J range of the Cols
     *  @param v initial value of the container
     **/
    Array2D( Range const& I, Range const& J, TYPE const& v)
           : _ITArray2DType(I, J)
    {
      // initialize vertically the container
      this->initializeCols(J);
      // initialize with v
      for (Integer j=J.first(); j<=J.last(); j++)
      {
        TYPE* p(this->data(j));
        for (Integer i=I.first(); i<=I.last(); i++) p[i]= v;
      }
    }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array2D( const Array2D &T, bool ref=false)
           : _ITArray2DType(T, ref)
    {
      if (!ref)
      {
        // initialize the Cols and Rows
        this->initializeCols(T.rangeHo());
        for (Integer j=T.firstCol(); j<=T.lastCol(); j++)
        {
          // ptr on the rows
          TYPE* p = this->data(j);
          TYPE* pt= T.data(j);
          for (Integer i=T.firstRow(); i<=T.lastRow(); i++) p[i]=pt[i];
        }
      }
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    Array2D( const Array2D& T, Range const& I, Range const& J)
           : _ITArray2DType(T, I, J)
    { ;}

    /** Wrapper constructor Contruct a reference container.
     *  @param q pointer on the data
     *  @param I range of the  Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    Array2D( TYPE** q, Range const& I, Range const& J)
           : _ITArray2DType(q, I, J)
    { ;}

    /** virtual destructor. */
    virtual ~Array2D()
    { ;}

    /** access to a sub-array.
     *  @param I range of the rows
     *  @param J range of the columns
     **/
    inline Array2D sub(Range const& I, Range const& J) const
    { return Array2D(*this, I, J);}

    /** pseudo virtual method required by ITArray2D.
     *  compute the range of the effectively stored elements in the column icol.
     *  @param icol the index of the column we want to know the range
     **/
    inline Range compRangeVe( Integer const& icol) const
    { return this->rangeVe();}

    /** operator = : overwrite the Array2D with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline Array2D& operator=(const Array2D &T)
    {
      // Resize if necessary.
      if ( (this->sizeVe() != T.sizeVe())
         ||(this->sizeHo() != T.sizeHo())
         )
        this->resize(T.rangeVe(), T.rangeHo());
      // Copy without overlapping
      if (T.firstRow()>=this->firstRow())
      {
        if (T.firstCol()>=this->firstCol())
        {
          for ( Integer jt=T.firstCol(), j=this->firstCol()
              ; jt<=T.lastCol()
              ; j++, jt++
              )
          {
            TYPE *p =this->data(j), *pt =T.data(jt);
            for ( Integer it=T.firstRow(), i=this->firstRow()
                ; it<=T.lastRow()
                ; i++, it++
                )
              p[i] = pt[it];
          }
          return *this;
        }
        // T.firstCol()<this->firstCol()
        for ( Integer jt=T.lastCol(), j=this->lastCol()
            ; jt>=T.firstCol()
            ; j--, jt--
            )
        {
          TYPE *p =this->data(j), *pt =T.data(jt);
          for ( Integer it=T.firstRow(), i=this->firstRow()
              ; it<=T.lastRow()
              ; i++, it++
              )
            p[i] = pt[it];
        }
        return *this;
      }
      // T.firstRow()<this->firstRow()
      if (T.firstCol()>=this->firstCol())
      {
        for ( Integer jt=T.firstCol(), j=this->firstCol()
            ; jt<=T.lastCol()
            ; j++, jt++
            )
        {
          TYPE *p =this->data(j), *pt =T.data(jt);
          for ( Integer it=T.lastRow(), i=this->lastRow()
              ; it>=T.firstRow()
              ; i--, it--
              )
            p[i] = pt[it] ;
        }
        return *this;
      }
      // T.firstCol()<this->firstCol()
      for ( Integer jt=T.lastCol(), j=this->lastCol()
          ; jt>=T.firstCol()
          ; j--, jt--
          )
      {
        TYPE *p =this->data(j), *pt =T.data(jt);
        for ( Integer it=T.lastRow(), i=this->lastRow()
            ; it>=T.firstRow()
            ; i--, it--
            )
          p[i] = pt[it];
      }
      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline Array2D& operator=(TYPE const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { TYPE *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] = v;
      }
      return *this;
    }

    /** Swapping the pos1 row and the pos2 row.
     *  @param pos1 position of the first row
     *  @param pos2 position of the second row
     **/
    void swapRows( Integer const& pos1, Integer pos2)
    {
#ifdef STK_BOUNDS_CHECK
      // check conditions
      if (this->firstRow() > pos1)
      { throw out_of_range("ITArray2D::swapRows(pos1, pos2) "
                           "this->firstRow() > pos1");
      }
      if (this->lastRow() < pos1)
      { throw out_of_range("ITArray2D::swapRows(pos1, pos2) "
                           "this->lastRow() < pos1");
      }
      if (this->firstRow() > pos2)
      { throw out_of_range("ITArray2D::swapRows(pos1, pos2) "
                           "this->firstRow() > pos2");
      }
      if (this->lastRow() < pos2)
      { throw out_of_range("ITArray2D::swapRows(pos1, pos2) "
                           "this->lastRow() < pos2");
      }
#endif
      // swap
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      {
        TYPE aux(this->asLeaf().elt(pos1, j));
        this->asLeaf().elt(pos1, j) = this->asLeaf().elt(pos2, j);
        this->asLeaf().elt(pos2, j) = aux;
      }
    }
};

/** @ingroup Arrays
 *  ostream for Array2D.
 *  @param s the output stream
 *  @param V the Array2D to write
 **/
template<class TYPE>
ostream& operator<<(ostream& s, Array2D<TYPE> const& V)
{ return out2D(s,V);}


} // namespace STK

#endif
// STK_ARRAY2D_H
