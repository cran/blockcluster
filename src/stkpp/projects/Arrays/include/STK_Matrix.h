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
 * Purpose:  Define the Matrix classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Matrix.h
  * @brief In this file we specialize the Array2D class for Real values.
 **/

#ifndef STK_MATRIX_H
#define STK_MATRIX_H

#include "STK_Vector.h"
#include "STK_Point.h"

#include "STK_Array2D.h"
#include "../../Sdk/include/STK_ContainerTraits.h"

namespace STK
{

/** @ingroup Arrays
  * @brief Specialization of the Array2D class for Real values.
  * 
  * A Matrix is a column oriented 2D container of Real.
 **/
typedef Array2D<Real> Matrix;

/** @ingroup Arrays
 *  @brief  Specialization of the Array2D class for Real values.
 *
 *  The Matrix class (Array2D<Real> specialization) is a final class.
 *  It would not possible to derive this class and to use the method
 *  @code
 *    clone()
 *    asLeaf()
 *    asPtrLeaf()
 *  @endcode
 *  as the curious recursive template paradigm would be break.
 **/
template<>
class Array2D<Real> : public ITArray2D< Real, Array2D<Real> >
{
  /** Type for the Interface base Class. */
  typedef ITArray2D<Real, Array2D<Real> > _ITArray2DType;

  public:
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
    Array2D( Range const& I, Range const& J, Real const& v)
           : _ITArray2DType(I, J)
    {
      // initialize vertically the container
      this->initializeCols(J);
      // initialize with v
      for (Integer j=J.first(); j<=J.last(); j++)
      {
        Real* p(this->data(j));
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
          Real* p = this->data(j);
          Real* pt= T.data(j);
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
    Array2D( Real** q, Range const& I, Range const& J)
           : _ITArray2DType(q, I, J)
    { ;}

    /** virtual destructor. */
    virtual ~Array2D()
    { ;}

    /** access to one element. */
    inline Real& elt(Integer i, Integer j)
    { return this->data(j)[i];}

    /** access to one element const.  */
    inline Real const& elt(Integer i, Integer j) const
    { return this->data(j)[i];}

    /** access to a sub-array.      */
    inline Array2D sub(Range const& I, Range const& J) const
    { return Array2D(*this, I, J);} 

    /** access to a column.        */
    inline Vector col( Integer j) const
    { return Vector(this->data(j), this->rangeVe(), j);}

    /** access to a sub-col.        */
    inline Vector col(Range const& I, Integer j) const
    { return Vector(this->data(j), I, j);}

    /** access to a row.
     *  @param i index of the row
     *  @return A reference on the row i
     **/
    inline Point row( Integer const& i) const
    { return Point(*this, this->rangeHo(), i);}

    /** access to a sub-row.        */
    inline Point row(Integer i, Range const& J) const
    { return Point(*this, J, i);}

    /** pseudo virtual method required by ITArray2D.
     *  compute the range of the effectively stored elements in the col
     *  icol.
     *  @param icol the index of the column to compute the range
     **/
    inline Range compRangeVe( Integer const& icol) const
    { return this->rangeVe();}

    /** Swapping the pos1 row and the pos2 row.
     *  @param pos1 position of the first row
     *  @param pos2 position of the second row
     **/
    void swapRows( Integer const& pos1, Integer pos2)
    {
#ifdef STK_BOUNDS_CHECK
      // check conditions
      if (this->firstRow() > pos1)
      { throw out_of_range("Array2D::swapRows(pos1, pos2) "
                                "this->firstRow() > pos1");
      }
      if (this->lastRow() < pos1)
      { throw out_of_range("Array2D::swapRows(pos1, pos2) "
                                "this->lastRow() < pos1");
      }
      if (this->firstRow() > pos2)
      { throw out_of_range("Array2D::swapRows(pos1, pos2) "
                                "this->firstRow() > pos2");
      }
      if (this->lastRow() < pos2)
      { throw out_of_range("Array2D::swapRows(pos1, pos2) "
                                "this->lastRow() < pos2");
      }
#endif
      // swap
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      {
        Real aux(this->asLeaf().elt(pos1, j));
        this->asLeaf().elt(pos1, j) = this->asLeaf().elt(pos2, j);
        this->asLeaf().elt(pos2, j) = aux;
      }
    }

    /** operator = : overwrite the Array2D with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline Array2D& operator=(const Array2D &T)
    {
//      // clear any allocated memory
//      this->clear();
//      // set AllocatorBase part
//      this->setRef(true);
//      this->setPtrData(T.ptrData(), T.rangeData());
//      // ITContainer2D part
//      this->setRange(T.rangeVe_, T.rangeHo_);
//      // Array2D part
//      rangeCols_ = T.rangeCols();
//      this->setCapacityHo( T.capacityHo());
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
            Real *p =this->data(j), *pt =T.data(jt);
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
          Real *p =this->data(j), *pt =T.data(jt);
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
          Real *p =this->data(j), *pt =T.data(jt);
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
        Real *p =this->data(j), *pt =T.data(jt);
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
    inline Array2D& operator=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { Real *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] = v;
      }
      return *this;
    }

    /** operator += : add to the container a constant value.
     *  @param v the value to set
     **/
    inline Array2D& operator+=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { Real *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] += v;
      }
      return *this;
    }

    /** operator -= : Subtract to the container a constant value.
     *  @param v the value to set
     **/
    inline Array2D& operator-=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { Real *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] += v;
      }
      return *this;
    }

    /** operator *= : multiply to the container a constant value.
     *  @param v the value to set
     **/
    inline Array2D& operator*=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { Real *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] += v;
      }
      return *this;
    }

    /** operator /= : divide to the container a constant value.
     *  @param v the value to set
     **/
    inline Array2D& operator/=(Real const& v)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      { Real *p =this->data(j);
        for (Integer i=this->firstRow(); i<=this->lastRow(); i++)
          p[i] /= v;
      }
      return *this;
    }

    /** operator = overwriting a Array2D with An Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array2D& operator=(const Exp& rhs)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); ++j)
        for (Integer i=this->firstRow(); i<=this->lastRow(); ++i)
        {
          this->elt(i,j) = rhs(i,j);
        }
      return *this;
    }

    /** operator += Adding an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array2D& operator+=(const Exp& rhs)
    {
      for (Integer j=this->firstCol(); j<=this->lastCol(); ++j)
        for (Integer i=this->firstRow(); i<=this->lastRow(); ++i)
        {
          this->elt(i,j) += rhs(i,j);
        }
      return *this;
    }

    /** operator -= decreasing an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array2D& operator-=(const Exp& rhs)
    {
        for (Integer j=this->firstCol(); j<=this->lastCol(); ++j)
          for (Integer i=this->firstRow(); i<=this->lastRow(); ++i)
          {
            this->elt(i,j) -= rhs(i,j);
          }
      return *this;
    }

    /** operator /= Dividing an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array2D& operator/=(const Exp& rhs)
    {
        for (Integer j=this->firstCol(); j<=this->lastCol(); ++j)
          for (Integer i=this->firstRow(); i<=this->lastRow(); ++i)
          {
            this->elt(i,j) /= rhs(i,j);
          }
      return *this;
    }

    /** operator *= Multiplying an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array2D& operator*=(const Exp& rhs)
    {
        for (Integer j=this->firstCol(); j<=this->lastCol(); ++j)
          for (Integer i=this->firstRow(); i<=this->lastRow(); ++i)
          {
            this->elt(i,j) *= rhs(i,j);
          }
      return *this;
    }
};

} // namespace STK

#endif
// STK_MATRIX_H
