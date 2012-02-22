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
 * Project:  Astkpp::rrays
 * Purpose:  Define the Vector class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/


/** @file STK_Vector.h
  * @brief In this file we specialize the templated class Array1D
  * for the type Real.
 **/

#ifndef STK_VECTOR_H
#define STK_VECTOR_H

#include "../../STKernel/include/STK_Real.h"
#include "STK_Array1D.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Real one dimensional Vertical Array.
 *
 *  A Vector is a column oriented 1D container of Real. It supports
 *  templates expression or partial evaluation at compile-time.
 **/
typedef Array1D<Real> Vector;

/** @ingroup Arrays
 *  @brief Specialization of the templated class @c Array1D for Real. 
 * 
 *  A @c Vector is a final Array1D class. It is a vertical (column) container.
 **/
template<>
class Array1D<Real> : public RecursiveArray1D< Real, Array1D<Real> >
{
  /** Type for the Base reference Class.                            */
  typedef AllocatorBase<Real*> _AllocatorBaseType_;
  /** Type for the IArray1DBase Class.                              */
  typedef RecursiveArray1D<Real, Array1D<Real> > _RecArray1DType_;

  public:
    /** Default constructor : first_ =1 and last_ =0.
     *  @param I range of the container
     **/
    Array1D( Range const& I = Range())
          : _RecArray1DType_(I)
    { ;}

    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    Array1D( Range const& I, Real const& v)
           : _RecArray1DType_(I)
    {
      for (Integer i=first(); i<=last(); i++)
        setData(i, v);
    }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array1D( const Array1D &T, bool ref =false)
           : _RecArray1DType_(T, ref)
    {
      // check if we want just a reference
      if (!ref)
      {
        for (Integer j=first(); j<=last(); j++)
          setData(j, T.data(j));
      }
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I range of the data to wrap
     **/
    Array1D( Array1D const& T, Range const& I)
           : _RecArray1DType_(T, I)
    { ;}

    /** Wrapper constructor : the container is a reference.
     *  @param q pointer on data
     *  @param I range of the data
     *  @param index index (row or col) of the data
     **/
    Array1D( Real* q, Range const& I, Integer const& index =0)
           : _RecArray1DType_(q, I, index)
    { ;}

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param index the index of the column to wrap
     **/
     Array1D( const _AllocatorBaseType_& T, Range const& I, Integer const& index)
            : _RecArray1DType_(T.data(index), I, index)
    { ;}

    /** Virtual destructor: allocated memory is liberated by AllocatorBase
     *  base class.
     **/
    virtual ~Array1D() { ;}

    /** access to one element.
     *  @param pos index of the element
     **/
    inline Real& elt(Integer const& pos)
    { return this->data(pos);}

    /** access to one element const.
     *  @param pos index of the const element
     **/
    inline Real const& elt(Integer const& pos) const
    { return this->data(pos);}

    /** access to many elements.
     *  @param J the range of the elements
     **/
    inline Array1D elt(Range const& J) const
    {
#ifdef STK_BOUNDS_CHECK
      if ((J.first()<this->first()))
      { throw out_of_range("Array1D::elt(J) "
                                "J.first()<this->first()");
      }
      if ((J.last()>this->last()))
      { throw out_of_range("Array1D::elt(J) "
                                "J.last()>this->last()");
      }
#endif
      return Array1D(*this, J);
    }

    /** operator = : overwrite the Array1D with T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline Array1D& operator=(Array1D const& T)
    {
      // check size
      if (size()!=T.size()) resize(T.range());
      // copy without ovelapping.
      if (first() < T.first())
      { for (Integer i=first(), j=T.first(); j<=T.last(); i++, j++)
          setData(i, T.data(j));
      }
      else
      { for (Integer i=last(), j=T.last(); j>=T.first(); i--, j--)
          setData(i, T.data(j));
      }
      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline Array1D& operator=(Real const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        setData(i, v);
      return *this;
    }

    /** operator+= Adding a constant.
     *  @param v the value to add
     **/
    inline Array1D& operator+=(Real const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) += v;
      return *this;
    }

    /** operator -= decreasing a constant.
     *  @param v the value to decrease
     **/
    inline Array1D& operator-=(Real const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) -= v;
      return *this;
    }

    /** operator /= Dividing a constant.
     *  @param v the value to divide
     **/
    inline Array1D& operator/=(Real const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) /= v;
      return *this;
    }

    /** operator *= Multiplying a constant.
     *  @param v the value to multiply
     **/
    inline Array1D& operator*=(Real const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) *= v;
      return *this;
    }

    /** operator = : overwrite the Array1D with the container T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    template <class LEAF>
    inline Array1D& operator=(const ITContainer1D<Real, LEAF> &T)
    { // We have to resize if this and T does not have the same size
      // but if they have the same size, we don't scale the index
      if (size()!=T.size()) resize(T.range());

      // copy without ovelapping
      if (first() < T.first())
      { for (Integer j=first(), i=T.first(); i<=T.last(); j++, i++)
          setData(j, T[i]);
      }
      else
      { for (Integer j=last(), i=T.last(); i>=T.first(); j--, i--)
          setData(j, T[i]);
      }
      return *this;
    }

    /** operator = overwriting a Array1D with An Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array1D& operator=(const Exp& rhs)
    {
      for (Integer i=first(); i<=last(); i++)
        setData(i, rhs[i]);
      return *this;
    }

    /** operator += Adding an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array1D& operator+=(const Exp& rhs)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) += rhs[i];
      return *this;
    }

    /** operator -= decreasing an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array1D& operator-=(const Exp& rhs)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) -= rhs[i];
      return *this;
    }

    /** operator /= Dividing an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array1D& operator/=(const Exp& rhs)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) /= rhs[i];
      return *this;
    }

    /** operator *= Multiplying an Expression.
     *  @param rhs the right hand side expression
     **/
    template< class Exp>
    inline Array1D& operator*=(const Exp& rhs)
    {
      for (Integer i=first(); i<=last(); i++)
        elt(i) *= rhs[i];
      return *this;
    }
};

} // namespace STK

#endif
// STK_VECTOR_H
