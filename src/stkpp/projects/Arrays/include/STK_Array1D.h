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
 * created on: 26 nov. 2007
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Array1D.h
  * @brief Implementation of the final class Array1D
  * 
  * An Array1D is a non-oriented templated uni-dimensional Array.
 **/

#ifndef STK_ARRAY1D_H
#define STK_ARRAY1D_H

#include "STK_RecursiveArray1D.h"
#include "STK_Display1D.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Templated one dimensional Array.
 * 
 * An Array1D is a templated non-oriented container implementing the interface
 * base class RecursiveArray1D. It contains objects of the type @c TYPE. It is a
 * final class for the curious recursive paradigm.
 **/
template<class TYPE >
class Array1D : public RecursiveArray1D< TYPE, Array1D<TYPE> >
{
    /** Type for the Base reference Class. */
    typedef AllocatorBase<TYPE*> _AllocatorBaseType_;
    /** Type for the IArray1DBase Class. */
    typedef RecursiveArray1D<TYPE, Array1D<TYPE> > _RecArray1DType_;

  public:
    /** Default constructor : beg_ =1 and end_ =0.
     *  @param I range of the container
     **/
    Array1D( Range const& I = Range())
           : _RecArray1DType_(I)
    { ;}

    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    Array1D( Range const& I, TYPE const& v)
           : _RecArray1DType_(I, v)
    { }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array1D( const Array1D &T, bool ref =false)
           : _RecArray1DType_(T, ref)
    { }

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the data to wrap
     **/
    Array1D( Array1D const& T, Range const& I)
           : _RecArray1DType_(T, I)
    { ;}

    /** Wrapper constructor : the container is a reference.
     *  @param q pointer on data
     *  @param I range of the data
     *  @param index index (row or column) of the data
     **/
    Array1D( TYPE* q, Range const& I, Integer const& index =0)
           : _RecArray1DType_(q, I, index)
    { ;}

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param index the index of the columnto wrap
     **/
     Array1D( _AllocatorBaseType_ const& T, Range const& I, Integer const& index)
            : _RecArray1DType_(T, I, index)
    { ;}

   /** Virtual destructor: allocated memory is liberated by AllocatorBase base class. */
    virtual ~Array1D() { ;}

    /** access to one element.
     *  @param pos index of the element
     **/
    inline TYPE& elt(Integer const& pos)
    { return this->data(pos);}

    /** access to one element const.
     *  @param pos index of the const element
     **/
    inline TYPE const& elt(Integer const& pos) const
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
      return Array1D(this->asLeaf(), J);
    }

    /** operator = : overwrite the Array1D with t.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline Array1D<TYPE>& operator=(const Array1D<TYPE> &T)
    {
      // check size
      if (this->size()!=T.size()) this->resize(T.range());

      // copy without ovelapping.
      const Integer first = this->first(), last = this->last();
      if (first < T.first())
      { for (Integer i=first, j=T.first(); i<=last; i++, j++)
          this->setData(i, T.data(j));
      }
      else
      { for (Integer i=last, j=T.last(); i>=first; i--, j--)
          this->setData(i, T.data(j));
      }
      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline Array1D& operator=(TYPE const& v)
    {
      const Integer first = this->first(), last = this->last();
      for (Integer i=first; i<=last; i++)
        this->setData(i, v);
      return *this;
    }

    /** operator = : overwrite the ArrayHo with the ITContainer1D T.
     *  We resize @c this if this and T don't have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    template <class LEAF>
    inline Array1D& operator=(const ITContainer1D<TYPE, LEAF> &T)
    { // We have to resize if this and T does not have the same size
      // but if they have the same size, we don't scale the index
      if (this->size()!=T.size()) this->resize(T.range());

      // copy without ovelapping
      const Integer first = this->first(), last = this->last();
      if (first < T.first())
      {
        for (Integer j=first, i=T.first(); j<=last; j++, i++)
          this->setData(j, T[i]);
      }
      else
      {
        for (Integer j=last, i=T.last(); j>=first; j--, i--)
          this->setData(j, T[i]);
      }
      return *this;
    }
};

/** @ingroup Arrays
 *  ostream for Array1D.
 *  @param s the output stream
 *  @param V the Array1D to write
 **/
template<class TYPE>
ostream& operator<<(ostream& s, const Array1D<TYPE>& V)
{ return out1D(s,V);}

} // namespace STK

#endif
// STK_ARRAY1D_H
