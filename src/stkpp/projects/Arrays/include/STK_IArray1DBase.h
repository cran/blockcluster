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

    Contact : Serge.Iovleff@stkpp.or
*/

/*
 * Project: stkpp::Arrays
 * Purpose:  Define the Interface Base class for the Array1D classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IArray1DBase.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 * 
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array1D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY1DBASE_H
#define STK_IARRAY1DBASE_H

#include "../../STKernel/include/STK_Misc.h"

#include "../../Sdk/include/STK_ITContainer1D.h"
#include "STK_AllocatorBase.h"

namespace STK
{

/** @ingroup Arrays
 * @brief Templated Interface base class for uni-dimensional arrays.
 *
 * AIArray1DBase is an abstract non-oriented uni-dimensional
 * container. It inherits of all pseudo virtual methods defined in
 * ITContainer1D and IContainer1D Interfaces.
 * 
 * Description of the template parameters :
 * - @c TYPE is the type of the data stored in the container.
 * - @c PTRELT allow to use indirection when accessing to the data.
 * - @c TArray1D is the type of the class which derives from
 *  IArray1DBase. For exemple :
 * <code>
 * template<class TYPE>
 * class TArray1D : public IArray1DBase<TYPE, PTRELT, TArray1D<TYPE> >
 * {...}
 * </code>
 **/
template< class TYPE
        , class PTRELT
        , class TArray1D
        >
class IArray1DBase : public ITContainer1D<TYPE, TArray1D>
                   , public AllocatorBase<PTRELT>
{                                                                      
  /** Type for the Array Base Class. */
  typedef AllocatorBase<PTRELT> _AllocatorBaseType_;
  /** Type for the Interface Array1D Class. */
  typedef ITContainer1D<TYPE, TArray1D> _ITContainer1DType_;

  protected:
    /** Default constructor : beg_ =1 and end_ =0.
     *  @param I range of the container
     **/
    IArray1DBase( Range const&  I = Range())
                : _ITContainer1DType_(I)
                , _AllocatorBaseType_()
                , capacity_(0)
                , index_(0)
    { init1D(I);}

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    IArray1DBase( const IArray1DBase& T, bool ref =false)
                : _ITContainer1DType_(T)
                , _AllocatorBaseType_(T, ref)
                , capacity_(0)
                , index_(T.index_)
    {
      if (!ref)
        init1D(T.range());
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I range of the data to wrap
     **/
    IArray1DBase( const IArray1DBase& T, Range const& I)
                : _ITContainer1DType_(I)
                , _AllocatorBaseType_(T, true)
                , capacity_(0)
                , index_(T.index_)
    {
#ifdef STK_BOUNDS_CHECK
      if (I.first() < T.first())
      { throw out_of_range("IArray1DBase::IArray1DBase(T, I) "
                           "I.first() < T.first()");
      }
      if (I.last() > T.last())
      { throw out_of_range("IArray1DBase::IArray1DBase(T, I) "
                           "I.last() > T.last()");
      }
#endif
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I range of the data to wrap
     *  @param index index (row or col) of the container
     **/
    IArray1DBase( const _AllocatorBaseType_& T, Range const& I, Integer const& index)
                : _ITContainer1DType_(I)
                , _AllocatorBaseType_(T, true)
                , capacity_(0)
                , index_(index)
    { ;}

    /** Wrapper constructor
     *  @param q pointer on data
     *  @param I range of the data
     *  @param index index (row or col) of the data
     **/
    IArray1DBase( PTRELT* q, Range const& I, Integer const& index=0)
           : _ITContainer1DType_(I)
           , _AllocatorBaseType_(q, I)
           , capacity_(0)
           , index_(index)
    { ;}

  public:
    /** Virtual destructor.                                                 */
    virtual ~IArray1DBase() { ;}

    /** return the maximum possible number of elements without
     *  reallocation. */
    inline Integer const& capacity() const { return capacity_;}

    /** Get the index of the oriented container */
    inline Integer const& getIndex() const { return index_;}

    /** reserve internal memory for at least size elements.
     *  @param size number of elements to reserve
     **/
    void reserve(Integer const& size)
    {
      // nothing to do
      if (size < this->capacity()) return;
      // is this structure a ptr ?
      if (this->isRef())
      { throw runtime_error("Array1D::reserve(size) "
                            "can't operate on references.");
      }
      // reserve
      TArray1D Taux;    // auxiliary empty container
      this->swap(Taux); // save old elts in Taux
      try
      {
        this->mallocPtrData(size, Taux.first());
      }
      catch (runtime_error & error) // if an alloc error occur
      {
        // go back
        this->swap(Taux);
        // and send the Exception
        throw error;
      }
      // if no alloc error
      this->setCapacity(size);         // update size
      this->setRange(Taux.range()); // set range
      // copy elts
      for (Integer j=this->first(); j<=this->last(); j++)
      {  
        this->setData(j, Taux.data(j));
        Taux.setData(j, PTRELT());
      }
    }

    /** swap this Container with T.
     *  @param T the Array to swap with T
     **/
    void swap(IArray1DBase &T)
    {
      // swap AllocatorBase part
      this->AllocatorBase<PTRELT>::swap(T);
 
      // swap IArray1DBase part
      STK::swap(capacity_, T.capacity_);
      STK::swap(index_, T.index_);

      // swap ITContainer1D part
      Range auxRange(this->range());
      this->setRange(T.range());
      T.setRange(auxRange);
    }

    /** New beginning index for the object.
     *  @param beg the index of the first column to set
     **/
    void shift(Integer const& beg =1)
    {
      // compute increment
      Integer inc = beg - this->first();
      // if there is something to do
      if (inc != 0)
      {
        // is this structure just a pointer?
        if (this->isRef())
        { throw runtime_error("IArray1DBase::shift(beg) "
                                   "can't operate on references.");
        }
        // translate rangeHo_
        this->incRange(inc);
        // translate data
        this->shiftPtrData(beg);
      }
    }

    /** Swapping the pos1 elt and the pos2 elt.
     *  @param pos1 position of the first elt
     *  @param pos2 position of the second elt
     **/
    void swap(Integer const& pos1, Integer const& pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->first() > pos1)
      { throw out_of_range("IArray1DBase::swap(pos1, pos2) "
                           "this->first() >pos1");
      }
      if (this->last() < pos1)
      { throw out_of_range("IArray1DBase::swap(pos1, pos2) "
                           "this->last() <pos1");
      }
      if (this->first() > pos2)
      { throw out_of_range("IArray1DBase::swap(pos1, pos2) "
                           "this->first() >pos2");
      }
      if (this->last() < pos2)
      { throw out_of_range("IArray1DBase::swap(pos1, pos2) "
                           "this->last() <pos2");
      }
#endif
      // swap
      PTRELT qaux(this->data(pos1));
      this->setData(pos1, this->data(pos2));
      this->setData(pos2, qaux);
    }

    /** Insert element @c v in the range @c I of the Array.
     *  @param I range of the index where to insert elements
     *  @param v the value to insert
     **/
    void insert( Range const& I, TYPE const& v)
    {
      this->insertElt(I.first(), I.size());
      for (Integer i=I.first(); i<=I.last(); i++)
        this->elt(i) = v;
    }

    /** STL compatibility : push front an element.
     *  @param v value to append
     **/
    inline void push_front(TYPE const& v)
    { insert(Range(this->first(), this->first()), v);}

  protected:
    /** Set the maximum possible number of elements without
     *  reallocation.
     *  @param capacity capacity of the container
     **/
    inline void setCapacity(Integer const& capacity =0)
    { capacity_ = capacity;}

    /** Set the index (row or col) of the container.
     *  @param index index of the container
     **/
    inline void setIndex(Integer const& index =0)
    { index_ = index;}

    /** function for memory allocation and initialization.
     *  The range is not set in this method. If an
     *  error occur, we set the range of the container to default.
     *  @param I range of the container
     **/
    void init1D(Range const& I)
    {
      // compute the size necessary (can be 0)
      Integer size = _AllocatorBaseType_::evalCapacity(I.size());
      // try to allocate memory
      try
      {
        // initialize Elts
        this->mallocPtrData(size, I.first());
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0) 
        setCapacity();
        // set default range
        this->setRange(); 
        // throw the error
        throw error;
      }
      // set new capacity if no error occur
      setCapacity(size);
    }

    /** Method for memory deallocation. If the derived class
     *  use indirection, we have to free the mem if necessary prior
     *  to this method. The beginning of the Container is not modified
     **/
    void free1D()
    { 
      // Nothing to do for ref
      if (this->isRef()) return;
      // free allocated memory
      this->freePtrData();
      // set capacity to default
      this->setCapacity(); 
      // set range of the Cols to default
      this->setRange(Range(this->first(), this->first()-1));
    }

    /** Move the column pos2 from the container T to the column pos1 of
     *  this. set the column pos2 in T to default value.
     *  @param pos1 the column to initialize
     *  @param T the container with the column to move
     *  @param pos2 the column in the container T to move in this
     **/
    void moveElt( Integer const& pos1, IArray1DBase& T
                , Integer const& pos2
                )
    {
      // copy column pos2 of T in pos1 of this
      this->setData(pos1, T.data(pos2));
      // set column of T to default
      T.setData(pos2, PTRELT());
    }

  private:
    /** capacity of the array. */
    Integer capacity_;
    /** row or column of the container (if any). */
    Integer index_;
};

} // namespace STK

#endif
// STK_IARRAY1DBASE_H
