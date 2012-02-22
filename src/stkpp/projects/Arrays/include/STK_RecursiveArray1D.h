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
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_RecursiveArray1D.h
  * @brief Pre-Implementation of the interface IArray1DBase.
  * 
  * An RecursiveArray1D is a non-oriented templated one dimensional Array
  * which is not final and can be used by any sub-class that need to be
  * a final class.
 **/

#ifndef STK_RECURSIVEARRAY1D_H
#define STK_RECURSIVEARRAY1D_H

#include "STK_IArray1DBase.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Templated one dimensional Array.
 * 
 * An RecursiveArray1D is a templated non-oriented container implementing the interface
 * base class IArray1DBase. It contains objects of the type @c TYPE.
 **/
template<class TYPE, class Container1D >
class RecursiveArray1D : public IArray1DBase<TYPE, TYPE, Container1D >
{
  protected:
    /** Type for the Base reference Class of a 2D Array. */
    typedef AllocatorBase<TYPE*> _AllocatorBaseType_;
    /** Type for the IArray1DBase Class.                              */
    typedef IArray1DBase<TYPE, TYPE, Container1D > _IArray1DType;

  public:
    /** Default constructor : beg_ =1 and end_ =0.
     *  @param I range of the container
     **/
    RecursiveArray1D( Range const& I = Range())
                    : _IArray1DType(I)
    { ;}

    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    RecursiveArray1D( Range const& I, TYPE const& v)
                   : _IArray1DType(I)
    {
      for (Integer i=this->first(); i<=this->last(); i++)
        this->setData(i, v);
    }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    RecursiveArray1D( const RecursiveArray1D &T, bool ref =false)
                   : _IArray1DType(T, ref)
    {
      // check if we want just a reference
      if (!ref)
      {
        for (Integer j=this->first(); j<=this->last(); j++)
          this->setData(j, T.data(j));
      }
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the data to wrap
     **/
    RecursiveArray1D( RecursiveArray1D const& T, Range const& I)
                   : _IArray1DType(T, I)
    { ;}

    /** Wrapper constructor : the container is a reference.
     *  @param q pointer on data
     *  @param I range of the data
     *  @param index index (row or column) of the data
     **/
    RecursiveArray1D( TYPE* q, Range const& I, Integer const& index =0)
                   : _IArray1DType(q, I, index)
    { ;}

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param index the index of the column to wrap
     **/
     RecursiveArray1D( _AllocatorBaseType_ const& T, Range const& I, Integer const& index)
                    : _IArray1DType(T.data(index), I, index)
    { ;}

   /** Virtual destructor: allocated memory is liberated by AllocatorBase base
    *  class.
    **/
    virtual ~RecursiveArray1D() { }

    /** Clear the object. Memory is liberated and the
     *  range of the Container is set to 0:-1.
     **/
    void clear()
    {
      if (this->isRef()) return;   // Nothing to do for ref
      this->freeMem();  // Free Mem
      this->setRange(); // Set dimension to default
    }

    /** Add n Elts to the container.
     *  @param n number of elements to add
     **/
    virtual void pushBack( Integer const& n=1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
#ifdef STK_DEBUG
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("Array1D::pushBack(n) "
                                 "can't operate on references.");
      }
#endif
      // If the container is empty : create it
      if (this->empty())
        this->initialize(Range(this->first(), this->first()+n-1));
      else
        this->insertElt(this->last()+1, n);
    }

    /** Delete last elts of the container.
     *  @param n number of elts to delete
     **/
    virtual void popBack(Integer const& n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
#ifdef STK_DEBUG
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("Array1D::popBack() "
                            "can't operate on references.");
      }
#endif
#ifdef STK_BOUNDS_CHECK
      // if there is elts to erase
      if (this->size()<n)
      { throw out_of_range("Array1D::popBack(n) "
                           "this->size() < n");
      }
#endif
      // update range
      this->decLast(n);
      // if there is no more elts
      if (this->size() == 0) this->freeMem();
    }

    /** Delete n elements at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    virtual void erase(Integer const& pos, Integer const& n=1)
    {
      // if n==0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("Array1D::erase(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->first() > pos)
      { throw out_of_range("Array1D::erase(pos, n) "
                           "this->first() > pos");
      }
      if (this->last() < pos)
      { throw out_of_range("Array1D::erase(pos, n) "
                           "this->last() < pos");
      }
      if (this->last() < pos+n-1)
      { throw out_of_range("Array1D::erase(pos, n) "
                           "this->last() < pos+n-1");
      }
#endif
      // translate remaining elts
      const Integer last = this->last()-n;
      for (Integer k=pos; k<=last; k++)
        this->setData(k, this->data(k+n));
      // update dimensions
      this->decLast(n);
      // if there is no more cols, free mem
      if (this->size() == 0) this->freeMem();
    }

  protected:
    /** function for memory allocation and initialization.
     *  This method will free all allocated memory owned by this
     *  container before initialization.
     *  @param I range of the container
     **/
    void initialize(Range const& I)
    {
      // check if there is memory allocated
      this->clear();
      // if we initialize the memory the container is not a reference
      this->setRef(false);
      // try to allocate memory
      this->init1D(I);
      // set the range of the container if init1D is successful
      this->setRange(I);
    }

    /** Method for memory deallocation. Memory is liberated and the
     *  range of the Container is set to begin:begin-1.
     **/
    void freeMem()
    {
      if (this->isRef()) return;   // Nothing to do for ref
      this->free1D();              // free the elts
    }

    /** Insert n elts at the position pos of the container. The bound
     *  last_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    virtual void insertElt( Integer const& pos, Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("Array1D::insertElt(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check indices
      if (this->first() > pos)
      { throw out_of_range("Array1D::insertElt(pos, n) "
                                "this->first() > pos");
      }
      if (this->last()+1 < pos)
      { throw out_of_range("Array1D::insertElt(pos, n) "
                                "this->last()+1 < pos");
      }
#endif
      // allocate, if necessary, the mem for the elts
      if (this->capacity() < this->size()+n)
      {
        // temporary empty container
        RecursiveArray1D Taux;
        // save elts
        this->swap(Taux);
        // compute range of the container after insertion
        Range range(Taux.range());
        range.incLast(n);
        // initialize
        try
        {
          this->init1D(range);
        }
        catch (runtime_error & error)   // if an error occur
        {
          this->swap(Taux); // restore elts
          throw error;      // and send again the Exception
        }
        // reset initial stored in range
        this->setRange(Taux.range());
        // copy first elts
        for (Integer k=this->first(); k<pos; k++)
          this->setData(k, Taux.data(k));
        // translate and copy last elts
        for (Integer k=this->last(); k>=pos; k--)
          this->setData(k+n, Taux.data(k));
      }
      else // enough space -> shift the last elts
      {
        // translate data
        for (Integer k=this->last(); k>=pos; k--)
          this->setData(k+n, this->data(k));
      }
      // update range_
      this->incLast(n);
    }
};

} // namespace STK

#endif
// STK_RECURSIVEARRAY1D_H
