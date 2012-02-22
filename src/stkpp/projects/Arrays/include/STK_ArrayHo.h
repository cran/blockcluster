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
 * Purpose:  Define the ArrayHo class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ArrayHo.h
  * @brief A ArrayHo is a one dimensional horizontal container
  * 
  * An ArrayHo is an implementation of the interface IArray1DBase.
  * It's a one dimensional horizontal container.
 **/

#ifndef STK_ARRAYHO_H
#define STK_ARRAYHO_H

#include "STK_IArray1DBase.h"
#include "STK_Display1D.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Templated one dimensional horizontal Array. 
 * 
 *  A ArrayHo is a horizontal container. In STK, the objects are
 *  vertical (column) or horizontal (rows).
 * 
 *  By default the index of the first element is 1 but this can be
 *  modified using the appropriate constructor or using the method @c shift.
 **/
template<class TYPE>
class ArrayHo : public IArray1DBase<TYPE, TYPE*, ArrayHo<TYPE> >
{
  public:
    /** Type for the Base reference Class.                            */
    typedef AllocatorBase<TYPE*> _AllocatorBaseType_;

    /** Type for the Implementation Class.                            */
    typedef IArray1DBase<TYPE, TYPE*, ArrayHo<TYPE> > _IArrayHoType;

    /** Default constructor : beg_ =1 and end_ =0.
     *  @param I range of the container
     **/
    ArrayHo( Range const& I = Range())
           : _IArrayHoType(I)
    { this->initElts(I);}

    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    ArrayHo( Range const& I, TYPE const& v)
           : _IArrayHoType(I)
    {
      this->initElts(I);
      // set value v
      for (Integer i=this->first(); i<=this->last(); i++)
        this->elt(i) = v;
    }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if this is a wrapper of T
     **/
    ArrayHo( const ArrayHo &T, bool ref =false)
           : _IArrayHoType(T, ref)
    {
      // check if we want just a reference
      if (!ref)
      {
        // if this is not a reference, initialize the container
        this->initElts(T.range());
        // this and T have the same row : save time
        const Integer row(T.getIndex());
        // and copy the data
        for (Integer j=this->first(); j<=this->last(); j++)
          this->data(j)[row] = T.data(j)[row];
      }
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     **/
    ArrayHo( const ArrayHo<TYPE>& T, Range const& I)
           : _IArrayHoType(T, I)
    { ;}

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param row the index of the row to wrap
     **/
     ArrayHo( const _AllocatorBaseType_& T, Range const& I, Integer row)
            : _IArrayHoType(T, I, row)
    { ;}

    /** virtual destructor.                                                 */
    virtual ~ArrayHo()
    { if (!this->isRef())
      { freeElts(this->range());}
    }

    /** Get one element const.
     *  @param pos index of the element (const)
     **/
    inline TYPE const & elt(Integer const& pos) const
    { return this->data(pos)[this->getIndex()];}

    /** Get one element.
     *  @param pos index of the element
     **/
    inline TYPE& elt(Integer const& pos)
    { return this->data(pos)[this->getIndex()];}

    /** access to many elements.
     *  @param J Range of the elements
     **/
    inline ArrayHo elt(Range const& J) const
    { return ArrayHo(*this, J, this->getIndex());}

    /** clear the object.                                             */
    void clear()
    {
      this->freeMem();  // Free Mem
      this->setRange(); // Set dimension to default
    }

    /** Method for memory deallocation. Memory is liberated but the
     *  range of the container is not updated.
     **/
    void freeMem()
    {
      if (this->isRef()) return;        // Nothing to do for ref
      this->freeElts(this->range()); // free the ptr of elts
      this->free1D();                   // free the elts
    }

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
      // initialize Horizontally the container
      this->initElts(I);
      // set the range of the container if init1D is successful
      this->setRange(I);
    }

    /** Add n elements to the container.
     *  @param n number of elements to add
     **/
    void pushBack( Integer const& n=1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ArrayHo::pushBack(n) "
                                 "can't operate on references.");
      }
      // If the container is empty : create it
      if (this->empty())
        this->initialize(Range(this->first(), this->first()+n-1));
      else
        this->insertElt(this->last()+1, n);
    }

    /** Delete last elts of the container.
     *  @param n number of elts to delete
     **/
    void popBack(Integer const& n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ArrayHo::popBack() "
                            "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // if there is elts to erase
      if (this->size()<n)
      { throw out_of_range("ArrayHo::popBack(n) "
                           "this->size() < n");
      }
#endif
      // delete each elt
      this->freeElts(Range(this->last()-n+1, this->last()));
      // update range
      this->decLast(n);
      // if there is no more elements
      if (this->size() == 0) this->freeMem();
    }

    /** Delete n elts at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void erase( Integer const& pos, Integer const& n=1)
    {
      // if n==0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ArrayHo::erase(pos, n) "
                            "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->first() > pos)
      { throw out_of_range("ArrayHo::erase(pos, n) "
                           "this->first() > pos");
      }
      if (this->last() < pos)
      { throw out_of_range("ArrayHo::erase(pos, n)"
                           " this->last() < pos");
      }
      if (this->last() < pos+n-1)
      { throw out_of_range("ArrayHo::erase(pos, n)"
                           " this->last() < pos+n-1");
      }
#endif
      // delete each col
      freeElts(Range(pos, pos+n-1));
      // shift elements
      const Integer last = this->last()-n;
      for (Integer k=pos; k<=last; k++)
        this->setData(k, this->data(k+n));
      // update rangeHo_
      this->decLast(n);
      // if there is no more elements
      if (this->size() == 0) this->freeMem();
    }
 
    /** operator = : overwrite the ArrayHo with the ArrayHo T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    inline ArrayHo& operator=(const ArrayHo &T)
    { // We have to resize if this and T does not have the same size
      // but if they have the same size, we don't scale the index
      if (this->size()!=T.size()) this->resize(T.range());

      // copy without ovelapping
      if (this->first() < T.first())
      { for (Integer j=this->first(), i=T.first(); i<=T.last(); j++, i++)
          elt(j) = T.elt(i);
      }
      else
      { for (Integer j=this->last(), i=T.last(); i>=T.first(); j--, i--)
          elt(j) = T.elt(i);
      }
      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    inline ArrayHo& operator=(TYPE const& v)
    {
      for (Integer i=this->first(); i<=this->last(); i++)
        this->elt(i) = v;
      return *this;
    }

    /** operator = : overwrite the ArrayHo with the container T.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
    template <class LEAF>
    inline ArrayHo& operator=(const ITContainer1D<TYPE, LEAF> &T)
    { // We have to resize if this and T does not have the same size
      // but if they have the same size, we don't scale the index
      if (this->size()!=T.size()) this->resize(T.range());

      // copy without ovelapping
      if (this->first() < T.first())
      { for (Integer j=this->first(), i=T.first(); i<=T.last(); j++, i++)
          elt(j) = T[i];
      }
      else
      { for (Integer j=this->last(), i=T.last(); i>=T.first(); j--, i--)
          elt(j) = T[i];
      }
      return *this;
    }

  protected:
    /** Insert n elts at the position pos of the container. The bound
     *  last_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    void insertElt(Integer const& pos, Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ArrayHo::insertElt(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check indices
      if (this->first() > pos)
      { throw out_of_range("ArrayHo::insertElt(pos, n) "
                                "this->first() > pos");
      }
      if (this->last()+1 < pos)
      { throw out_of_range("ArrayHo::insertElt(pos, n) "
                                "this->last()+1 < pos");
      }
#endif
      // allocate, if necessary, the mem for the elts
      if (this->capacity() < this->size()+n)
      {
        // compute range of the container after insertion
        Range range(this->range());
        range.incLast(n);
        // temporary empty container
        ArrayHo Taux;
        // save elements in Taux
        this->swap(Taux);
        // initialize elements of the container
        try
        {
          this->init1D(range);
        }
        catch (runtime_error & error)   // if an error occur
        {
          this->swap(Taux);   // restore container
          throw error;        // and send again the Exception
        }
        // reset initial stored in range
        this->setRange(Taux.range());
        // move first elements from Taux to this
        for (Integer k=this->first(); k<pos; k++)
          this->moveElt(k,  Taux, k);
        // translate and copy last elements from Taux to this
        for (Integer k=this->last(); k>=pos; k--)
          this->moveElt(k+n, Taux, k);
      }
      else // enough space -> initialize and shift the last elts
      {
        // translate data
        for (Integer k=this->last(); k>=pos; k--)
          this->moveElt(k+n, *this, k);
      }
      // initialize the elements in the range pos:pos+n-1
      initElts(Range(pos, pos+n-1));
      // update range
      this->incLast(n);
    }
  private:
    /** Method for elements memory allocation.
     *  @param J the range of the elements to initialize
     **/
    void initElts(Range const& J)
    {
      // for each col
      for (Integer j=J.first(); j<=J.last(); j++)
      {
        // try to Allocate mem for the jth elt
        try
        {
          this->initElt(j);
        }
        catch (runtime_error & error)   // if an error occur
        {
          // for each column allocated
          for (Integer k=J.first(); k<j; k++)
            this->freeElt(k);
          // put default parameters for the elts j to end
          for (Integer k=j; k<=J.last(); k++)
            this->setData(k, 0);
          // and throw an Exception
          throw error;
        }
      }
    }
    /** function for the the allocation of memory of element pos.
     *  @param pos the number of the elt to initialize
     **/
    void initElt(Integer const& pos)
    {
      // try to Allocate mem for each col
      try
      {
        this->setData(pos, new TYPE);
        this->data(pos) -= this->getIndex();
      }
      catch (std::bad_alloc & error)        // if an alloc error occur
      {
        // set default
        this->setData(pos, 0);
        // throw an Exception
        throw runtime_error("ArrayHo::initElt(pos) "
                            "memory allocation failed.");
      }
    }

    /** Method for elements memory deallocation.
     *  @param pos the position of the elements to liberate
     **/
    void freeElt(Integer const& pos)
    {
      if (this->data(pos))  // if there is an elt
      {
        this->data(pos) += this->getIndex();
        delete this->data(pos);  // erase
        this->setData(pos, 0);         // set default
      }
    }

    /** protected function for the the deallocation of memory of the
     *  element in the given Range.
     *  @param J the range of the elements to liberate
     **/
    void freeElts(Range const& J)
    { 
      // for all elts
      for (Integer j=J.first(); j<=J.last(); j++)
        this->freeElt(j);
    }
};

/** @ingroup Arrays
 *  ostream for ArrayHo.
 *  @param s the output stream
 *  @param V the ArrayHo to write
 **/
template<class TYPE>
ostream& operator<<(ostream& s, const ArrayHo<TYPE>& V)
{ return out1D(s,V);}


} // namespace STK

#endif // STK_ARRAYHO_H
