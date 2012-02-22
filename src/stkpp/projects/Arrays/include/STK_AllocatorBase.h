/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004  Serge Iovleff

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
 * Purpose:  Define the Base Interface for the Array classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_AllocatorBase.h
  * @brief In this file we define the class AllocatorBase
 **/

#ifndef STK_ARRAYBASE_H
#define STK_ARRAYBASE_H

#include "../../Sdk/include/STK_IContainerRef.h"
#include "../../STKernel/include/STK_Range.h"
#include "../../STKernel/include/STK_Exceptions.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Templated base class for all *Array* classes.
 *
 *  The AllocatorBase class is the base class of all arrays
 *  stored in memory : it manages the main pointer on the data.
 *  It derive from the IContainerRef class as an array stored in
 *  memory can always be wrapped in some way or be a wrapper of
 *  data stored in memory.
 * 
 *  This class can also be used as a concrete class.
 *
 *  The class @c DATA can be any type of data that can be stored in memory.
 **/
template<typename DATA>
class AllocatorBase : public IContainerRef
{
  protected:
    /** Copy constructor. It is set as protected as we don't know, how the
     *  end-user want to copy the data.
     *  @param T : the array to copy
     *  @param ref : is this a wrapper of T ?
     **/
    AllocatorBase( const AllocatorBase& T, bool ref = false)
                 : IContainerRef(ref)
    {
      if (ref) // set T data
      {
        p_data_   = T.p_data_;
        rangeData_ = T.rangeData_;
      }
      else // Allocate memory in base class
      {
        setDefault();
      }
    }

  public:
    /** Default constructor
     *  @param I range of the data
     **/
    AllocatorBase( Range const& I = Range()) : IContainerRef(false)
                                          , p_data_(0)
                                          , rangeData_(I)
    { mallocPtrData( I.size(), I.first());}

    /** Wrapper constructor
     *  @param q ptr to the data to wrap
     *  @param I range of the data wrapped
     **/
    AllocatorBase(DATA* q, Range const& I) : IContainerRef(true)
                                       , p_data_(q)
                                       , rangeData_(I)
    { ;}

    /** Wrapper constructor : second form. This constructor assume that
     *  the data as a C-like array. Thus the first index is 0 and the last
     *  index is size-1.
     *  @param q ptr to the data to wrap
     *  @param size of the data to wrap
     **/
    AllocatorBase(DATA* q, Integer const& size) : IContainerRef(true)
                                            , p_data_(q)
                                            , rangeData_(Range(0,size-1))
    { ;}

    /** Virtual destructor. */
    ~AllocatorBase()
    { if (!this->isRef()) this->freePtrData(); }

    /** Get the first index of the data.
     * @return first index of the data
     **/
    inline Integer const& firstData() const { return rangeData_.first();}
   /** Get the last index of the data.
    * @return the last index of the data
    **/
    inline Integer const& lastData() const { return rangeData_.last();}
    /** Get the size of the data.
     *  @return the size of the data
     **/
    inline Integer const& sizeData() const { return rangeData_.size();}
    /** Get the range of the data.
     * @return the range of the data
     **/
    inline Range const& rangeData() const { return rangeData_;}
    /** Get the main ptr on the data.
     *  @return a pointer on the data set
     **/
    inline DATA* const& ptrData() const { return p_data_;}

    /** Get the const element number pos.
     *  @param pos the position of the element we get 
     **/
    inline DATA const& data( Integer const& pos) const
    {
#ifdef STK_BOUNDS_CHECK
      if (pos < firstData())
      { throw out_of_range(_T("AllocatorBase::data(pos) const "
                              "AllocatorBase::firstData() > pos"));
      }
      if (pos > lastData())
      { throw out_of_range(_T("AllocatorBase::data(pos) const "
                              "AllocatorBase::lastData() < pos"));
      }
#endif
      return p_data_[pos];
    }
    
    /** Get the element number pos.
     *  @param pos the position of the element we get 
     **/
    inline DATA& data(Integer const& pos)
    {
#ifdef STK_BOUNDS_CHECK
      if (pos < firstData())
      { throw out_of_range(_T("AllocatorBase::data(pos) "
                              "AllocatorBase::firstData() > pos"));
      }
      if (pos > lastData())
      { throw out_of_range(_T("AllocatorBase::data(pos) "
                              "AllocatorBase::lastData() < pos"));
      }
#endif
      return p_data_[pos];
    }
        
    /** @brief Set the element number pos.
     *  No check about the given position is performed in this method.
     *  @param pos the position of the element we set
     *  @param data the value we set 
     **/
   inline void setData( Integer const& pos, DATA const& data)
   {
#ifdef STK_BOUNDS_CHECK
      if (pos < firstData())
      { throw out_of_range(_T("AllocatorBase::setData(pos, data) "
                              "AllocatorBase::firstData() > pos"));
      }
      if (pos > lastData())
      { throw out_of_range(_T("AllocatorBase::setData(pos, data) "
                              "AllocatorBase::lastData() < pos"));
      }
#endif
     p_data_[pos] = data;
   }

    /** Return n+m, where n is the first number such that m < 2^n.
     *  if m <=0 return 0
     *  @param m the size of the container
     **/
    static Integer evalCapacity(Integer const& m)
    {
      Integer n = 0;
      Integer b = m;
      for (Integer k=1 ; k <= b; n++, k <<= 1);
      return m+n;
    }

  protected:
    /** Get the main ptr on the data. */
    inline DATA*& ptrData() { return p_data_;}

    /** swap this with T.
     *  @param T the container to swap
     **/
    void swap(AllocatorBase &T)
    {
      // copy main ptr data and range of this
      DATA*   auxMainPtr(p_data_);
      Range   auxRange(rangeData_);
      bool    auxRef(this->isRef());

      // overwrite main ptr data of this
      this->setPtrData(T.p_data_);
      this->setRangeData(T.rangeData_);
      this->setRef(T.isRef());

      // overwrite main PTRCOL of T
      T.setPtrData(auxMainPtr);
      T.setRangeData(auxRange);
      T.setRef(auxRef);
    }

    /** shift the first index of the data to first.
     *  @param first the index of the first data to set
     **/
    void shiftPtrData(Integer const& first)
    {
      // check if there is something to do
      if (first == firstData()) return;
      // check for reference
      if (this->isRef())
      { throw runtime_error(_T("In AllocatorBase::shiftPtrData(first)"
                               " can't operate on reference."));
      }
      // compute increment
      Integer inc = first - firstData();
      // translate data
      decPtrData(inc);
    }
    /** protected function for main ptr memory allocation.
     *  @param size the size to reserve in memory
     *  @param inc the increment to apply to p_data_ after allocation
     **/
    void mallocPtrData( Integer const& size, Integer const& inc  = 0)
    {
      // check for reference
      if (this->isRef())
      { throw runtime_error(_T("In AllocatorBase::mallocPtrData(first)"
                               " can't operate on reference."));
      }
      // delete any memory allocated
      freePtrData();
      // check size
      if (size <= 0)
      {
        // set initial values
        setRangeData(Range(0,size-1));
        // first index is 0
        return;
      }
      // allocate memory
      try
      {
        setPtrData(new DATA[size]);
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        // initialize to default
        setDefault();
        // and throw an Exception
        throw runtime_error(_T("AllocatorBase::mallocPtrData(size, inc) "
                               "memory allocation failed."));
      }
      // set initial values
      setRangeData(Range(0,size-1));  // first index is 0
      // apply increment
      decPtrData(inc);
    }

    /** @brief protected function for main ptr memory reallocation.
     *
     *  If the size requested is greater than the allocated size,
     *  the DATA stored are saved and copied using the operator=. the DATA
     *  class have to provide this operator.
     *
     *  If the size requested is lesser than the allocated size, only
     *  the first elements fitting in the container are copied.
     *
     *  @param size the size to reserve in memory
     *  @param inc the increment to apply to p_data_ after allocation
     **/
    void reallocPtrData( Integer const& size, Integer const& inc  = 0)
    {
      // check for reference
      if (this->isRef())
      { throw runtime_error(_T("In AllocatorBase::reallocPtrData(first)"
                               " can't operate on reference."));
      }
      // if there is no memory allocated, we use mallocPtrData
      if (!p_data_)
      {
        mallocPtrData(size, inc);
        return;
      }
      // if the new size is empty we can safely remove existing data
      if (size <= 0)
      {
        freePtrData();
        return;
      }
      // reset to zero increment and get data adress
      const Integer first = firstData();
      incPtrData(first);
      DATA* p_data = p_data_;
      // allocate memory
      try
      {
        p_data_ = new DATA[size];
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        // reset array to previous state
        p_data_ = p_data;
        // initialize members to default
        decPtrData(first);
        // and throw an Exception
        throw runtime_error(_T("AllocatorBase::reallocPtrData(size, inc) "
                                 "memory allocation failed."));
      }
      // copy existing data
      Integer nbData = min(size, sizeData());
      for (Integer i=0; i< nbData; i++)
      { p_data_[i] = p_data[i];}
      // delete existing data
      delete [] p_data;
      // set initial values
      setRangeData(Range(0,size-1));     // first index is 0
      // apply increment
      decPtrData(inc);
    }

    /** protected function for main ptr memory deallocation. */
    void freePtrData()
    {
      // nothing to do for ref
      if (this->isRef()) return;
      // if there is elts
      if (p_data_)
      {
        incPtrData(firstData());  // translate
        delete [] p_data_;   // erase
        setPtrData();        // set default values
      }
    }

    /** Set the address of the data : this method is not destined
     *  to the end-user.
     *  @param p_data the address to set
     *  @param rangeData the range of the data
     **/
    inline void setPtrData( DATA* p_data = 0, Range const& rangeData = Range())
    {
      p_data_    = p_data;
      rangeData_ = rangeData;
    }

  private:
    /** Main pointer on the data. */
    DATA* p_data_;
    /** Range of the data */
    Range rangeData_;

    /** Set the index of the first data : this method is not destined
     *  to the end-user.
     *  @param rangeData the range of the data to set
     **/
    inline void setRangeData( Range const& rangeData = Range())
    { rangeData_ = rangeData;}

    /** Increment the address of the data : this method is not destined
     *  to the end-user.
     *  @param inc the increment to apply
     **/
    inline void incPtrData( Integer const& inc)
    {
      // translate p_data_ only if there exists data
      if (p_data_) { p_data_ += inc;}
      rangeData_.dec(inc);
    }

    /** Decrement the address of the data : this method is not destined
     *  to the end-user.
     *  @param dec the increment to apply
     **/
    inline void decPtrData( Integer const& dec)
    {
      // translate p_data_ only if there exists data
      if (p_data_) { p_data_  -= dec;}
      rangeData_.inc(dec);
    }

    /** Set array members to default values. */
    inline void setDefault()
    {
      setPtrData();
      setRangeData();
    }
};

} // namespace STK

#endif /* STK_ARRAYBASE_H */
