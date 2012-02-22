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
 * Purpose:  Define the Interface for the Array classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ITArray2DBase.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ITARRAY2DBASE_H
#define STK_ITARRAY2DBASE_H

#include "../../Sdk/include/STK_ITContainer2D.h"

#include "STK_AllocatorBase.h"
#include "STK_Array1D.h"

namespace STK
{
/** @ingroup Arrays
  * @brief Templated interface base class for two-dimensional arrays.
  *
  * A ITArray2DBase is an interface class for two-dimensional containers
  * stored in columns.
  *
  * Each column has a range stored in the array @c rangeCols_ and a
  * capacity stored in the array @c capacityCols_. It should be worth
  * noting that we should have
  * @code
  *   (rangeCols_[j].size() <= capacityCols_[j]) == true;
  *   (rangeCols_[j].isIn(this->rangeVe()) == true;
  *@endcode
  *
  * Template parameter description:
  * - @c TYPE is the type of the data stored in the container.
  * - @c PTRCOL is the type of the ptr of column in a two-dimensional
  * array: for exemple @c TYPE*, @c Array1D<TYPE>*, @c DBACCESS*....
  * - @c TArrayHo is the name of the class that return a row of the array.
  * - @c TArrayVe is the name of the class that return a column of the array.
  * - @c TArray2D is the name of the class that implements @c ITArray2DBase.
  *
 **/
template < class TYPE
         , class PTRCOL
         , class TArray2D
         >
class ITArray2DBase : public ITContainer2D< TYPE, TArray2D>
                   , public AllocatorBase<PTRCOL>
{
  protected:
	/** capacity of the columns of the container (for each column: number of
   *  available Rows without reallocation in this column)
	 **/
	Array1D<Integer> capacityCols_;

	/** range of the Cols of the container.
	 **/
	Array1D<Range> rangeCols_;

  public:
    /** Type for the Base reference Class. */
    typedef AllocatorBase<PTRCOL> _AllocatorBaseType_;

    /** type of the Base Container Class. */
    typedef ITContainer2D<TYPE, TArray2D> _ITContainer2DType;

  protected:
    /** Default constructor
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    ITArray2DBase( Range const& I = Range(), Range const& J = Range())
                : _ITContainer2DType(I, J)
                , _AllocatorBaseType_()
                , capacityCols_()
                , rangeCols_()
                , capacityHo_(0)
    { mallocHo(J);}

    /** Copy constructor If we want to wrap T, the main ptr will be wrapped
     *  in AllocatorBase class. If we want to copy  T, AllocatorBase is
     *  initialized to default values.
     *  @param T the container to copy
     *  @param ref true if we wrap T
     **/
    ITArray2DBase( const ITArray2DBase& T, bool ref =false)
                : _ITContainer2DType(T)
                , _AllocatorBaseType_(T, ref)
                , capacityCols_(T.capacityCols_)
                , rangeCols_(T.rangeCols_)
                , capacityHo_(T.capacityHo_)
    {
      if (!ref)
        mallocHo(T.rangeHo());
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I range of the Rows to wrap
     *  @param J range of the column to wrap
     **/
    ITArray2DBase( const ITArray2DBase& T, Range const& I, Range const& J)
                : _ITContainer2DType(I, J)
                , _AllocatorBaseType_(T, true)
                , capacityCols_(J)
                , rangeCols_(J)
                , capacityHo_(T.capacityHo_)
    {
#ifdef STK_BOUNDS_CHECK
      if (I.first() < T.firstRow())
      { throw out_of_range("ITArray2DBase::ITArray2DBase(T, I, J) "
                           "I.first() < T.firstRow()");
      }
      if (I.last() > T.lastRow())
      { throw out_of_range("ITArray2DBase::ITArray2DBase(T, I, J) "
                           "I.last() > T.lastRow()");
      }
      if (J.first() < T.firstCol())
      { throw out_of_range("ITArray2DBase::ITArray2DBase(T, I, J) "
                           "J.first() < T.firstCol()");
      }
      if (J.last() > T.lastCol())
      { throw out_of_range("ITArray2DBase::ITArray2DBase(T, I, J) "
                           "J.last() > T.lastCol()");
      }
#endif
      // adjust capacity and range of each Cols
      for (Integer j=J.first(); j<=J.last(); j++)
      {
        // copy capacity of the column j (is it necessary ?)
        capacityCols_[j] = T.capacityCols_[j];
        // compute available range of the column j
        rangeCols_[j] = Range::inf(I, T.rangeCols_[j]);
      }
    }

    /** Wrapper constructor We get a reference of the data.
     *  @param q pointer on data
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    ITArray2DBase( PTRCOL* q, Range const& I, Range const& J)
                : _ITContainer2DType(I, J)
                , _AllocatorBaseType_(q, J)
                , capacityCols_(J, I.size())
                , rangeCols_(J, I)
                , capacityHo_(0)
    { ;}

  public:
    /** Virtual destructor. Allocated horizontal memory is liberated by
     * the base class AllocatorBase.
     **/
    virtual ~ITArray2DBase()
    { ;}

    /** give the maximum possible number of Cols without
     *  reallocation.
     **/
    inline Integer const& capacityHo() const
    { return capacityHo_;}

    /** give the maximum possible number of Rows without
     *  reallocation for all Cols.
     **/
    inline const Array1D<Integer> & capacityCols() const
    { return capacityCols_;}

    /** give the capacity of a Col.
     *  @param col index of the column we want the range
     **/
    inline Integer const& capacityCol(Integer const& col) const
    { return capacityCols_[col];}

    /** give the range of all Cols.
     **/
    inline Array1D<Range> const& rangeCols() const
    { return rangeCols_;}

    /** give the range of a Col.
     *  @param col index of the column we want the range
     **/
    inline Range const& rangeCol(Integer const& col) const
    { return rangeCols_[col];}

    /** New beginning index for the Cols of the object.
     *  @param cbeg the index of the first column to set
     **/
    void shiftHo(Integer const& cbeg =1)
    {
      // compute increment
      Integer cinc = cbeg - this->firstCol();
      // if there is something to do
      if (cinc != 0)
      {
        // is this structure just a pointer?
        if (this->isRef())
        { throw runtime_error("ITArray2DBase::shiftHo(cbeg) "
                                   "can't operate on references.");
        }
        // translate rangeHo_
        this->incRangeHo(cinc);
        // translate data
        this->shiftPtrData(cbeg);
        // tranlate capacityCols_
        capacityCols_.shift(cbeg);
        // translate rangeCols_
        rangeCols_.shift(cbeg);
      }
    }

    /** Swapping the pos1 column and the pos2 column.
     *  @param pos1 position of the first col
     *  @param pos2 position of the second col
     **/
    void swapCols(Integer const& pos1, Integer const& pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->firstCol() > pos1)
      { throw out_of_range("ITArray2D::swapCols(pos1, pos2) "
                           "this->firstCol() >pos1");
      }
      if (this->lastCol() < pos1)
      { throw out_of_range("ITArray2D::swapCols(pos1, pos2) "
                           "this->lastCol() <pos1");
      }
      if (this->firstCol() > pos2)
      { throw out_of_range("ITArray2D::swapCols(pos1, pos2) "
                           "this->firstCol() >pos2");
      }
      if (this->lastCol() < pos2)
      { throw out_of_range("ITArray2D::swapCols(pos1, pos2) "
                           "this->lastCol() <pos2");
      }
#endif
      // swap
      PTRCOL qaux(this->data(pos1));
      this->setData(pos1, this->data(pos2));
      this->setData(pos2, qaux);
      // update capacityCols_
      STK::swap(capacityCols_[pos1],  capacityCols_[pos2]);
      // update rangeCols_
      STK::swap(rangeCols_[pos1],  rangeCols_[pos2]);
    }

    /** swap this container with T.
     * @param T the container to swap with this
     **/
    void swap(ITArray2DBase &T)
    {
      // swap AllocatorBase part
      this->AllocatorBase<PTRCOL>::swap(T);

      // swap ITContainer2D part
      this->_ITContainer2DType::swap(T);

      // swap ITArray2DBase part
      STK::swap(this->capacityHo_, T.capacityHo_);
      capacityCols_.swap(T.capacityCols_);
      rangeCols_.swap(T.rangeCols_);
    }

    /** Append the container @c T to @c this without copying the data
     *  explicitly. The column of @c T are appended to this and
     *  @c T will become a reference container. Observe that the @c const
     *  keyword is not respected in this method: but it is useful to
     *  define this method even for constant objects. The data in itself are not
     *  altered, the Array2D become a reference on its own data.
     *  @param T the container to append to this
     **/
    void pushBackByTransfer(ITArray2DBase const& T)
    {
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("In ITArray2DBase::pushBackByTransfer(T) "
                                 "(*this) is a reference.");
      }
      // is T just a pointer?
      if (T.isRef())
      { throw runtime_error("In ITArray2DBase::pushBackByTransfer(T) "
                                 "T is a reference.");
      }
      // if there is no columns, we can safely modify the vertical range
      if (this->sizeHo() <= 0)
        this->setRangeVe(T.rangeVe());
      // Are ranges corrects ?
      if (this->rangeVe() != T.rangeVe())
      { throw runtime_error("In ITArray2DBase::pushBackByTransfer(T) "
                                 "this->rangeVe() != T.rangeVe().");
      }
      // break const reference
      ITArray2DBase& Tref = const_cast<ITArray2DBase&>(T);
      // compute horizontal range of the container after insertion
      Range rangeHo(this->rangeHo());
      // compute first index of the first column added
      const Integer first = rangeHo.last() + 1;
      // reallocate memory
      rangeHo.incLast(Tref.sizeHo());
      reallocHo(rangeHo);
      this->setRangeHo(rangeHo);
      // align T range
      const Integer last = rangeHo.last();
      Tref.shiftHo(first);
      // copy data from T
      for (Integer j=first; j<= last; j++)
      {
        copyColumn(Tref, j, j);
      }
      // release memory allocated for the columns
      Tref.freePtrData();
      Tref.setPtrData(this->ptrData(), this->rangeData());
      Tref.setRef(true);
    }

  protected:
    /** set the maximum possible number of Cols without
     *  reallocation.
     *  @param capacity the maximum number of columns
     **/
    inline void setCapacityHo(Integer const& capacity = 0)
    { capacityHo_ = capacity;}

    /** copy the column pos2 of the container T to the columnthe column of
     *  pos1 of this. On of the container (either this or T but not both)
     *  have to be a reference otherwise, user will experiment a memory leak.
     *
     *  @param T the container with the column to transfer
     *  @param pos1 index of the column to initialize
     *  @param pos2 the column in the container T to transfer in this
     **/
    void copyColumn( ITArray2DBase const& T, Integer const& pos1, Integer const& pos2)
    {
      // copy column pos2 of T in pos1 of this
      this->setData(pos1, T.data(pos2));
      // set capacityCols_
      capacityCols_[pos1] = T.capacityCols_[pos2];
      // set rangeCols_
      rangeCols_[pos1] = T.rangeCols_[pos2];
    }

    /** Transfer the column pos2 of the container T to the column
     *  pos1 of this. Set the column pos2 in T to a default value.
     *  The column pos1 should not exists or should be deleted
     *  in this otherwise user will experiment a memory leak.
     *
     *  @param T the container with the column to transfer
     *  @param pos1 index of the column to initialize
     *  @param pos2 the column in the container T to transfer in this
     **/
    void transferColumn( ITArray2DBase& T, Integer const& pos1, Integer const& pos2)
    {
      // copy column pos2 of T in pos1 of this
      this->setData(pos1, T.data(pos2));
      // set capacityCols_
      capacityCols_[pos1] = T.capacityCols_[pos2];
      // set rangeCols_
      rangeCols_[pos1] = T.rangeCols_[pos2];
      // set column of T to default
      T.setDefaultCol(pos2);
    }

    /** Method for memory allocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the rangeHo_ of the container to default.
     *  @param J horizontal range
     **/
    void mallocHo(Range const& J)
    {
      // compute the size necessary (can be 0)
      Integer size = this->evalCapacity(J.size());
      // try to allocate memory
      try
      {
        // initialize this->capacityCols_
        capacityCols_.resize(J);
        // initialize this->rangeCols_
        rangeCols_.resize(J);
        // allocate memory for the columns
        this->mallocPtrData(size, J.first());
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0)
        setCapacityHo();
        // set default range
        this->setRangeHo();
        // clear this->capacityCols_
        capacityCols_.clear();
        // clear this->rangeCols_
        rangeCols_.clear();
        // throw the error
        throw error;
      }
      // set new capacity if no error occur
      this->setCapacityHo(size);
    }

    /** Method for memory reallocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the rangeHo_ of the container to default.
     *  @param J horizontal range
     **/
    void reallocHo(Range const& J)
    {
      // compute the size necessary (can be 0)
      Integer size = this->evalCapacity(J.size());
      // try to allocate memory
      try
      {
        // allocate memory for the columns
        this->reallocPtrData(size, J.first());
        // initialize this->capacityCols_
        capacityCols_.resize(J);
        // initialize this->rangeCols_
        rangeCols_.resize(J);
      }
      catch (runtime_error & error)   // if an error occur
      {
        // set default capacity (0)
        this->setCapacityHo();
        // set default range
        this->setRangeHo();
        // clear this->capacityCols_
        this->capacityCols_.clear();
        // clear this->rangeCols_
        this->rangeCols_.clear();
        // throw the error
        throw error;
      }
      // set new capacity if no error occur
      this->setCapacityHo(size);
    }

    /** Horizontal Memory deallocation.
     *  This method clear all allocated memory. The range of the columns
     *  is set to (firstCol_:firstCol_-1). The range of the Rows remain
     *  unmodified. If there is allocated memory for the columns, it
     *  should be liberated prior to this method.
     **/
    void freeHo()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free memory allocated in AllocatorBase
      this->freePtrData();
      // set capacity size to default
      this->setCapacityHo();
      // set range of the Cols to default
      this->setRangeHo(Range(this->firstCol(), this->firstCol()-1));
      // set capacityCols_ to default
      capacityCols_.clear();
      // set rangeCols_ to default
      rangeCols_.clear();
    }

  private:
    /** Horizontal capacity of the container (number of available
     *  Cols without reallocation)
     **/
    Integer capacityHo_;

    /** set the default parameters and dimension to a column of
     *  the container.
     *  @param pos the position of the column to initialize to a default
     *  value.
     **/
    inline void setDefaultCol(Integer const& pos)
    {
      // set column of T to default
      this->setData(pos, (PTRCOL)NULL);
      // set capacityCols_
      this->capacityCols_[pos] = 0;
      // set rangeCols_
      this->rangeCols_[pos] = Range();
    }
};

} // namespace STK

#endif
// STK_ITARRAY2DBASE_H
