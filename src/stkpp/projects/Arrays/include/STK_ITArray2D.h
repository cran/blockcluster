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
 * Project:  stkpp::Arrrays
 * Purpose:  Define the Interface for the Array classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ITArray2D.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ITARRAY2D_H
#define STK_ITARRAY2D_H

#include "../../STKernel/include/STK_Range.h"
#include "STK_ITArray2DBase.h"

#include "STK_Array1D.h"
#include "STK_ArrayHo.h"

namespace STK
{

/** @ingroup Arrays
  * @brief Templated interface base class for two-dimensional arrays.
  *
  * A ITArray2D is a specialized interface class for two-dimensional
  * containers stored in columns. All derived class from @c ITArray2D
  * access to the column using a @c TYPE* ptr.
  *
  * Template parameter description:
  * - The template parameter @c TYPE is the fundamental type of the
  *  data stored in the container.
  * - The template parameter @c TArray2D is the name of the class
  *  implementing @c ITArray2D.
  *
  * The derived class have to implement the following public method
  * @code
  *    Range compRangeVe( Integer const& icol) const
  * @endcode
 **/
template < class TYPE, class  TArray2D  >
class ITArray2D : public ITArray2DBase< TYPE, TYPE*, TArray2D>
{
    /** type of the Base Container Class. */
    typedef ITArray2DBase< TYPE, TYPE*, TArray2D> _AllocatorBaseType_;

  protected:
    /** Default constructor
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    ITArray2D( Range const& I = Range(), Range const& J = Range())
            : _AllocatorBaseType_(I, J)
    { ;}

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if we wrap T
     **/
    ITArray2D( const ITArray2D& T, bool ref =false)
            : _AllocatorBaseType_(T, ref)
    { ;}

    /** constructor by reference, ref_=1.
     *  @param T the container to copy
     *  @param I range of the Rows to wrap
     *  @param J range of the Col to wrap
     **/
    ITArray2D( const ITArray2D& T, Range const& I, Range const& J)
            : _AllocatorBaseType_(T, I, J)
    { ;}

    /** Wrapper constructor The Container is a ref.
     *  @param q pointer on data
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    ITArray2D( TYPE** q, Range const& I, Range const& J)
            : _AllocatorBaseType_(q, I, J)
    { ;}

  public:
    /** Virtual destructor.
     *  free the vertically allocated memory (the columns). The horizontally
     *  allocated memory is handled by the AllocatorBase class.
     **/
    virtual ~ITArray2D()
    { if (!this->isRef())
        this->freeCols(this->rangeHo());
    }

    /** access to one element.
     *  @param i index of the row
     *  @param j index of the col
     *  @return a reference on the (i,j) element
     **/
    inline TYPE& elt(Integer const& i, Integer const& j)
    { return this->data(j)[i];}

    /** access to one element const.
     *  @param i index of the row
     *  @param j index of the col
     *  @return a constant reference on the (i,j) element
     **/
    inline TYPE const& elt(Integer const& i, Integer const& j) const
    { return this->data(j)[i];}

    /** access to a part of a column.
     *  @param I range of the rows
     *  @param j index of the column
     *  @return A reference with range I on the column j of this
     **/
    inline Array1D<TYPE> col(Range const& I, Integer const& j) const
    { return Array1D<TYPE>(this->data(j), I, j);}

    /** access to a column.
     *  @param j index of the column
     *  @return A reference on the column j of this
     **/
    inline Array1D<TYPE> col( Integer const& j) const
    { return Array1D<TYPE>(this->data(j), this->rangeVe(), j);}

    /** access to a part of a row.
     *  @param i index of the row
     *  @param J range of the columns
     *  @return A reference with range J on the row i
     **/
    inline ArrayHo<TYPE> row( Integer const& i, Range const& J) const
    { return ArrayHo<TYPE>(*this, J, i);}

    /** access to a row.
     *  @param i index of the row
     *  @return A reference on the row i
     **/
    inline ArrayHo<TYPE> row( Integer const& i) const
    { return ArrayHo<TYPE>(*this, this->rangeHo(), i);}

    /** clear the object.
     *  This will free all allocated memory and reset all range
     *  to 1:0.
     **/
    void clear()
    {
      // free allocated mem
      this->freeMem();
      // Set dimensions to default
      this->setRange();
    }

    /** New first indexes for the object.
     *  @param rbeg the index of the first row to set
     *  @param cbeg the index of the first column to set
     **/
    void shift(Integer const& rbeg =1, Integer const& cbeg =1)
    {
      // move begin of the col
      this->shiftHo(cbeg);
      // move begin of the row
      this->shiftVe(rbeg);
    }

    /** New beginning index for the Rows of the object.
     *  @param rbeg the index of the first row to set
     **/
    void shiftVe(Integer const& rbeg =1)
    {
       // compute increment
      Integer rinc = rbeg - this->firstRow();
      // if there is something to do
      if (rinc != 0)
      {
        // is this structure just a pointer?
        if (this->isRef())
        { throw runtime_error("ITArray2D::shiftVe(rbeg) "
                              "can't operate on references.");
        }
        // translate rangeVe_()
        this->incRangeVe(rinc);
        // For all cols, move begin
        for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
        {
          shiftCol(j, this->rangeCols_[j].first()+rinc);
        }
      }
    }

    /** Add n Rows to the container.
     *  @param n number of Rows to add
     **/
    void pushBackRows(Integer const& n=1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::pushBackRows(n) "
                            "can't operate on references.");
      }
      // If the container have no rows : create its
      if (this->sizeVe() <=0)
      {
        // update the range of the container
        this->incLastVe(n);
        // initialize the container
        this->initializeCols(this->rangeHo());
      }
      else
      {
        // update the range of the rows
        this->incLastVe(n);
        // allocate new Rows for each Col
        for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
        {
          // compute range from the leaf
          Range range(this->asLeaf().compRangeVe(j));
          // if there is no column or the end is less than the container
          // end
          if ((range.size()>0)&&(range.last()>this->lastRow()-n))
          {
            // if the column is empty create it
            if (this->rangeCols_[j].size()<=0)
            {
              this->initializeCol(j, range);
            }
            else
            {
              // compute position
              Integer pos(this->lastRow()-n+1);
              // add elts
              insertRowsToCol(j, pos, range.last() - pos +1);
            }
          }
        }
      }
    }

    /** Insert n Rows at the position pos of the container.
     *  If pos is outside the range of a column, then the
     *  method do nothing.
     *  @param pos index where to insert Rows
     *  @param n number of elements to insert (default 1)
     **/
    void insertRows(Integer const& pos, Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
#ifdef STK_BOUNDS_CHECK
      // check indices
      if (this->firstRow() > pos)
      { throw out_of_range("ITArray2D::insertRows(pos, n) "
                           "this->firstRow() > pos");
      }
      if (this->lastRow()+1 < pos)
      { throw out_of_range("ITArray2D::insertRows(pos, n) "
                           "this->lastRow()+1 < pos");
      }
#endif
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::insertRows(pos, n) "
                            "can't operate on references.");
      }
      // update the range of the rows
      this->incLastVe(n);
      // allocate new Rows for each Col
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
      {
        // check position
        if ( (pos >= this->rangeCols_[j].first())
           ||(pos <= this->rangeCols_[j].last()+1)
           )
        {
          insertRowsToCol(j, pos, n);
        }
      }
    }

    /** Delete last Rows of the container.
     *  @param n number of Rows to delete
     **/
    void popBackRows(Integer const& n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::popBackRows() "
                            "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // if there is Rows to erase
      if (this->sizeVe()<n)
      { throw out_of_range("ITArray2D::popBackRows(n) "
                           "this->sizeVe() < n");
      }
#endif
      // update range of the container
      this->decLastVe(n);
      // decrease range of each Col
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
        eraseRowsToCol(j, this->lastRow()+1, n);
    }

    /** Delete n Rows at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void eraseRows(Integer const& pos, Integer const& n=1)
    {
      // if n==0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::eraseRows(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstRow() > pos)
      { throw out_of_range("ITArray2D::eraseRows(pos, n) "
                           "this->firstRow() > pos");
      }
      if (this->lastRow() < pos)
      { throw out_of_range("ITArray2D::eraseRows(pos, n)"
                           " this->lastRow() < pos");
      }
      if (this->lastRow() < pos+n-1)
      { throw out_of_range("ITArray2D::eraseRows(pos, n)"
                           " this->lastRow() < pos+n-1");
      }
#endif
      // save posistion and size
      Integer posRow = pos, nbRow = n;
      // update dimensions
      this->decLastVe(n);
      // update each Col
      for (Integer j=this->firstCol(); j<=this->lastCol(); j++)
        eraseRowsToCol(j, posRow, nbRow);
    }

    /** Add n Cols to the container.
     *  @param n the number of Cols to add
     **/
    void pushBackCols(Integer const& n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::pushBackCols(n) "
                            "can't operate on references.");
      }
      // If the container have no Cols : create its
      if (this->sizeHo() <=0)
      {
        // update end col
        this->incLastHo(n);
        // initialize Horizontally the container
        this->mallocHo(this->rangeHo());
        // initialize Vertically the container
        initializeCols( this->rangeHo());
      }
      else // else insert to the end of the container
        insertCols(this->lastCol()+1, n);
    }

    /** Insert n Columns at the index pos to the container.
     *  @param pos the position of the inserted Cols
     *  @param n the number of column to insert
     **/
    void insertCols(Integer const& pos, Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::insertCols(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > pos)
      { throw out_of_range("ITArray2D::insertCols(pos, n) "
                                "this->firstCol() > pos");
      }
      if (this->lastCol()+1 < pos)
      { throw out_of_range("ITArray2D::insertCols(pos, n) "
                                "this->lastCol()+1 < pos");
      }
#endif
      // compute horizontal range of the container after insertion
      Range range_ho(this->rangeHo());
      range_ho.incLast(n);
      // allocate, if necessary, the mem for the Cols
      if (this->capacityHo() < range_ho.size()) //  not enough space
      {
        // temporary empty container (the number of Rows is the same)
        // but there is no Cols
        TArray2D Taux(this->rangeVe(), Range());
        // swap with Taux
        this->swap(Taux);
        // initialize columns of the container
        try
        {
          this->mallocHo(range_ho);
        }
        catch (runtime_error & error)   // if an error occur
        {
          this->swap(Taux);   // restore container
          throw error;        // and send again the stdexcept
        }
        // set the range of the Columns
        this->setRangeHo(range_ho);
        // move first Columns from Taux to this
        for (Integer k=this->firstCol(); k<pos; k++)
          this->transferColumn(Taux, k, k);
        // translate and copy last Cols from Taux to this
        for (Integer k=Taux.lastCol(); k>=pos; k--)
          this->transferColumn(Taux, k+n, k);
      }
      else // enough space -> shift the last Cols
      {
        Range addRange(this->lastCol()+1, this->lastCol()+n);
        // insert capacity for the new Cols
        this->capacityCols_.insert(addRange, 0);
        // insert range for the new Cols
        this->rangeCols_.insert(addRange, Range());
        // update range_
        this->incLastHo(n);
        // translate data
        for (Integer k=this->lastCol()-n; k>=pos; k--)
          this->transferColumn( this->asLeaf(), k+n, k);
      }
      // initialize the rows for the Cols, this->capacityCols_, this->rangeCols_
      // in the range pos:pos+n-1
      this->initializeCols(Range(pos, pos+n-1));
    }

    /** Delete last Cols of the container
     *  @param n the number of Cols to delete
     **/
    void popBackCols(Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::popBackCols(n) "
                            "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check range
      if (this->sizeHo() < n)
      { throw out_of_range("ITArray2D::popBackCols(n) "
                           "this->sizeHo() < n");
      }
#endif
      // delete each col
      this->freeCols(Range(this->lastCol()-n+1, this->lastCol()));
      // update this->capacityCols_
      this->capacityCols_.popBack(n);
      // update this->rangeCols_
      this->rangeCols_.popBack(n);
      // update rangeHo
      this->decLastHo(n);
      // if there is no more Cols
      if (this->sizeHo() == 0) this->freeMem();
    }

    /** Delete n Cols at the specified position of the container.
     *  @param pos the position of the deleted Cols
     *  @param n the number of column to delete
     **/
    void eraseCols(Integer const& pos, Integer const& n = 1)
    {
      if (n<=0) return;        // if n<=0 nothing to do
#ifdef STK_BOUNDS_CHECK
      // check range
      if (this->firstCol() > pos)
      { throw out_of_range("ITArray2D::eraseCols(pos, n) "
                           "this->firstCol() > pos");
      }
      if (this->lastCol() < pos)
      { throw out_of_range("ITArray2D::eraseCols(pos, n) "
                           "this->lastCol() < pos");
      }
      if (this->lastCol() < pos+n-1)
      { throw out_of_range("ITArray2D::eraseCols(pos, n) "
                           "this->lastCol() < pos+n-1");
      }
#endif
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("ITArray2D::eraseCols(pos, n) "
                            "can't operate on references.");
      }
      // delete each col
      this->freeCols(Range(pos, pos+n-1));
      // update rangeHo_
      this->decLastHo(n);
      // shift Cols
      for (Integer k=pos; k<=this->lastCol(); k++)
        this->setData(k, this->data(k+n));
      // update this->capacityCols_
      this->capacityCols_.erase(pos, n);
      // update this->rangeCols_
      this->rangeCols_.erase(pos, n);
      // if there is no more Cols
      if (this->sizeHo() == 0) this->freeMem();
    }

    /** Update the cols of the container in the specified range.
     *  @param J range of the column to udpate
     **/
    void update(Range const& J)
    {
#ifdef STK_BOUNDS_CHECK
      // check range
      if (this->firstCol() > J.first())
      { throw out_of_range("ITArray2D::update(J) "
                                "this->firstCol() > J.first()");
      }
      if (this->lastCol() < J.last())
      { throw out_of_range("ITArray2D::update(J) "
                                "this->lastCol() < J.last()");
      }
#endif
      Integer firstCol(J.first());
      Integer lastCol(J.last());

      for ( Integer icol = firstCol; icol <= lastCol ; ++icol)
      {
         if (this->asLeaf().compRangeVe(icol) != this->rangeCol(icol))
         { this->resizeCol(icol, this->asLeaf().compRangeVe(icol));}
      }
    }

    /** Update the cols of the container in the specified position.
     *  @param pos index of the column to update
     **/
    void update(Integer const& pos)
    {
#ifdef STK_BOUNDS_CHECK
      // check range
      if (this->firstCol() > pos)
      { throw out_of_range("ITArray2D::update(pos) "
                           "this->firstCol() > pos");
      }
      if (this->lastCol() < pos)
      { throw out_of_range("ITArray2D::update(pos) "
                           "this->lastCol() < pos");
      }
#endif
      if (this->asLeaf().compRangeVe(pos) != this->rangeCol(pos))
      { this->resizeCol(pos, this->asLeaf().compRangeVe(pos));}
    }

  protected:
    /** Memory deallocation.
     *  This method clear all allocated memory. The range of the Cols
     *  is set to (beginHo_:beginHo_-1). The range of the Rows remain
     *  unmodified.
     **/
    void freeMem()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free the Rows memory
      this->freeCols(this->rangeHo());
      // liberate horizontally
      this->freeHo();
    }

    /** Function for memory allocation and initialization.
     *  This method will free all allocated memory owned by this
     *  container before initialization. If you don't want to free
     *  the allocatd memory, you should use @c resize(I, J)
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    void initialize(Range const& I, Range const& J)
    {
      // check if there is memory allocated
      this->clear();
      // if we initialize the memory the container is not a reference
      this->setRef(false);
      // create this->capacityCols_
      this->capacityCols_.resize(J);
      // create this->rangeCols_
      this->rangeCols_.resize(J);
      // set the Horizontal range of the container
      this->setRangeHo(J);
      // initialize Horizontally the container
      this->mallocHo(J);
      // set the vertical range of the container
      this->setRangeVe(I);
      // initialize Vertically the container
      this->initializeCols(J);
    }

    /** Function for memory allocation and initialization.
     *  The capacity for the Rows have to be set before calling this
     *  method.
     *  @param J vertical range of the Cols to initialize
     **/
    void initializeCols(Range const& J)
    {
    // for each col
      for (Integer j=J.first(); j<=J.last(); j++)
      {
        // try to Allocate mem for the jth col
        try
        {
          // initialize the cols with the computed range
          // specific to the container
          this->initializeCol(j, this->asLeaf().compRangeVe(j));
        }
        catch (runtime_error & error)   // if an error occur
        {
          // free each column allocated
          for (Integer k=J.first(); k<j; k++)
            this->freeCol(k);
          // put default for the other Cols
          for (Integer k=j; k<=J.last(); k++) this->setData(k, 0);
          // and throw an stdexcept
          throw error;
        }
      }
    }

    /** @brief internal method for initializing a column.
     *
     *  Method for the the allocation of memory of the col
     *  pos with the given range.
     *  @param pos the index of the column to initialize
     *  @param I   range of the Col
     **/
    void initializeCol(Integer const& pos, Range const& I)
    {
      if (I.size() <=0)
      {
        // set default for ptr
        this->setData(pos, (TYPE*)NULL);
        // set default value for this->capacityCols_[pos]
        this->capacityCols_[pos] = 0;
        // set default value for this->rangeCols_[pos]
        this->rangeCols_[pos] = I;
        // return
        return;
      }
      // compute the size necessary (cannot be 0)
      Integer size = this->evalCapacity(I.size());
      // try to allocate memory
      try
      {
        this->setData(pos, new TYPE[size]);
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        // set default for ptr
        this->setData(pos, (TYPE*)NULL);
        // set default value for this->capacityCols_[pos]
        this->capacityCols_[pos] = 0;
        // set default value for this->rangeCols_[pos]
        this->rangeCols_[pos] = Range();
        // and throw an stdexcept
        throw runtime_error("ITArray2D::initializeCol(pos, J) "
                                 "memory allocation failed.");
      }
      // increment ptr of the column
      this->data(pos) -= I.first();
      // set size for this->capacityCols_[pos]
      this->capacityCols_[pos] = size;
      // set value for this->rangeCols_[pos]
      this->rangeCols_[pos] = I;
    }

    /** vertical memory deallocation.
     *  @param J the range of the Cols to liberate.
     **/
    void freeCols(Range const& J)
    {
      // for all Cols
      for (Integer j=J.first(); j<=J.last(); j++)
        this->freeCol(j);
    }

    /** Method for memory deallocation.
     *  @param pos the number of the column to free
     **/
    void freeCol(Integer const& pos)
    {
      if (this->data(pos)) // if there is cols
      {
        // increment the ptr
        this->data(pos) += this->rangeCols_[pos].first();
        // delete allocated mem for the column pos
        delete [] this->data(pos);
        // set default value for ptr
        this->setData(pos, (TYPE*)NULL);
        // set default value for this->capacityCols_[pos]
        this->capacityCols_[pos] = 0;
        // set default value for this->rangeCols_[pos]
        this->rangeCols_[pos] = Range();
      }
    }

    /** @brief internal method for translating a column.
     *
     *  Method for the the allocation of memory of the col
     *  pos with the given range.
     *  @param pos the index of the column to translate
     *  @param beg new begin ofthe col
     **/
    void shiftCol( Integer const& pos, Integer const& beg)
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > pos)
      { throw out_of_range("ITArray2D::shiftCol(pos, n) "
                           "this->firstCol() > pos");
      }
      if (this->lastCol() < pos)
      { throw out_of_range("ITArray2D::shiftCol(pos, n) "
                           "this->lastCol() < pos");
      }
#endif
      // compute increment
      Integer rinc = beg - this->rangeCols_[pos].first();
      // check if there is data
      if (this->data(pos))
      {
        // transate ptr
        this->data(pos) -= rinc;
      }
      // translate this->rangeCols_
      this->rangeCols_[pos].inc(rinc);
    }

    /** @brief Internal method for resizing a column with a specified
     *  range.
     *
     *  This method resize the column @c pos to the desired range
     *  using:
     * - @c shiftCol
     * - either @c popBackRowsToCol or @c pushBackRowsToCol if needed.
     *  @param pos index of the column
     *  @param I range to set to the column
    **/
    void resizeCol( Integer const& pos, Range const& I)
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > pos)
      { throw out_of_range("ITArray2D::resizeCol(pos, I) "
                           "this->firstCol() > pos");
      }
      if (this->lastCol() < pos)
      { throw out_of_range("ITArray2D::resizeCol(pos, I) "
                           "this->lastCol() < pos");
      }
#endif
      // check if there is something to do
      if (this->rangeCol(pos) == I) return;
      // shift to the desired first index
      shiftCol(pos, I.first());
      // compute difference of size
      Integer inc = this->rangeCol(pos).size() - I.size();
      // nothing to do
      if (inc == 0) return;
      // add row
      if (inc < 0)
      {
        pushBackRowsToCol(pos, -inc);
      }
      else // delete rows
      {
        popBackRowsToCol(pos, inc);
      }
    }

    /** @brief Internal method for inserting rows to a specified column.
     *
     *  Insert n Rows at the position pos to the column column of the
     *  container. No check is done about the index.
     *  @param col column index
     *  @param pos index where to insert Rows
     *  @param n number of elements to insert (default 1)
     **/
    void insertRowsToCol( Integer const& col
                        , Integer const& pos
                        , Integer const& n =1
                        )
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > col)
      { throw out_of_range("ITArray2D::insertRowsToCol(col, n) "
                           "this->firstCol() > col");
      }
      if (this->lastCol() < col)
      { throw out_of_range("ITArray2D::insertRowsToCol(col, n) "
                           "this->lastCol() < col");
      }
      if (this->rangeCols_[col].first() > pos)
      { throw out_of_range("ITArray2D::insertRowsToCol(col, pos, n) "
                           "this->rangeCols_[col].first() > pos");
      }
      if (this->rangeCols_[col].last()+1 < pos)
      { throw out_of_range("ITArray2D::insertRowsToCol(col, pos, n) "
                           "this->rangeCols_[col].last()+1 < pos");
      }
#endif
      // wrap old Col
      TYPE* ptr_old_col(this->data(col));
      // get vertical range of the Col
      Range range_ve(this->rangeCols_[col]);
      // update range
      this->rangeCols_[col].incLast(n);
      // allocate if necessary the Col
      if (this->capacityCols_[col] < this->rangeCols_[col].size())
      {
        // create new Col
        this->initializeCol(col, this->rangeCols_[col]);
        // if there was data, copy and liberate
        if (ptr_old_col)
        {
          // get ptr on the new col
          TYPE* ptr_new_col(this->data(col));
          // copy first Elts
          for (Integer k=range_ve.first(); k<pos; k++)
            ptr_new_col[k] = ptr_old_col[k];
          // translate and copy last Elts
          for (Integer k=range_ve.last(); k>=pos; k--)
            ptr_new_col[k+n] = ptr_old_col[k];
          // increment ptr_col
          ptr_old_col += range_ve.first();
          // and free old col
          delete [] ptr_old_col;
        }
      }
      else // enough space
      {
        // translate last Elts
        for (Integer k=range_ve.last(); k>=pos; k--)
          ptr_old_col[k+n] = ptr_old_col[k];
      }
    }

    /** @brief Internal method for appending rows to a specified column.
     *
     *  Push back n Rows at the end of the column column of the
     *  container.
     *  @param col column index
     *  @param n number of elements to append (default 1)
     **/
    void pushBackRowsToCol( Integer const& col
                          , Integer const& n =1
                          )
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > col)
      { throw out_of_range("ITArray2D::pushBackRowsToCol(col, n) "
                           "this->firstCol() > col");
      }
      if (this->lastCol() < col)
      { throw out_of_range("ITArray2D::pushBackRowsToCol(col, n) "
                           "this->lastCol() < col");
      }
#endif
      // wrap old Col
      TYPE* ptr_old_col(this->data(col));
      // get vertical range of the Col
      Range range_ve(this->rangeCol(col));
      // compute vertical range of the Col after insertion
      //range_ve.incLast(n);
      this->rangeCols_[col].incLast(n);
      // allocate if necessary the Col
      if (this->capacityCols_[col] < range_ve.size())
      {
        // create new Col
        this->initializeCol(col, this->rangeCols_[col]);
        // ger ptr on the new col
        TYPE* ptr_new_col(this->data(col));
        // copy first Elts
        for (Integer k=range_ve.first(); k<=range_ve.last(); k++)
          ptr_new_col[k] = ptr_old_col[k];
        // if there is data
        if (ptr_old_col)
        {
          // increment ptr_col
          ptr_old_col += range_ve.first();
          // and free old col
          delete [] ptr_old_col;
        }
      }
    }

    /** @brief Internal method for deleting rows to a specified column.
     *
     *  Delete n Rows at the position pos to the column col of the container.
     *  No check is done about indexes. It is possible to remove data
     *  outside the range of the column. In this case it is assumed
     *  that the data are known and there was no necessity to store
     *  them inside the container.
     *
     *  TODO save memory after deletion if capacity is really
     *  too large ?
     *
     *  @param col index of the Col
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void eraseRowsToCol( Integer const& col
                       , Integer const& pos
                       , Integer const& n=1)
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > col)
      { throw out_of_range("ITArray2D::eraseRowsToCol(col, pos, n) "
                                "this->firstCol() > col");
      }
      if (this->lastCol() < col)
      { throw out_of_range("ITArray2D::eraseRowsToCol(col, pos, n) "
                                "this->lastCol() < col");
      }
#endif
      // ptr on data
      TYPE* p_col(this->data(col));
      // range of the data (can be different from the range of the container)
      Range rangeCol(this->rangeCols_[col]);
      // number of rows to delete
      Integer nRowsToDelete = n;
      // position of the row to delete
      Integer firstRowToDelete = pos;
      // if pos is after the data there is nothing to delete
      if (rangeCol.first()>rangeCol.last()) return;
      // if pos is before the data
      if (firstRowToDelete<rangeCol.first())
      {
        // remove the not existing data !
        nRowsToDelete -= nRowsToDelete -rangeCol.first();
        // if there is no more rows to delete return
        if (nRowsToDelete<=0) return;
        // else shift data
        firstRowToDelete = rangeCol.first();
      }
      // if pos+n is after the data
      if (firstRowToDelete+nRowsToDelete-1>rangeCol.last())
      {
        // remove the not existing data !
        nRowsToDelete = rangeCol.last()-firstRowToDelete+1;
      }
      // update range
      this->rangeCols_[col].decLast(nRowsToDelete);
      rangeCol.decLast(nRowsToDelete);
      // get range of the data to delete
      const Integer last = rangeCol.last();
      // translate remaining Rows
      for ( Integer k=firstRowToDelete; k<=last; k++)
      {  p_col[k]   = p_col[k+nRowsToDelete];}
      // free mem if necessary
      if (this->rangeCols_[col].size()==0) freeCol(col);
    }

    /** @brief Internal method for deleting last rows to a
     * specified column.
     *
     *  Delete n last Rows to the container.
     *  TODO save memory if capacity is really too large ?
     *
     *  @param col index of the Col
     *  @param n number of elements to delete (default 1)
    **/
    void popBackRowsToCol( Integer const& col, Integer const& n=1)
    {
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->firstCol() > col)
      { throw out_of_range("ITArray2D::popBackRowsToCol(col, n) "
                           "this->firstCol() > col");
      }
      if (this->lastCol() < col)
      { throw out_of_range("ITArray2D::popBackRowsToCol(col, n) "
                           "this->lastCol() < col");
      }
#endif
      // check if there is something to do
      if (n <= 0) return;
      // check argument error
      if (n > this->rangeCol(col).size())
      {
        throw runtime_error("ITArray2D::popBackRowsToCol(col, n) "
                            "n > rangeCol(col).size().");
      }
      // update range
      this->rangeCols_[col].decLast(n);
      // free mem if necessary
      if (this->rangeCols_[col].size()==0)
        freeCol(col);
    }
};

} // namespace STK

#endif
// STK_ITARRAY2D_H
