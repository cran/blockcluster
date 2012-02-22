/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * created on: 25 nov. 2011
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_CArray2D.h
 *  @brief In this file we implement the final class CArray2D<TYPE>.
 **/

#ifndef STK_CARRAY2D_H
#define STK_CARRAY2D_H

#include "../../Sdk/include/STK_ContainerTraits.h"
#include "../../Sdk/include/STK_ITContainer2D.h"
#include "STK_AllocatorBase.h"
#include "STK_Array1D.h"
#include "STK_CArrayHo.h"

#include "STK_Display2D.h"

namespace STK
{
template<typename> class CArray2D;
/** @ingroup Arrays
 *  Specialization of the traits class ContainerTraits for Matrix
 *  class.
 */
template<class TYPE>
struct ContainerTraits<TYPE, CArray2D<TYPE> >
{
   typedef CArrayHo<TYPE> TContainerHo;
   typedef Array1D<TYPE> TContainerVe;
};

/** @brief A CArray2D is an two dimensional array
 *
 */
template < class TYPE>
class CArray2D: public ITContainer2D<TYPE, CArray2D<TYPE> >
              , public AllocatorBase<TYPE>
{
  public:
    /** Column type. */
    typedef Array1D<TYPE> TContainerVe;
    /** Row type. */
    typedef CArrayHo<TYPE> TContainerHo;
    /** 2D Interface type. */
    typedef ITContainer2D<TYPE, CArray2D<TYPE> > _ITContainer2DType_;
    /** AllocatorBase type. */
    typedef AllocatorBase<TYPE> _AllocatorBaseType_;

    /** Default constructor
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    CArray2D( Range const& I = Range(), Range const& J = Range())
            : _ITContainer2DType_(I, J)
    { initialize();}

    /** constructor with rbeg, rend, cbeg and cend specified,
     *  initialization with a constant.
     *  @param I range of the Rows
     *  @param J range of the Cols
     *  @param v initial value of the container
     **/
    CArray2D( Range const& I, Range const& J, TYPE const& v)
            : _ITContainer2DType_(I, J)
    {
      // initialize vertically the container
      initialize();
      // get dimension
      const Integer firstCol = J.first(), lastCol = J.last();
      const Integer firstRow = I.first(), lastRow = I.last();
      // initialize with v
      for (Integer j=firstCol; j<=lastCol; j++)
      {
        for (Integer i=firstRow; i<=lastRow; i++)
        { this->setData( j*index_ + i, v);}
      }
    }

    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    CArray2D( const CArray2D &T, bool ref=false)
            : _ITContainer2DType_(T)
            , _AllocatorBaseType_(T, ref)
            , index_(T.index_)
    {
      // if this is not a reference, AllocatorBase will not have allocated memory
      if (!ref)
      {
        // allocate memory
        initialize();
        // get dimension
        const Integer firstCol = T.firstCol(), lastCol = T.lastCol();
        const Integer firstRow = T.firstRow(), lastRow = T.lastRow();
        // copy data
        for (Integer j=firstCol; j<=lastCol; j++)
        {
          for (Integer i=firstRow; i<=lastRow; i++)
          { this->setData( j*index_ + i, T(i,j) );}
        }
      }
    }

    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    CArray2D( const CArray2D& T, Range const& I, Range const& J)
            : _ITContainer2DType_(I, J)
            , _AllocatorBaseType_(T, true)
            , index_(T.index_)
    { ;}

    /** Wrapper constructor for a C-like array.
     *  @param q pointer on the data
     *  @param nbRow number of rows
     *  @param nbCol number of columns
     **/
    CArray2D( TYPE* q, Integer const& nbRow, Integer const& nbCol)
            : _ITContainer2DType_(Range(0, nbRow-1), Range(0, nbCol-1))
            , _AllocatorBaseType_(q, nbRow*nbCol)
            , index_(nbRow)
    { ;}

    /** destructor. */
    virtual ~CArray2D()
    { ;}

    /** get a single element
     * @param i the index of the row
     * @param j the index of the column
     **/
    inline TYPE& elt( Integer const& i, Integer const& j)
    { return this->data(j*index_ + i);}

    /** get a single constant element
     * @param i the index of the row
     * @param j the index of the column
     **/
    inline TYPE const& elt( Integer const& i, Integer const& j) const
    { return this->data(j*index_ + i);}

    /** get a column of the array
     * @param j the index of the column
     **/
    Array1D<TYPE> col( Integer const& j) const
    { return Array1D<TYPE>(this->ptrData()+j*index_, this->rangeVe() ,j);}

    /** get a column of the array in a given range
     * @param I the range in the column we want to wrap
     * @param j the index of the column
     **/
    Array1D<TYPE> col( Range const& I, Integer const& j) const
    { return Array1D<TYPE>(this->ptrData()+j*index_, I ,j);}

//    CArrayHo<TYPE> row( Integer const& i) const;
//    CArrayHo<TYPE> row( Integer const& i, Range const& J) const;

    CArray2D<TYPE> sub( Range const& I, Range const& J) const
    { return CArray2D<TYPE>(*this, I,J);}

    /** New first indexes for the object.
     *  @param rbeg the index of the first row to set
     *  @param cbeg the index of the first column to set
     **/
    void shift(Integer const& rbeg =1, Integer const& cbeg =1)
    {
      { throw runtime_error("CArray2D::shift(n) "
          "not implemented.");
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
      { throw runtime_error("CArray2D::pushBackRows(n) "
          "not implemented.");
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
      // is this structure just a pointer?
      { throw runtime_error("CArray2D::insertRows(pos, n) "
          "not implemented.");
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
      { throw runtime_error("CArray2D::popBackRows() "
          "not implemented.");
      }
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
      { throw runtime_error("CArray2D::eraseRows(pos, n) "
          "not implemented.");
      }
    }

    /** Add n Cols to the container.
     *  @param n the number of Cols to add
     **/
    void pushBackCols(Integer const& n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      { throw runtime_error("CArray2D::pushBackCols(n) "
          "not implemented.");
      }
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
      { throw runtime_error("CArray2D::insertCols(pos, n) "
          "not implemented.");
      }
    }

    /** Delete last Cols of the container
     *  @param n the number of Cols to delete
     **/
    void popBackCols(Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      { throw runtime_error("CArray2D::popBackCols(n) "
          "not implemented.");
      }
    }

    /** Delete n Cols at the specified position of the container.
     *  @param pos the position of the deleted Cols
     *  @param n the number of column to delete
     **/
    void eraseCols(Integer const& pos, Integer const& n = 1)
    {
      if (n<=0) return;        // if n<=0 nothing to do
      { throw runtime_error("CArray2D::eraseCols(pos, n) "
                                 "not implemented.");
      }
    }


  private:
    /** index to apply in order to get the column. Will certainly be equal to
     *  the capacity_
     **/
    Integer index_;
    /** initialize the container with the dimension rangeHo_ and rangeVe_  */
    void initialize()
    {
      index_ = evalCapacity(this->rangeVe_.size());
      Integer const size = index_ * this->rangeHo_.size();
      if (this->ptrData())
      {
        this->reallocPtrData(size, compFirstData());
      }
      else
      {
        this->mallocPtrData(size, compFirstData());
      }
    }
    /** Compute the Increment we have to apply to the main pointer in order to
     * access to an element of the array using the easy to understand code
     * @code
     *   TYPE& elt(i,j)
     *   { return this->data(j*index_ + i);}
     * @endcode
     * The formula depend of the Horizontal range and of the Vertical range
     * of the container and is given by the formula.
     **/
    Integer compFirstData() const
    { return this->rangeVe_.size()*this->rangeHo_.first()+this->rangeVe_.first();}
};

/** @ingroup Arrays
 *  ostream for CArray2D.
 *  @param s the output stream
 *  @param V the CArray2D to write
 **/
template<class TYPE>
ostream& operator<<(ostream& s, const CArray2D<TYPE>& V)
{ return out2D(s,V);}

}

#endif /* STK_CARRAY2D_H_ */
