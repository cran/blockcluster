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
 * Project:  stkpp::Sdk::TContainer
 * Purpose:  Define the interface 2D Container classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IContainer2D.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ICONTAINER2D_H
#define STK_ICONTAINER2D_H

#include "../../STKernel/include/STK_Range.h"

namespace STK
{
/** @ingroup TContainer
 * @brief Interface class for 2D containers.
 *
 * The IContainer2D class is the base class for all two-dimensional
 * containers.
 *
 * The pure virtual function defined in this interface have the
 * following definition:
 * @code
 *   virtual void shift(Integer const& rbeg =1, Integer const& cbeg =1);
 *   virtual void pushBackCols(Integer const& n =1);
 *   virtual void popBackCols(Integer const& n =1);
 *   virtual void pushBackRows(Integer const& n =1);
 *   virtual void popBackRows(Integer const& n =1);
 *   virtual void eraseCols(Integer const& pos, Integer const& n = 1);
 *   virtual void eraseRows(Integer const& pos, Integer const& n=1);
 * @endcode
 * and they are used in the @c resize(I,J) method.
 **/
class IContainer2D
{
  protected:
    /** Default constructor. rangeHo_ = 1:0 and rangeVe_ = 1:0.
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    IContainer2D( Range const& I = Range(), Range const& J =  Range())
                : rangeHo_(J)
                , rangeVe_(I)
    { ;}

    /** Copy constructor
     *  @param T the container to copy
     **/
    IContainer2D( const IContainer2D &T)
                : rangeHo_(T.rangeHo_)
                , rangeVe_(T.rangeVe_)
    { ;}

    /** destructor. **/
    ~IContainer2D() { ;}

  public:
    /** Index of the first column.
     *  @return the index of the first column
     **/
    inline Integer const& firstCol() const { return rangeHo_.first();}
    /** Index of the last column.
     *  @return the index of the last column
     **/
    inline Integer const& lastCol() const { return rangeHo_.last();}
    /** get the number of columns.
     * @return the Horizontal size (the number of column)
     **/
    inline Integer const& sizeHo() const { return rangeHo_.size();}
    /** Range of the column of the container.
     * @return the Horizontal range
     **/
    inline Range const& rangeHo() const { return rangeHo_;};
    /** Index of the first row.
     * @return the index of the first row
     **/
    inline Integer const& firstRow() const { return rangeVe_.first();}
    /** Index of the last row.
     * @return the index of the last row
     **/
    inline Integer const& lastRow() const { return rangeVe_.last();}
    /** get the number of rows.
     * @return the Vertical size (the number of rows)
     **/
    inline Integer const& sizeVe() const { return rangeVe_.size();}
    /**  Range of the rows of the container.
     *  @return the Vertical range
     **/
    inline Range const& rangeVe() const { return rangeVe_;}
    /** Is there some data ?
     *  @return @c true if the container is empty, @c false otherwise
     **/
    inline bool empty() const { return (rangeHo_.empty() || rangeVe_.empty());}

    /** resize the container:
     * - call @c shift(I.first(), J.first()
     * - call @c popBackCols() (@c insertRows()) and/or @c popBackCols()
     *  (@c popBackRows()).
     *  The implicit assumption made by this method is that it is easier and
     *  faster to add column than add rows to the 2D container.
     * @param I the new range for the rows of the container
     * @param J the new range for the columns of the container
     **/
    void resize( Range const& I, Range const& J)
    {
      // check if there is something to do
      if ((rangeVe_ == I) && (rangeHo_ == J)) return;
      //  translate beg
      shift(I.first(), J.first());
      // number of rows to del or add
      Integer rinc = I.last() - lastRow();
      // number of cols to del or add
      Integer cinc = J.last() - lastCol();
      // check if we add cols
      if (cinc >=0)   // work first on rows
      {
        if (rinc < 0) popBackRows(-rinc); // less rows
        else          pushBackRows(rinc); // more rows
        pushBackCols(cinc); // add columns
      }
      else // work first on columns
      {
        popBackCols(-cinc); // remove columns
        if (rinc < 0) popBackRows(-rinc); // less rows
        else          pushBackRows(rinc); // more rows
      }
    }

    /** Delete n Cols at the specified position to the container.
      *  @param pos the position of the columns to delete
      *  @param n the number of column to delete
      **/
     virtual void eraseCols(Integer const& pos, Integer const& n = 1) =0;

     /** Delete n Rows at the pos index to the container.
      *  @param pos the position of the rows to delete
      *  @param n number of rows to delete (default 1)
     **/
     virtual void eraseRows(Integer const& pos, Integer const& n=1) =0;

     /** remove the n last columns to the container.
      * @param n number of columns to remove to the container
      **/
     virtual void popBackCols( Integer const& n =1) =0;

   /** remove the n last rows to the container.
     * @param n number of rows to delete to the container
     **/
    virtual void popBackRows( Integer const& n =1) =0;

  protected:
    Range rangeHo_;    ///< Range of the cols
    Range rangeVe_;    ///< Range of the rows

    /** Set range of the container.
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline void setRange(Range const& I = Range(), Range const& J = Range())
    {
      setRangeVe(I);
      setRangeHo(J);
    }
    /** Set the range of the number of rows.
     *  @param I the range of the rows number
     **/
    inline void setRangeVe( Range const& I = Range())
    { rangeVe_ = I;}
    /** Increment the range of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incRangeVe( Integer const& inc)
    { rangeVe_.inc(inc);}
    /** Decrement range of the number of rows.
     *  @param dec the decrement to apply
     **/
    inline void decRangeVe( Integer const& dec)
    { rangeVe_.dec(dec);}
    /** Set the beginning of the rows.
     *  @param beg the first index of the rows
     **/
    inline void setFirstVe( Integer const& beg)
    { rangeVe_.shift(beg);}
    /** Increment the beginning of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incFirstVe( Integer const& inc)
    { rangeVe_.incFirst(inc);}
    /** Decrement the beginning of the number of rows.
     *  @param dec the decrement to apply
     **/
    inline void decFirstVe( Integer const& dec)
    { rangeVe_.decFirst(dec);}
    /** Increment the end of the number of rows.
     *  @param inc the increment to apply
     **/
    inline void incLastVe( Integer const& inc)
    { rangeVe_.incLast(inc);}
    /** Decrement the end of the number of rows.
     *  @param dec the decrement to apply
     **/
    inline void decLastVe( Integer const& dec)
    { rangeVe_.decLast(dec);}
    /** Set the range of the columns.
     * @param J the range of the cols number
     **/

    inline void setRangeHo( Range const& J = Range())
    { rangeHo_ = J;}
    /** Increment the range of the number of columns.
     *  @param inc the increment to apply
     **/
    inline void incRangeHo( Integer const& inc)
    { rangeHo_.inc(inc);}
    /** Decrement range of the number of columns.
     *  @param dec the decrement to apply
     **/
    inline void decRangeHo( Integer const& dec)
    { rangeHo_.dec(dec);}
    /** Set the beginning of the number of columns.
     *  @param beg the new first index
     **/
    inline void setFirstHo( Integer const& beg)
    { rangeHo_.shift(beg);}
    /** inc the beginning of the number of columns.
     *  @param inc the increment to apply
     **/
    inline void incFirstHo( Integer const& inc)
    { rangeHo_.incFirst(inc);}
    /** Decrement the beginning of the number of columns.
     *  @param dec the decrement to apply
     **/
    inline void decFirstHo( Integer const& dec)
    { rangeHo_.decFirst(dec);}
    /** Increment the end of the number of columns.
     *  @param inc the increment to apply
     **/
    inline void incLastHo( Integer const& inc)
    { rangeHo_.incLast(inc);}
    /** Decrement the end of the number of columns.
     *  @param dec the decrement to apply
     **/
    inline void decLastHo( Integer const& dec)
    { rangeHo_.decLast(dec);}

    /** swap this container with T
     * @param T the container to swap with T
     **/
     inline void swap(IContainer2D& T)
     {
       // get ranges of this
       Range range_ho(rangeHo());
       Range range_ve(rangeVe());
       // set range of T to this
       setRange(T.rangeVe(), T.rangeHo());
       // set range of this to T
       T.setRange(range_ve, range_ho);
     }

  private:
    /** Pure virtual method. Shift the ranges of the container.
     * @param rbeg the first index for the rows of the container
     * @param cbeg the first index for the columns of the container
     **/
    virtual void shift( Integer const& rbeg, Integer const& cbeg) =0;

    /** Pure virtual method. Add n columns to the container.
     * @param n number of columns to add to the container
     **/
    virtual void pushBackCols( Integer const& n =1) =0;

     /** Pure virtual method. Add n columns to the container.
     * @param n number of columns to add to the container
     **/
    virtual void pushBackRows( Integer const& n =1) =0;
};

} // namespace STK

#endif // STK_ICONTAINER2D_H
