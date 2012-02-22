/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2009  Serge Iovleff

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
 * Project:  stkpp::STKernel::Base
 * Purpose:  Define the Range class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Range.h
 *  @brief In this file we define the class Range, a range of Index
 *  denoting a sub-vector region.
 **/

#ifndef STK_RANGE_H
#define STK_RANGE_H

#include "STK_Integer.h"
#include "STK_Stream.h"

namespace STK
{
/** @ingroup TContainer
 *  @brief Index sub-vector region.
 *
 *  An Range is an ordered pair [first,last] denoting a sub-vector
 *  region, similar to a Fortran 90 or Matlab colon notation.
 *  For example :
 *  @code
 *  Vector A(10), B(0,20);
 *
 *  Range I(2,4);
 *
 *  A(I) = B(Range(0,2));
 *  @endcode
 *  overwrite the elements 2, 3 and 4 of A by the elements 0, 1
 *  and 2 of B.  There is no stride argument, only contiguous regions
 *  are allowed.
 **/
class Range
{
  private:
    Integer first_;    ///< First index
    Integer last_;     ///< Last index
    Integer size_;     ///< Theoretic Dimension size_ = last_- first_ +1

  public:
    /** @brief Default constructor. The first index is 1.
     *  @param last the last index of the range
     *  **/
    Range(Integer const& last =0);

    /** @brief Complete constructor.
     *  @param first is the beginning of the sub-region
     *  @param last is the end of the sub-region.
     **/
    Range(Integer const& first, Integer const& last);

    /** @brief Copy constructor
     *  @param I The index to copy
     **/
    Range(Range const& I);

    /** destructor. */
     ~Range();

    /** get the first index of the Range
     *  @return the first index of the range
     **/
    inline Integer const& first()  const { return first_;};
    /** get the last index of the Range
     *  @return the last index of the range
     **/
    inline Integer const& last()    const { return last_;};
    /** get the size of the Range (the number of elements)
     *  @return the size of the range
     **/
    inline Integer const& size()   const { return size_;};
    /** check if the range is empty or not
     *  @return @c true if size <=0, @c false otherwise
     */
    inline bool    empty() const { return size_<=0;};
    /** check if this Range in include in an other Range
     *  @param I the index to compare
     *  @return @c true if this is include in @c I, @c false otherwise
     **/
    inline bool isIn(Range const& I) const
    { return ((first_>= I.first_)&&(last_<=I.last_));}
    /** check if the Range I is include in the this Range
     *  @param I the index to compare
     *  @return @c true if this contain I, @c false otherwise
     **/
    inline bool isContaining(Range const& I) const
    { return ((first_<= I.first_)&&(last_>=I.last_));}
    /** Return true if i is in this Range
     *  @param i the integer to compare
     *  @return @c true if i is in this, @c false otherwise
     **/
    inline bool isContaining(Integer const& i) const
    { return ((first_<=i)&&(last_>=i));}
    /** compare this range with range @c I
     *  @param I the Index to compare
     *  @return @c true if the range are equals, @c false otherwise
     **/
    inline bool operator==(Range const& I) const
    { return ((first_ == I.first()) && (last_ == I.last()));}
    /** compare this range with range @c I
     *  @param I the Index to compare
     *  @return @c true if the range are different, @c false otherwise
     **/
    inline bool operator!=(Range const& I) const
    { return ((first_ != I.first()) || (last_ != I.last()));}

    /** @brief Set dimensions of the Index. **/
    Range& set(Integer const& first =1, Integer const& last =0);
    /** @brief create the Range [first, first+size_]. **/
    Range& shift(Integer const& first =1);
    /** @brief create the Range [first_+inc, last_+inc_]. **/
    Range& inc(Integer const& inc =1);
    /** @brief create the Range [first_+inc, last_]. **/
    Range& incFirst(Integer const& inc =1);
    /** @brief create the Range [first_, last_+inc]. **/
    Range& incLast(Integer const& inc =1);
    /** @brief create the Range [first_-dec, last_-dec]. **/
    Range& dec(Integer const& dec =1);
    /** @brief create the Range [first_-dec, last_]. **/
    Range& decFirst(Integer const& dec =1);
    /** @brief create the Range [first_, last_-dec]. **/
    Range& decLast(Integer const& dec =1);
    /** @brief compute sup(this,J) **/
    Range& sup(Range const& I);
    /** @brief compute inf(this,J) **/
    Range& inf(Range const& I);
    /** @brief same as @c Range::inc() **/
    Range& operator+=(Integer const& inc);
    /** @brief same as @c Range::dec() **/
    Range& operator-=(Integer const& dec);

    /** @brief create the Range [first_+inc, last_+inc_]. **/
    Range plus(Integer const& inc =1);
    /** @brief create the Range [first_-dec, last_-dec]. **/
    Range minus(Integer const& dec =1);

    /** @brief compute sup(I,J) **/
    static Range sup(Range const& I, Range const& J);
    /** @brief compute inf(I,J)I **/
    static Range inf(Range const& I, Range const& J);
    /** @brief Write a Range in the form first:last (MATLAB-like form) in an
     * output stream. **/
    friend ostream& operator<< (ostream& s, Range const& I);
    /** @brief Read a Range in the form first:last (MATLAB-like form) from
     * an input stream. **/
    friend istream& operator>> (istream& s, Range& I);
};

/** @brief return an Index applying @c Range::dec() to I
 *  @param I the original range
 *  @param dec the decrement to apply
 **/
Range operator-(Range const& I, Integer const& dec);

/** @brief return an Index applying @c Range::inc() to I
 *  @param I the original range
 *  @param inc the increment to aply
 **/
Range operator+(Range const& I, Integer const& inc);

} // Namespace STK

#endif // STK_RANGE_H
