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
 * Purpose:  implement the Range class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Range.cpp
 *  @brief In this file we implement the class Range, a range of Index
 *  denoting a sub-vector region.
 **/

#include "../include/STK_Range.h"
#include "../include/STK_Misc.h"

#include "../include/STK_String_Util.h"

namespace STK
{

/** The Default constructor Assume the beginning of the sub-region is 1.
 * @param last the end of the sub-region.
 **/
Range::Range( Integer const& last)
        : first_(1)
        , last_(last)
        , size_(last)
{ ;}

/** Complete constructor We have to give the beginning and the end of
 *  the sub-region.
 *  @param first the beginning of the sub-region
 *  @param last the end of the sub-region.
**/
Range::Range( Integer const& first, Integer const& last)
        : first_(first)
        , last_(last)
        , size_(last-first+1)
{ ;}


/** Copy constructor Create a copy of an existing Range.
 *  @param I The index to copy
**/
Range::Range( Range const& I)
        : first_(I.first())
        , last_(I.last())
        , size_(I.size())
{ ;}

Range::~Range() { ;}

/*--------------------------------------------------------------------*/
/* Manipulation of the index.                                         */
/** Set new range for the index.
 * @param first first element
 * @param last  last element
**/
Range& Range::set(Integer const& first, Integer const& last)
{
  first_ = first;
  last_  = last;
  size_     = last - first + 1;
  return *this;
}

/** Shift the Range giving the first element : the size is not modified.
 *  @param first new value of the first element
**/
Range& Range::shift(Integer const& first)
{
  return inc(first - first_); // use the inc method
}

/** Increase first_ and last_
 *  @param inc the increment to apply
**/
Range& Range::inc(Integer const& inc)
{
  first_ +=inc;
  last_  +=inc;
  return *this;
}
/** Increase first_
  * @param inc the increment to apply
 **/
Range& Range::incFirst(Integer const& inc)
{
  first_ +=inc;
  size_     -=inc;
  return *this;
}
/** Increase last_
  * @param inc the increment to apply
 **/
Range& Range::incLast(Integer const& inc)
{
  last_ +=inc;
  size_    +=inc;
  return *this;
}

/** Decrease first_ and last_
 *  @param dec the decrement to apply
**/
Range& Range::dec(Integer const& dec)
{
  first_ -=dec;
  last_  -=dec;
  return *this;
}
/** Decrease first_
  * @param dec the decrement to apply
 **/
Range& Range::decFirst(Integer const& dec)
{
  first_ -=dec;
  size_     +=dec;
  return *this;
}
/** Decrease last_
  * @param dec the decrement to apply
 **/
Range& Range::decLast(Integer const& dec)
{
  last_ -=dec;
  size_    -=dec;
  return *this;
}

/** Take the lowest value of first_ and I.first_ for first_
 *  and the largest value of last_ and I.last_ for last_.
 *  @param I the index to apply
**/
Range& Range::sup(Range const& I)
{
  first_ = min(first_, I.first_);
  last_  = max(last_, I.last_);
  size_     = last_ - first_ +1;
  return *this;
}

/** Take the largest value of first_ and I.first_ for first_
 *  and the lowest value of last_ and I.last_ for last_.
 *  @param I the index to apply
**/
Range& Range::inf(Range const& I)
{
  first_ = max(first_, I.first_);
  last_  = min(last_, I.last_);
  size_     = last_ - first_ +1;
  return *this;
}
/** Take the lowest value of I.first_ and J.first_ for first_
 *  and the largest value of I.last_ and J.last_ for last_.
 *  @param I first the index to apply
 *  @param J second the index to apply
*/
Range Range::sup(Range const& I, Range const& J)
{
  return Range(min(I.first_, J.first_), max(I.last_, J.last_));
}

/** Take the largest value of I.first_ and J.first_ for first_
 *  and the lowest value of I.last_ and J.last_ for last_.
 *  @param I first the index to apply
 *  @param J second the index to apply
 */
Range Range::inf(Range const& I, Range const& J)
{
  return Range(max(I.first_, J.first_), min(I.last_, J.last_));
}

/** @brief create the Range [first_+inc, last_+inc_]. **/
Range Range::plus(Integer const& inc)
{
  return *this+inc;
}
/** @brief create the Range [first_-dec, last_-dec]. **/
Range Range::minus(Integer const& dec)
{
  return *this-dec;
;
}

/** Increase first_ and last_
 * @param inc the increment to apply
**/
Range& Range::operator+=(Integer const& inc)
{ first_ += inc; last_ += inc; return *this;}

/** Decrease first_ and last_
 *  @param dec the decrement to apply
 **/
Range& Range::operator-=(Integer const& dec)
{ first_ -= dec; last_ -= dec; return *this;}

/** Print a Range using the colon notation.
 *  @param s output stream
 *  @param I the Range to write
 **/
ostream& operator<< (ostream& s, Range const& I)
{
  s << I.first_ << _T(":") << I.last_;
  return s;
}

/** Read a Range in the form first:last (MATLAB-like form) from
 * an input stream. The input stream can also be a number (say n). In this case
 * the range will be 1:n. If the range cannot be read the method return the
 * default range 1:0. FIXME? should be NA value, but NA value is not defined
 * for Range.
 * @param  s the input stream
 * @param I the range to set
 **/
istream& operator>> (istream& s, Range& I)
{
  String num;
  s >> std::skipws;
  // get first number
  std::getline(s, num, _T(':'));
  // check if the istream is exhausted
  if (s.eof())
  {
    I.first_ = 1;
    if (!stringToType(I.last_, num)) I.last_ =0;
    return s;
  }
  // otherwise we encounter a ":", thus skip the current char
  if (!stringToType(I.first_, num)) I.first_ =1;
  s.peek();
  if ((s >> I.last_).fail())
  {
    I.first_ =1; I.last_ =0;
  }
  return s;
}

/** @brief return an Index applying @c Range::dec() to I **/
Range operator-(Range const& I, Integer const& dec)
{
  return(Range(I.first() - dec, I.last() - dec));
}

/** @brief return an Index applying @c Range::inc() to I **/
Range operator+(Range const& I, Integer const& inc)
{
  return(Range(I.first() + inc, I.last() + inc));
}


} // Namespace STK
