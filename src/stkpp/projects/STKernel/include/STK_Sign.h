/*--------------------------------------------------------------------*/
/*     Copyright (C) 2007  Serge Iovleff

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
 * Project:  Base
 * Purpose:  Define the fundamental type Sign.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Sign.h
 *  @brief In this file we define the fundamental type Sign.
 **/

#ifndef STK_SIGN_H
#define STK_SIGN_H

#include "STK_Arithmetic.h"
#include "STK_IdTypeImpl.h"
#include "STK_Stream.h"
#include "STK_Proxy.h"

namespace STK
{
/** @ingroup Base
 *  @brief STK fundamental type of a sign.
 *
 *  The type Signed is an other representation of dichotomic variables
 **/
 enum Sign
 { negative =-1, ///< negative sign
   positive = 1, ///< positive sign
   SignNA     ///< Not Available value
 };

/** @ingroup Arithmetic
 *  @brief Specialization for Sign.
 *
 *  NA (not available) numbers is part of the @c union.
 */
template<>
struct Arithmetic<Sign> : public std::numeric_limits<Sign>
{
   /** Adding a Non Avalaible (NA) special number.
    **/
   static Sign NA() throw()
   { return SignNA;}

   /** True if the type has a representation for a "Not Available".
    **/
   static const bool hasNA = true;

   /** Test if x is a Non Avalaible (NA) special number.
    *  @param x the Sign number to test.
    **/
   static bool isNA(const Sign& x) throw()
   { return (x==SignNA);}
     
   /** Test if x is  infinite.
    *  @param x the Sign number to test.
    **/
   static bool isInfinite(const Sign& x) throw()
   { return false; }

   /** Test if x is  finite.
    *  @param x the Sign number to test.
    **/
   static bool isFinite(const Sign& x) throw()
   { return (!isNA(x) && !isInfinite(x));}
 };
 
/** @ingroup RTTI 
 *  @brief Specialization of the IdTypeImpl for the Type Sign.
 * 
 *  Return the IdType of a Sign.
 **/
template<>
struct IdTypeImpl<Sign>
{
  /** Give the IdType of the type Sign. */
  static IdType returnType()
  { return(signof);}
};
  
/** @ingroup stream
 *  @brief Overloading of the ostream << for the type Sign.
 *  @param os the output stream
 *  @param output the value to send to the stream
 **/
ostream& operator << (ostream& os, const Sign& output);

/** @ingroup stream
 *  @brief Overloading of the istream >> for the type Sign.
 *  @param is the input stream
 *  @param input the value to get from the stream
 **/
istream& operator >> (istream& is, Proxy<Sign>& input);

} // namespace STK

#endif /*STK_SIGN_H*/
