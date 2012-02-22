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
 * Purpose:  Define the fundamental type Binary.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Binary.h
 *  @brief In this file we define the fundamental type Binary.
 **/

#ifndef STK_BINARY_H
#define STK_BINARY_H

#include "STK_Arithmetic.h"
#include "STK_IdTypeImpl.h"
#include "STK_Stream.h"
#include "STK_Proxy.h"

namespace STK
{

/** @ingroup Base
 *  @brief STK fundamental type of a binary.
 *
 *  The type Binary is a representation of dichotomic variables.
 **/
 enum Binary
 { zero =0, ///< zero value
   one  =1, ///< one value
   binaryNA ///< Not Available value
 };

 /** @ingroup Arithmetic
  *  @brief Specialization for Binary.
  * 
  *  NA (not available) numbers is part of the @c union.
  */
 template<>
 struct Arithmetic<Binary> : public std::numeric_limits<Binary>
 {
   /** Adding a Non Avalaible (NA) special number.
    **/
   static Binary NA() throw()
   { return binaryNA;}

   /** True if the type has a representation for a "Not Available".
    **/
   static const bool hasNA = true;

   /** Test if x is a Non Avalaible (NA) special number
    *  @param x the Binary number to test.
    **/
   static bool isNA(const Binary& x) throw()
   { return (x==binaryNA);}
     
   /** test if x is  infinite.
    *  @param x the Binary number to test.
    **/
   static bool isInfinite(const Binary& x) throw()
   { return false; }

   /** test if x is  finite.
    *  @param x the Binary number to test.
    **/
   static bool isFinite(const Binary& x) throw()
   { return (!isNA(x) && !isInfinite(x));}
 };

/** @ingroup RTTI 
 *  @brief Specialization of the IdTypeImpl for the Type Binary.
 * 
 *  Return the IdType of a Binary.
 **/
template<>
struct IdTypeImpl<Binary>
{
  /** Give the IdType of the type Binary. */
  static IdType returnType()
  { return(binary);}
};
  
/** @ingroup stream
 *  @brief Overloading of the ostream << for the type Binary.
 *  @param os the output stream
 *  @param output the value to send to the stream
 **/
ostream& operator << (ostream& os, const Binary& output);

 /** @ingroup stream
  *  @brief Overloading of the istream >> for the type Binary.
  *  @param is the input stream
  *  @param input the value to get from the stream
  **/
istream& operator >> (istream& is, Proxy<Binary>& input);

} // namespace STK

#endif /*STK_BINARY_H*/
