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
 * Project:  Base
 * Purpose:  Define the fundamental type NotAvailable.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_NotAvailable.h
 *  @brief In this file we define the fundamental types
 *  @c NotAvailable.
 **/

#ifndef STK_NOTAVAILABLE_H
#define STK_NOTAVAILABLE_H

#include "STK_Arithmetic.h"
#include "STK_IdTypeImpl.h"
#include "STK_Stream.h"
#include "STK_Proxy.h"

namespace STK
{
/** @ingroup Base
  *  @brief STK fundamental type of a NA value.
  *
  *  The type NotAvailable is defined for the numerical computation and
  *  the internal representation of the not Available variables.
  *  It can be used for any unknown variable.
  **/
enum NotAvailable
{ missing };

/** @ingroup Arithmetic
 *  @brief Specialization for NotAvailable.
 * 
 *  The NA Type is for variable always Not Available. It is thus
 *  never available.
 */
template<>
struct Arithmetic<NotAvailable> : public std::numeric_limits<NotAvailable>
{
  /** Adding a Non Avalaible (NA) special number.
   **/
  static NotAvailable NA() throw()
  { return missing;}

  /** True if the type has a representation for a "Not Available".
   **/
  static const bool hasNA = true;

  /** Test if x is a Non Avalaible (NA) special number.
   *  @param x the NA value to test.
   **/
  static bool isNA(const NotAvailable& x) throw()
  { return true;}
    
  /** Test if x is  infinite.
   *  @param x the NA value to test.
   **/
  static bool isInfinite(const NotAvailable& x) throw()
  { return false; }

  /** Test if x is  finite.
   *  @param x the NA value to test.
   **/
  static bool isFinite(const NotAvailable& x) throw()
  { return false;}
};

/** @ingroup RTTI 
 *  @brief Specialization of the IdTypeImpl for the Type NotAvailable.
 * 
 *  Return the IdType of a NotAvailable.
 **/
template<>
struct IdTypeImpl<NotAvailable>
{
  /** Give the IdType of the type NotAvalailable. */
  static IdType returnType()
  { return(notavailable);}
};

/** @ingroup stream
 *  @brief Overloading of the ostream << for the type NotAvailable.
 *  @param os the output stream
 *  @param output the value to send to the stream
 **/
ostream& operator << (ostream& os, const NotAvailable& output);

/** @ingroup stream
 *  @brief Overloading of the istream >> for the type NotAvailable.
 *  @param is the input stream
 *  @param input the value to get from the stream
 **/
istream& operator >> (istream& is, Proxy<NotAvailable>& input);

} // namespace STK

#endif /* STK_NOTAVAILABLE_H */
