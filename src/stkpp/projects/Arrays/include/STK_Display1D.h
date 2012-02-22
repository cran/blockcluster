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
 * Project: stkpp::Arrays
 * Purpose:  Display on the standard output the one dimensional arrays.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Display1D.h
 *  @brief This file define very basic methods for displaying
 *  one dimensional Containers.
 **/

#ifndef STK_DISPLAY1D_H
#define STK_DISPLAY1D_H

//#include <iomanip>

#include "../../STKernel/include/STK_Proxy.h"
#include "../../STKernel/include/STK_Integer.h"
#include "../../STKernel/include/STK_Real.h"
#include "../../STKernel/include/STK_Stream.h"

namespace STK
{

/** @ingroup Arrays
 *  Write a Container1D.
 *  @param s the output stream
 *  @param V the Container to write
 **/
template<class TYPE, template<class> class Container1D>
ostream& out1D( ostream& s, const Container1D<TYPE>& V)
{
  Integer rbeg = V.first(), rend = V.last();

  s << std::right;
//    << std::setprecision(1);
  for (Integer i=rbeg; i<=rend; i++)
  { s << ConstProxy<TYPE>(V[i]) << _T(' ') ; }
  return s;
}

// Specialization
/** @ingroup Arrays
 * Specialization for Real.
 *  @param s the output stream
 *  @param V the Container to write
 **/
template<template<class> class Container1D>
ostream& out1D(ostream& s, const Container1D<Real>& V)
{
  Integer rbeg = V.first(), rend = V.last();

  s << std::right;
//    << std::setprecision(1);
  for (Integer i=rbeg; i<=rend; i++)
  { s << /* std::setw(12) << */ ConstProxy<Real>(V[i]) << _T(' ') ; }

  return s;

}

/** @ingroup Arrays
 * Specialization for Integer.
 *  @param s the output stream
 *  @param V the Array1D to write
 **/
template<template<class> class Container1D>
ostream& out1D(ostream& s, const Container1D<Integer> & V)
{
  Integer rbeg = V.first(), rend = V.last();

  s << std::right;
  for (Integer i=rbeg; i<=rend; i++)
  { s << /* std::setw(5) << */ ConstProxy<Integer> (V[i]) << _T(' ') ; }
  return s;
}

} // namespace STK

#endif // STK_DISPLAY_H
