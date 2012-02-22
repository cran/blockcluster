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
 * Project:  stkpp::Regress
 * created on: 23 juin 2011
 * Purpose:  num and other utilities for the Regress project..
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Regress_Util.h
 *  @brief In this file we declare the main .
 **/


#ifndef STK_REGRESS_UTIL_H_
#define STK_REGRESS_UTIL_H_

#include "../../STKernel/include/STK_String.h"

namespace STK
{

namespace Regress
{

/** @ingroup Regress
 * Regression method we will use. */
enum TypeRegression
{
  /** unknown regression */
  unknown_ = 0
  /** linear regression */
  , linear_
  /** additive BSpline regression */
  , additiveBSpline_
  /** adaptive BSpline regression */
  , adaptiveBSpline_
};

/** @ingroup Regress
 *  Convert a String to a TypeRegression.
 *  @param type the String we want to convert
 *  @return the TypeRegression represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
TypeRegression StringToTypeRegression( String const& type);

/** @ingroup Regress
 *  Convert a TypeRegression to a String.
 *  @param type the type of regression we want to convert
 *  @return the string associated to this type.
 **/
String TypeRegressionToString( TypeRegression const& type);

} // namespace Regress

} //namespace STK

#endif /* STK_REGRESS_UTIL_H_ */
