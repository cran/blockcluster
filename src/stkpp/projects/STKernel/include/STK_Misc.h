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
 * Purpose:  Define miscenaleous utility functions.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Misc.h
 *  @brief In this file we define miscellaneous utility (templated)
 *  functions.
 **/

#ifndef STK_MISC_H
#define STK_MISC_H

#include <cmath>
#include <cstdlib> // for rand

#include "STK_String.h"
#include "STK_Integer.h"
#include "STK_Real.h"
#include "STK_Stream.h"
#include "STK_Proxy.h"

namespace STK
{
/** @ingroup Base
 *
 *  templated sign value sign(x) * y:
 *  TYPE should be Integer , long, float or Real
 *  @param x the sign value
 *  @param y the signed value to return
 **/
template<class TYPE>
inline TYPE sign(TYPE const& x, TYPE const& y = 1)
{ return( (x<0) ? -y : y); }

/** @ingroup Base
 *
 *  templated swap method.
 *  @param x the first value to swap
 *  @param y the second value to swap
 **/
template<class TYPE>
inline void swap(TYPE& x, TYPE& y)
{
  TYPE aux(x);
  x = y;
  y = aux;
}

/** @ingroup Base
 *
 *  templated minimum for class having defined <.
 *  @param x the first value
 *  @param y the second value
 **/
template<class TYPE>
inline TYPE const& min(TYPE const& x, TYPE const& y)
{ return( (x<y) ? x : y); }

/** @ingroup Base
 *
 *  templated maximum for class having defined <.
 *  @param x the first value
 *  @param y the second value
 **/
template<class TYPE>
inline TYPE const& max(TYPE const& x, TYPE const& y)
{ return( (x<y) ? y : x); }

/** @ingroup Base
 *
 *  templated absolute value: TYPE should be Integer , long, float
 *  or Real.
 *  @param x the value
 **/
template<class TYPE>
inline TYPE abs(TYPE const& x)
{ return( (x<0) ? -x : x); }

/** @ingroup Base
 *
 *  frand() generate a Real uniform number. This is a very basic method
 *  and should only be used when speed is necessary only.
 **/
inline Real frand()
{ return (Real)rand() / (RAND_MAX+1.0);}

/** @ingroup Base
 *  @brief is x a odd number ?
 *
 *  This method return true if the rest of the euclidian division
 *  of x by 2 is 0.
 *  @param x the value to test
 *  @return @c true if x is odd, @c false otherwise
 **/
inline bool isOdd(Integer const& x)
{ return( (x%2) == 1 ); }

/** @ingroup Base
 *  @brief is x an even number ?
 *
 *  This method return true if the rest of the euclidian division
 *  of x by 2 is 1.
 *  @param x the value to test
 *  @return @c true if x is even, @c false otherwise
 **/
inline bool isEven(Integer const& x)
{ return( (x%2) == 0 ); }

/** @ingroup Base
 *
 *  Computation of round off : return an Integer value
 *  @param x the value to round
 *  @return the rouded value of x
 **/
inline Integer round(Real const& x)
{ return( x < 0.0 ? Integer(x-0.5) : Integer(x+0.5));}

/** @ingroup Base
 *
 *  Computation of sqrt(x^2 + y^2) without underflow or overflow.
 *  @param x first value
 *  @param y second value
 *  @return the value \f$ \sqrt(x^2 + y^2) \f$
 **/
inline Real norm(Real const& x, Real const& y)
{
  Real absx = abs(x), absy = abs(y);
  if (absx > absy)
  {
    return(absx * sqrt(double(1.0+(absy/absx)*(absy/absx))));
  }
  else
  {
    return(absy == 0.0 ? 0.0 : absy * sqrt(double(1.0+(absx/absy)*(absx/absy))));
  }
}

} // namespace STK

#endif /* STK_MISC_H */
