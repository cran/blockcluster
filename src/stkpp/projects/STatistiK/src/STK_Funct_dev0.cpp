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
 * Project:  Analysis
 * Purpose:  Compute the function dev0
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Funct_g0.cpp
 *  @brief In this file we implement the function g0.
 **/

#include <cmath>

// for Integer and Real
#include "../../STKernel/include/STK_Integer.h"
#include "../../STKernel/include/STK_Real.h"
#include "../../STKernel/include/STK_Misc.h"

#include "../include/STK_Funct_util.h"

namespace STK
{

namespace Funct
{
/** @ingroup Analysis
 *  Compute the partial deviance:
 * \f[
 *     d_0(a, b) = a\log(a/b)+ b - a
 * \f]
 * with good accuracy for a/b is near 1 using the Taylor series of
 * \f[ \log\left( \frac{a-b}{a+b}\right) \f]
 *  a should be >=0 and b should be > 0, no check is done.
 * @see http://www.herine.net/stat/software/dbinom.html
 * @param a first parameter
 * @param b second paramater
 * @return the evaluation of the function d_0
 **/
Real dev0(Real const& a, Real const& b)
{
  // special values
  if (a == 0.) return b;
  if (b == a)  return 0.;
  // ratio
  Real v = (a-b)/(a+b);
  // small v -> use dl
  if (abs(v)<0.125)
  {
    Real sum  = (a-b)*v;
    Real ej   = 2*a*v;
    Real term;
    v = v*v;
    Integer n =0;
    do
    {
      ej *= v;
      sum += (term = (ej/(2*(++n)+1)));
    }
    while (abs(term) > abs(sum) * Arithmetic<Real>::epsilon());
    return sum;
  }
  // else compute directly the result
  return (a*log(double(a/b))+b-a);
}

} // namespace Funct

} // namespace STK
