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
 * Purpose:  implementation of the g0 function
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
 *  Compute the function:
 * \f[
 *     g_0(x) = x\log(x)+ 1 - x
 * \f]
 * with good accuracy for x near 1. using the Taylor serie
 * of
 * \f[ \log\left( \frac{a-x}{a+x}\right) \f]
 * x should be >= 0, no check is done.
 * @see http://www.herine.net/stat/software/dbinom.html
 * @param[in] x of the partial deviance
 **/
Real g0(Real const& x)
{
  // special values
  if (x == 0.) return 1.;
  if (x == 1.) return 0.;
  // general case
  Real v = (x-1.)/(x+1.);
  // if 7/9 < x < 9/7  use Taylor serie of log((1+v)/(1-v))
  if (v < 0.125)
  {
    Real sum  = (x-1)*v;
    Real ej   = 2*x*v;
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
  // else compute directly
  return x*log(double(x))-x + 1.;
}

} // namespace Funct

} // namespace STK
