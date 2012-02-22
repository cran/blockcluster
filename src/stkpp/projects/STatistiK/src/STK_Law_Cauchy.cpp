/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2008  Serge Iovleff

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
 * Project: STatistiK
 * Purpose: Implementation of the Cauchy Distribution
 * Author:  Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Law_Cauchy.cpp
 *  @brief In this file we implement the Cauchy distribution.
 **/

#include "../include/STK_Law_Cauchy.h"

namespace STK
{

namespace Law
{

/* constructor */
Cauchy::Cauchy( Real const& location, Real const& scale)
              : ITUnivariate<Real>(String(_T("Cauchy")))
              , location_(location)
              , scale_(scale)
{
  // check parameters
  if (  !Arithmetic<Real>::isFinite(location)
     || !Arithmetic<Real>::isFinite(scale)
     || scale <= 0
     )
    throw domain_error("Cauchy::Cauchy(location, scale) "
                            "argument error");
}

/* destructor */
Cauchy::~Cauchy()
{ ;}

/*  Generate a pseudo Cauchy random variate. */
Real Cauchy::rand() const
{
  return location_ + scale_ * tan((double)(Const::_PI_ * generator.randUnif()));
}

/*
 *  Generate a pseudo Cauchy random variate with the specified parameters.
 *  (static)
 */
Real Cauchy::rand( Real const& location, Real const& scale)
{
  // check parameters
  if (  !Arithmetic<Real>::isFinite(location)
     || !Arithmetic<Real>::isFinite(scale)
     || scale <= 0
     )
    throw domain_error("Cauchy::Cauchy(location, scale) "
                            "argument error");
  return location + scale * tan(Const::_PI_ * generator.randUnif());
}

/*  Give the value of the pdf at x.*/
Real Cauchy::pdf( Real const& x) const
{
  // check NA value
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();

  // trivial case
  if (Arithmetic<Real>::isInfinite(x)) return 0.0;

  // general case
  Real y = (x - location_) / scale_;
  return 1. / (Const::_PI_ * scale_ * (1. + y * y));
}

/* Give the value of the log-pdf at x. */
Real Cauchy::lpdf( Real const& x) const
{
  // check NA value
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();

  // trivial case
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();

  // general case
  Real y = (x - location_) / scale_;
  return - log(Const::_PI_ * scale_ * (1. + y * y));
}

/* The cumulative distribution function at t.
 */
Real Cauchy::cdf( Real const& t) const
{
  // check NA value
  if (Arithmetic<Real>::isNA(t)) return Arithmetic<Real>::NA();

  // check parameter
  if (Arithmetic<Real>::isInfinite(t))
   return (t < 0.) ? 0.0 : 1.0;

  /* http://www.faqs.org/faqs/fr/maths/maths-faq-3/
   * arctan on [0, 1[:
   * if x<0 atan(x)= -atan(-x)
   * elseif x>1 atan(x)= Pi/2-atan(1/x).
   */
  if (abs(t) > 1)
  {
    double y = atan(1/t) / Const::_PI_;
    return (t > 0) ? (1 - y) : (-y);
  }
  else
    return 0.5 + atan(t) / Const::_PI_;
}
    
/*
 * The inverse cumulative distribution function at p.
 */
Real Cauchy::icdf( Real const& p) const
{
  // check NA value
  if (Arithmetic<Real>::isNA(p)) return Arithmetic<Real>::NA();

  // check parameter
  if ((p > 1.) || (p < 0.))
    throw domain_error("Cauchy::icdf(p) "
                            "argument error");
 // trivial cases
 if (p == 0.)  return -Arithmetic<Real>::infinity();
 if (p == 1.)  return  Arithmetic<Real>::infinity();

  // general case
  // tan(pi * (p - 1/2)) = -cot(pi * p) = -1/tan(pi * p) 
  return location_ - scale_ / tan(Const::_PI_ * p);
}

} // namespace Law

} // namespace STK
