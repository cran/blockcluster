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
 * Project:  stkpp::Algebra
 * Purpose:  Define 1D Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_LinAlgebra1D.h
 *  @brief In this file we implement Linear Algebra methods
 *  for Real one dimensional containers.
 **/
 
#ifndef STK_LINALGEBRA1D_H
#define STK_LINALGEBRA1D_H

#include "../../STKernel/include/STK_Real.h"
#include "../../STKernel/include/STK_Misc.h"

#include "../../Sdk/include/STK_ITContainer1D.h"

namespace STK
{

/** @ingroup Algebra
 *  @brief Sum the element of the container
 * 
 *  Sum of the element of the Container1D @c x
 *  \f[ s= \sum_{i=1}^n x_i \f]

 *  @param[in] x vector to treat
 *  @return the sum of the element of the container
 **/
template<class Container1D>
Real sum( ITContainer1D<Real, Container1D> const& x)
{
  // get dimensions
  const Integer first = x.first(), last = x.last();
  // compute the sum
  Real sum = 0.0;
  for (Integer i=first; i<=last; i++)
    sum += x[i];
  return (sum);
}

/** @ingroup Algebra
 *  @brief Weighted sum of the element of the container
 *
 *  Sum of the weighted elements of the Container1D x
 *  \f[ s= \sum_{i=1}^n w_i x_i \f]

 *  @param[in] x vector to treat
 *  @param w the weights to apply
 *  @return the weighted sum of the container
 **/
template<class Container1D1, class Container1D2>
Real weightedSum( ITContainer1D<Real, Container1D1> const& x
                , ITContainer1D<Real, Container1D2> const& w
                )
{
#ifdef STK_DEBUG
    if (!x.range().isIn(w.range()))
      throw runtime_error("In weightedSum(x, w) "
                               "x.range() not include in w.range()");
#endif
  // dimensions
  const Integer first = x.first(), last = x.last();
  // compute the weighted sum
  Real sum = 0.0;
  for (Integer i=first; i<=last; i++)
    sum += w[i] * x[i];
  return (sum);
}

/** @ingroup Algebra
 *  @brief Compute the infinite norm
 *
 *  Compute the maximal absolute value of the container @c x
 *  \f[ s= \max_{i=1}^n |x_i| \f]
 *  @param[in] x vector to treat
 *  @return the infinite norm f the container
 **/
template<class Container1D>
Real normInf( ITContainer1D<Real, Container1D> const& x)
{
  // get dimensions
  const Integer first = x.first(), last = x.last();
  // compute the maxmal value
  Real scale = 0.0;
  for (Integer i=first; i<=last; i++)
    scale = max(scale, abs(x[i]));
  return (scale);
}

/** @ingroup Algebra
 *  @brief Compute the weighted infinite norm
 *
 *  Compute the maximal absolute weighted value of the container @c x
 *  \f[ s= \max_{i=1}^n |w_i x_i| \f]
 *  @param[in] x vector to treat
 *  @param w the weights to apply
 *  @return the weighted infinite norm
 **/
template<class Container1D1, class Container1D2>
Real weightedNormInf( ITContainer1D<Real, Container1D1> const& x
                    , ITContainer1D<Real, Container1D2> const& w
                    )
{
#ifdef STK_DEBUG
    if (!x.range().isIn(w.range()))
      throw runtime_error("In weightedNormInf(x, w) "
                               "x.range() not include in w.range()");
#endif
  // get dimensions
  const Integer first = x.first(), last = x.last();
  // compute weighted norm inf
  Real scale = 0.0;
  for (Integer i=first; i<=last; i++)
    scale = max(scale, abs(w[i]*x[i]));
  return (scale);
}

/** @ingroup Algebra
 *  @brief compute the norm two
 *
 *  Compute the norm of the container @c x avoiding overflow
 *  \f[ \|x\| = \sqrt{\sum_{i=1}^n x^2_i } \f]
 *  @param[in] x vector to treat
 *  @return the norm two of the container
 **/
template<class Container1D>
Real normTwo( ITContainer1D<Real, Container1D> const& x)
{
  // compute the maximal value of x
  Real scale =normInf(x), norm =0.0;
  if (scale)
  {
    // get dimensions
    const Integer first = x.first(), last = x.last();
    // sum squared normalized values
    for (Integer i = first; i<=last; i++)
    {
      const Real aux = x[i]/scale;
      norm += aux * aux;
    }
  }
  // rescale sum
  return (Real(sqrt(double(norm)))*scale);
}

/** @ingroup Algebra
 *  @brief compute the weighted norm two
 * 
 *  Compute the weighted norm of the container @c x avoiding overflow
 *  \f[ \|x\| = \sqrt{\sum_{i=1}^n w_i x^2_i } \f]
 *  @param[in] x vector to treat
 *  @param w the weights to apply
 *  @return the weighted two norm
 **/
template<class Container1D1, class Container1D2>
Real weightedNormTwo( ITContainer1D<Real, Container1D1> const& x
                    , ITContainer1D<Real, Container1D2> const& w
                    )
{
#ifdef STK_DEBUG
    if (!x.range().isIn(w.range()))
      throw runtime_error("In weightedNormTwo(x, w) "
                               "x.range() not include in w.range()");
#endif
  // compute the maximal value of x
  Real scale = weightedNormInf(x, w), norm2 =0.0;
  if (scale)
  {
    // get dimensions
    const Integer first = x.first(), last = x.last();
    // compute norm2
    for (Integer i = first; i<=last; i++)
    {
      const Real aux = (w[i]*x[i])/scale;
      norm2 += aux * aux;
    }
  }
  // rescale sum
  return (Real(sqrt(double(norm2)))*scale);
}

/** @ingroup Algebra
 *  @brief Compute the squared norm two
 *
 *  Compute the square norm of the Container1D x avoiding overflow
 *  \f[ \|x\|^2 = \sum_{i=1}^n x^2_i \f]
 *
 *  @param[in] x vector to treat
 **/
template<class Container1D>
Real normTwo2( ITContainer1D<Real, Container1D> const& x)
{
  Real scale =normInf(x), norm =0.0;
  if (scale)
  {
    // get dimensions
    const Integer first = x.first(), last = x.last();
    // sum squared normalized values
    for (Integer i = first; i<=last; i++)
    {
      Real aux = x[i]/scale;
      norm += aux * aux;
    }
  }
  // scale result
  return (norm*scale*scale);
}

/** @ingroup Algebra
 *  @brief Compute the squared norm two
 * 
 *  Compute the weighted square norm two of the container x avoiding overflow
 *  \f[ \|x\|^2 = \sum_{i=1}^n w_i x^2_i \f]
 * 
 *  @param[in] x vector to treat
 *  @param w the weights to apply
 *  @return the weighted square norm two
 **/
template<class Container1D1, class Container1D2>
Real weightedNormTwo2( ITContainer1D<Real, Container1D1> const& x
                     , ITContainer1D<Real, Container1D2> const& w
                     )
{
#ifdef STK_DEBUG
  if (!x.range().isIn(w.range()))
    throw runtime_error("In weightedNormTwo2(x, w) "
                             "x.range() not include in w.range()");
#endif
  Real scale =weightedNormInf(x, w), norm =0.0;
  if (scale)
  {
    // get dimensions
    const Integer first = x.first(), last = x.last();
    // compute norm2
    for (Integer i = first; i <= last; i++)
    {
      const Real aux = x[i]/scale;
      norm += w[i] * aux * aux;
    }
  }
  // scale result
  return (norm*scale*scale);
}

/** @ingroup Algebra
 *  @brief Compute the dot product.
 * 
 *  Dot product of the vector x and the vector y: d = <x, y>.
 *  \f[ <x,y> = \sum_{i=1}^n x_i y_i \f]
 *  The common range of the vectors is used. The value outside the range
 *  are thus interpreted as zero.
 * 
 *  @param[in] x first vector
 *  @param[in] y second vector
 *  @return the dot product of the two vectors
 **/
template<class Container1D1, class Container1D2>
Real dot( ITContainer1D<Real, Container1D1> const& x
        , ITContainer1D<Real, Container1D2> const& y
        )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  // compute the sum product
  Real sum=0.0;
  Integer i;
  for (i = first; i<last; i+=2)
    sum += x[i] * y[i] + x[i+1] * y[i+1];
  // check if the number of element is odd
  if (i==last) sum+=x[last]*y[last];
  return (sum);
}

/** @ingroup Algebra
 *  @brief Compute the dot product.
 *
 *  Weighted dot product of the vector x and the vector y:
 *  \f[ <x,y> = \sum_{i=1}^n w_i x_i y_i. \f]
 *  The common range of the vectors is used. The value outside the range
 *  are thus interpreted as zero.
 *
 *  @param[in] x first vector
 *  @param[in] y second vector
 *  @param[in] w the weights to apply
 *  @return the weighted dot product of the two vectors
 **/
template<class Container1D1, class Container1D2, class Container1D3>
Real weightedDot( ITContainer1D<Real, Container1D1> const& x
                , ITContainer1D<Real, Container1D2> const& y
                , ITContainer1D<Real, Container1D3> const& w
                )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
#ifdef STK_DEBUG
  if (!Range(first,last).isIn(w.range()))
    throw runtime_error("In weightedDot(x, w) "
                             "Range(first,last) not include in w.range()");
#endif
  // compute the sum product
  Real sum=0.0;
  Integer i;
  for (i = first; i<last; i+=2)
    sum += w[i]*x[i] * y[i] + w[i+1]*x[i+1] * y[i+1];
  // check if last is odd
  if (i == last) sum += w[last]*x[last]*y[last];
  return (sum);
}

/** @ingroup Algebra
 *  @brief Compute the distance between two vectors.
 * 
 *  Compute the Euclidian distance between x and y without overflow.
 * 
 *  \f[ d(x,y) = || x - y|| = \sqrt{\sum_{i=1}^n |x_i - y_i|^2}. \f]
 *  The common range of the vectors is used. The value outside the range
 *  are thus interpreted as equal.
 * 
 *  @param[in] x first vector
 *  @param[in] y second vector
 *  @return the Euclidian distance between x and y
 **/
template<class Container1D1, class Container1D2>
Real dist( ITContainer1D<Real, Container1D1> const& x
         , ITContainer1D<Real, Container1D2> const& y
         )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  // compute the maximal difference
  Real scale = 0.;
  for (Integer i = first; i<=last; i++)
    scale = max(scale, abs(x[i] - y[i]));
  // Compute the norm
  Real norm2 = 0.;
  if (scale)
  { // comp the norm^2
    for (Integer i = first; i<=last; i++)
    {
      const Real aux = (x[i]-y[i])/scale;
      norm2 += aux * aux;
    }
  }
  // rescale sum
  return (Real(sqrt(double(norm2)))*scale);
}

/** @ingroup Algebra
 *  @brief Compute the weighted distance between two vectors.
 *
 *  Compute the weighted Euclidian distance between x and y without overflow.
 *
 *  \f[ d(x,y) = || x - y|| = \sqrt{\sum_{i=1}^n w_i |x_i - y_i|^2}. \f]
 *  The common range of the vectors is used. The value outside the range
 *  are thus interpreted as equal.
 *
 *  @param[in] x first vector
 *  @param[in] y second vector
 *  @param[in] w the weight of the data
 *  @return the weighted Euclidian distance between x and y
 **/
template<class Container1D1, class Container1D2, class Container1D3>
Real weightedDist( ITContainer1D<Real, Container1D1> const& x
                 , ITContainer1D<Real, Container1D2> const& y
                 , ITContainer1D<Real, Container1D3> const& w
                 )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
#ifdef STK_DEBUG
  if (!Range(first,last).isIn(w.range()))
    throw runtime_error("In weightedDist(x, w) "
                             "Range(first,last) not include in w.range()");
#endif
  // compute the maximal difference
  Real scale = 0., norm2= 0.;
  for (Integer i = first; i<=last; i++)
    scale = max(scale, abs(w[i]*(x[i] - y[i])));
  // Compute the norm
  if (scale)
  { // comp the norm^2
    for (Integer i = first; i<=last; i++)
    {
      const Real aux = (x[i]-y[i])/scale;
      norm2 += w[i]*aux * aux;
    }
  }
  // rescale sum
  return (Real(sqrt(double(norm2)))*scale);
}

/** @ingroup Algebra
 *  @brief add two vectors.
 *
 *  Add two vectors and store the result in the first vector.
 *
 *  \f[ y = y + x. \f]
 *
 *  @param y the vector to translate
 *  @param x the vector to add
 **/
template<class Container1D1, class Container1D2>
void add(ITContainer1D<Real, Container1D1> const& x, ITContainer1D<Real, Container1D2>& y)
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  for (Integer i = first; i<=last; i++)
  {
    y[i] += x[i];
  }
}

/** @ingroup Algebra
 *  @brief add two vectors.
 *
 *  Add two vectors and store the result in a third vector.
 *
 *  \f[ z = y + x. \f]
 *
 *  @param y first the vector to add
 *  @param x second vector to add
 *  @param z the result
 **/
template<class Container1D1, class Container1D2, class Container1D3>
void add( ITContainer1D<Real, Container1D1> const& x
        , ITContainer1D<Real, Container1D2> const& y
        , ITContainer1D<Real, Container1D3>& z
        )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  z.resize(Range(first, last));
  for (Integer i = first; i<=last; i++)
  {
    z[i] = x[i] + y[i];
  }
}

/** @ingroup Algebra
 *  @brief add two vectors.
 *
 *  soustract two vectors and store the result in the first vector.
 *
 *  \f[ y = y - x. \f]
 *
 *  @param y the vector to translate
 *  @param x the vector to add
 **/
template<class Container1D1, class Container1D2>
void diff(ITContainer1D<Real, Container1D1> const& x, ITContainer1D<Real, Container1D2>& y)
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  for (Integer i = first; i<=last; i++)
  {
    y[i] -= x[i];
  }
}

/** @ingroup Algebra
 *  @brief add two vectors.
 *
 *  Soustract two vectors and store the result in a third vector.
 *
 *  \f[ z = x - y. \f]
 *
 *  @param y first the vector to add
 *  @param x second vector to add
 *  @param z the result
 **/
template<class Container1D1, class Container1D2, class Container1D3>
void diff( ITContainer1D<Real, Container1D1> const& x
         , ITContainer1D<Real, Container1D2> const& y
         , ITContainer1D<Real, Container1D3>& z
         )
{
  // compute the valid range
  const Integer first = max(x.first(), y.first()) , last = min(x.last(), y.last());
  z.resize(Range(first, last));
  for (Integer i = first; i<=last; i++)
  {
    z[i] = x[i] - y[i];
  }
}


} // Namespace STK

#endif // STK_LINALGEBRA1D_H
