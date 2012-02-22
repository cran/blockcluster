/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2010  Serge Iovleff

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
 * created on: 25 juin 2010
 * Purpose:  implement the BSplineCoefficients class.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_BSplineCoefficients.cpp
 *  @brief In this file we implement the BSplineCoefficients class.
 **/

#include "../include/STK_BSplineCoefficients.h"
#include "../../DManager/include/STK_HeapSort.h"
#include "../../STKernel/include/STK_String_Util.h"

#ifdef STK_VERBOSE
#include "../../Arrays/include/STK_Display2D.h"
#endif


namespace STK
{
/* convert a String to a TypeReduction.
 *  @param type the type of reduction we want to define
 *  @return the TypeReduction represented by the String @c type. if the string
 *  does not match any known name, the @c unknown_ type is returned.
 **/
BSplineCoefficients::KnotsPosition BSplineCoefficients::StringToKnotsPosition( String const& type)
{
  if (toUpperString(type) == toUpperString(_T("uniform")))  return uniform_;
  if (toUpperString(type) == toUpperString(_T("periodic"))) return periodic_;
  if (toUpperString(type) == toUpperString(_T("density"))) return density_;
  return unknown_;
}

/* convert a TypeReduction to a String.
 *  @param type the type of reduction we want to convert
 *  @return the string associated to this type.
 **/
String BSplineCoefficients::KnotsPositionToString( KnotsPosition const& type)
{
  if (type == uniform_)  return String(_T("uniform"));
  if (type == periodic_) return String(_T("periodic"));
  if (type == density_) return String(_T("density"));
  return String(_T("unknown"));
}

/* constructor */
BSplineCoefficients::BSplineCoefficients( Vector const* p_data
                                        , Integer const& nbControlPoints
                                        , Integer const& degree
                                        , const KnotsPosition& position
                                        )
                                        : p_data_(p_data)
                                        , nbKnots_(nbControlPoints + degree +1)
                                        , lastKnot_(nbKnots_-1)
                                        , nbControlPoints_(nbControlPoints)
                                        , lastControlPoint_(nbControlPoints_-1)
                                        , degree_(degree)
                                        , position_(position)
                                        , knots_(Range(0,lastKnot_))
                                        , coefficients_(p_data->range(), Range(0,lastControlPoint_), 0.0)
                                        , minValue_( Arithmetic<Real>::max())
                                        , maxValue_(-Arithmetic<Real>::max())
{
}

/* constructor */
BSplineCoefficients::BSplineCoefficients( Vector const& data
                                        , Integer const& nbControlPoints
                                        , Integer const& degree
                                        , const KnotsPosition& position
                                        )
                                        : p_data_(&data)
                                        , nbKnots_(nbControlPoints + degree +1)
                                        , lastKnot_(nbKnots_-1)
                                        , nbControlPoints_(nbControlPoints)
                                        , lastControlPoint_(nbControlPoints_-1)
                                        , degree_(degree)
                                        , position_(position)
                                        , knots_(Range(0,lastKnot_))
                                        , coefficients_(p_data_->range(), Range(0,lastControlPoint_), 0.0)
                                        , minValue_( Arithmetic<Real>::max())
                                        , maxValue_(-Arithmetic<Real>::max())
{}

// destructor
BSplineCoefficients::~BSplineCoefficients() {}

/** run the computations. */
void BSplineCoefficients::run()
{
  computeKnots();
  computeCoefficients();
}


/*
 *  run the computations for the given value.
 *  @param p_data the input data values
 **/
void BSplineCoefficients::setData( Vector const* p_data
                                 , Integer const& nbControlPoints
                                 , Integer const& degree
                                 , KnotsPosition const& position
                                 )
{
  // set data
  p_data_ = p_data;
  // resize coeficients
  coefficients_.resize(p_data->range(), Range(0, lastControlPoint_));
  // initialize array of coefficient
  coefficients_ = 0.0;
}

/* compute the knots of the B-Spline curves.*/
void BSplineCoefficients::computeKnots()
{
  // get dimensions
  Integer first = p_data_->first(), last = p_data_->last();
  // compute min value
  for (Integer i=first; i<= last; i++)
  {
    minValue_ = min(minValue_, (*p_data_)[i]);
    maxValue_ = max(maxValue_, (*p_data_)[i]);
  }
  // if all value are equals, all the knots are equals to this value
  if (minValue_ == maxValue_)
  {
    knots_ = minValue_;
    return;
  }
  // set knots values
  switch (position_)
  {
    // uniform position
    case uniform_:
      computeUniformKnots();
      break;
    // periodic position
    case periodic_:
      computePeriodicKnots();
      break;
    // density position
    case density_:
      computeDensityKnots();
      break;
      // periodic position
    case unknown_:
      // check if there exists data
      throw runtime_error("Error In BSplineCoefficients::computeKnots():"
                               " unknowns positions");
      break;
  }
  // shift knots
  Real range = (maxValue_ - minValue_);
  for (Integer k = 0; k <= lastKnot_; k++)
    knots_[k] = minValue_ + range * knots_[k];
}

/* compute the position of the uniform knots.*/
void BSplineCoefficients::computeUniformKnots()
{
  // compute step
  Real step = 1.0/(nbControlPoints_ - degree_);
  // set internal knots
  const Integer first = degree_ + 1, last = lastControlPoint_;
  for (Integer k = first, j = 1; k <= last; j++, k++)
    knots_[k] = j * step;
  // set external knots
  for ( Integer k=0, j = last+1; k < first; j++, k++)
  {
    knots_[k] = 0;
    knots_[j] = 1;
  }
}
/* compute the position of the periodic knots.*/
void BSplineCoefficients::computePeriodicKnots()
{
  // compute step
  Real step = 1.0/(nbControlPoints_ - degree_);
  // set knots
  for (Integer k = 0, j = -degree_; k <= lastKnot_; j++, k++)
    knots_[k] = j * step;
;
}
/* compute the position of the density knots. */
void BSplineCoefficients::computeDensityKnots()
{
  // sorted data
  Vector xtri;
  // sort the data
  heapSort<Real, Vector>(*p_data_, xtri);

  // compute step
  Real step = xtri.size()/(Real)lastKnot_;
  Integer first = xtri.first(), last = xtri.last();
  // set knots
  for (Integer k = 0; k < lastKnot_; k++)
  {
    Integer cell = first + Integer(k* step);
    knots_[k] = (xtri[cell] + xtri[cell+1])/2.;
  }
  // set last knots
  knots_[lastKnot_] = (xtri[last-1] + xtri[last])/2.;
}

/* Compute the coefficients of the B-Spline curves.*/
void BSplineCoefficients::computeCoefficients()
{
#ifdef STK_VERBOSE
  stk_cout << _T("BSplineCoefficients::computeCoefficients()\n");
#endif
  // get dimensions
  Integer first = p_data_->first(), last = p_data_->last();
  // compute the coefficients
  for (Integer i=first; i<= last; i++)
  {
    computeCoefficientsRow(i);
  }
#ifdef STK_VERBOSE
  stk_cout << _T("BSplineCoefficients::computeCoefficients() done\n");
#endif
}

/* Compute a row of the coefficients
 * @param irow index of the row
 **/
void BSplineCoefficients::computeCoefficientsRow(Integer const& irow)
{
  // get current value
  const Real value = (*p_data_)[irow];
  // value outside the range of the knots case
  if (value <= minValue_)
  {
    coefficients_(irow, 0) = 1.0;
    return;
  }
  if (value >= maxValue_)
  {
    coefficients_(irow, lastControlPoint_) = 1.0;
    return;
  }
  // find interval
  Integer k, k1;
  for (k=0, k1=1; k<lastControlPoint_; k++, k1++)
  {
    if (value < knots_[k1]) break;
  }
  // begin recursion
  coefficients_(irow, k) = 1.0;
  for (Integer d=1; d<=degree_; d++)
  {
    // right (south-west corner) term only
    coefficients_(irow, k-d) = ( (knots_[k1] - value)/(knots_[k1] - knots_[k1-d]) )
                               * coefficients_(irow, k1-d);
    // compute internal terms
    for (Integer i = k1-d; i<k; i++)
    {
      const Real knots_i = knots_[i], knots_id1 = knots_[i+d+1];
      coefficients_(irow, i) = ( (value - knots_i)/(knots_[i+d] - knots_i) )
                               * coefficients_(irow, i)
                             + ( (knots_id1 - value)/(knots_id1 - knots_[i+1]) )
                               * coefficients_(irow, i+1);
    }
    // left (north-west corner) term only
    coefficients_(irow, k) *= (value - knots_[k])/(knots_[k+d] - knots_[k]);
  }
}

} // namespace STK
