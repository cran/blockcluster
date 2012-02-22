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
 * Project:  stkpp::aam
 * created on: 27 sept. 2010
 * Purpose:  .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_AdditiveBSplineRegressionPage.cpp
 *  @brief In this file we implement the AdditiveBSplineRegressionPage class
 **/

#include "../include/STK_AdditiveBSplineRegressionPage.h"

using namespace STK;

namespace STK
{

/** @brief some errors messages.
 **/
static const String ERRORCODES[] =
{
    _T("AdditiveBSplineRegressionPage::UNKNOWN_ERROR\n")
  , _T("AdditiveBSplineRegressionPage::INVALID_OPTION\n")
};

/* constructor. */
AdditiveBSplineRegressionPage::AdditiveBSplineRegressionPage(Integer const& level)
                            : IPage(_T("AdditiveBSplineRegression"), level, true)

{
  // reserve for the 4 options
  options_.reserve(3);
  // first option : standardize or no the data set
  options_.push_back(Option(_T("nbControlPoints"), Option::integer_, true));
  options_.back().setValue(_T("7"));
  // degree f the spline
  options_.push_back(Option(_T("degree"), Option::integer_, true));
  options_.back().setValue(_T("2"));
  // first option : type of model
  options_.push_back(Option(_T("knotsPositions"), Option::string_, true));
  options_.back().setValue(_T("uniform"));
}

/* copy constructor */
AdditiveBSplineRegressionPage::AdditiveBSplineRegressionPage( AdditiveBSplineRegressionPage const& page)
                                    : IPage(page)
                                    , nbControlPoints_(page.nbControlPoints_)
                                    , degree_(page.degree_)
                                    , position_(page.position_)
{ }

/* destructor. */
AdditiveBSplineRegressionPage::~AdditiveBSplineRegressionPage()
{ }


bool AdditiveBSplineRegressionPage::validate()
{
  // first option
  nbControlPoints_ = options_[0].get(Integer());
  if ( (nbControlPoints_ <= 0))
  {
    msg_error_ = ERRORCODES[1];
    msg_error_ += "nbControlPoints must be strictly positive\n";
    return false;
  }
  // second option
  degree_ = options_[1].get(Integer());
  if ( (degree_ <= 0))
  {
    msg_error_ = ERRORCODES[1];
    msg_error_ += "degree must be strictly positive\n";
    return false;
  }

  // third option
  position_ = BSplineCoefficients::StringToKnotsPosition(options_[2].get(String()));
  if (position_ == BSplineCoefficients::unknown_ )
  {
    msg_error_ = ERRORCODES[1];
    msg_error_ += "knots position is unknown.\n";
    return false;
  }
  // validation done
  return true;
}

} // namespace STK
