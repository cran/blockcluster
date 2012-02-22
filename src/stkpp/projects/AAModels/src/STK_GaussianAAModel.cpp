  /*--------------------------------------------------------------------*/
/*     Copyright (C) 2004  Serge Iovleff

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
 * Project:  stkpp::AAModels
 * Purpose:  Interface base class for AA models.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_GaussianAAModel.cpp
 *  @brief In this file we implement the class GaussianAAModel for the
 *  auto-associative models.
 **/

#include "../include/STK_GaussianAAModel.h"

#include "../../Algebra/include/STK_LinAlgebra1D.h"
#include "../../Algebra/include/STK_TExpAlgebra.h"

#include "../../Reduct/include/STK_IReduct.h"
#include "../../Regress/include/STK_IRegression.h"

#include "../../STatistiK/include/STK_Stat_Transform.h"
#include "../../STatistiK/include/STK_Stat_MultivariateReal.h"

#include "../../STatistiK/include/STK_Stat_MultivariateReal.h"
#include "../../STatistiK/include/STK_Law_MultivariateNormal.h"

#ifdef STK_VERBOSE
#include "../../Arrays/include/STK_Display2D.h"
#endif

namespace STK
{

// constructor
GaussianAAModel::GaussianAAModel( Matrix& workData)
                                : IAAModel(workData)
                                , IModel(workData.sizeVe(), workData.sizeHo())
{ }

// destructor
GaussianAAModel::~GaussianAAModel()
{ }

/* update the container when the data set is modified. **/
void GaussianAAModel::setWorkData(Matrix& workData)
{
  // update data set and flags for the IAAModel part
  IAAModel::setWorkData(workData);
  // set dimensions to new size for the IModel part
  setDefault(workData.sizeVe(), workData.sizeHo());
}

/* compute the log-likelihood of the model */
void GaussianAAModel::computeLogLikelihood()
{
#ifdef STK_VERBOSE
  stk_cout << "GaussianAAModel::computeLogLikelihood().\n";
#endif
  // get number of free parameters
  nbFreeParameter_ = p_regressor_->nbParameter();
  // range of the column to use
  Range cols = Range(dim())+(p_reduced_->firstCol()-1);
  // create a reference with the first columns of the reduced data
  Matrix reducedData(*p_reduced_, p_reduced_->rangeVe(), cols);
  // create a reference with the first columns of the reduced data
  MatrixSquare reducedCovariance(covProjected(), cols);

  // compute first part of the log-likehhod
  Point mean(reducedCovariance.range(), 0.);
  Law::MultivariateNormal<Point> normalP(mean, reducedCovariance);
  Real loglikehood1 = normalP.logLikelihood(reducedData);

  // compute second part of the log-likehhod
  Real loglikehood2 = (Const::_LNSQRT2PI_+log((double)varResiduals()))*(dim() - nbVar_);
  const Integer firstSample = p_residuals_->firstRow(), lastSample= p_residuals_->lastRow();
  for (Integer i=firstSample; i<=lastSample; i++)
  {
    loglikehood2 -= normTwo2(p_residuals_->row(i))/(2.*varResiduals());
  }
  logLikelihood_ = loglikehood1 + loglikehood2;

#ifdef STK_VERBOSE
  stk_cout << "GaussianAAModel::computeLogLikelihood() done.\n";
  stk_cout << _T("loglikehood1 = ") << loglikehood1 << _T("\n");
  stk_cout << _T("loglikehood2 = ") << loglikehood2 << _T("\n");
#endif
}


} // namespace STK

