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
 * Project:  stkpp::
 * created on: 13 août 2011
 * Purpose:  .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_GaussianModel.cpp
 *  @brief In this file we implement the GaussianModel class.
 **/

#include "../include/STK_GaussianModel.h"

#include "../../STatistiK/include/STK_Stat_BivariateRealReal.h"

#include "../../Algebra/include/STK_TExpAlgebra.h"
#include "../../STatistiK/include/STK_Law_MultivariateNormal.h"

namespace STK
{

/* constructor */
GaussianModel::GaussianModel( Matrix const* p_data)
                            : IGaussianModel(p_data)
                            , cov_(p_data_->rangeHo())
{
  nbFreeParameter_ = nbVar_ + (nbVar_* (nbVar_-1))/2;
}

/* destructor */
GaussianModel::~GaussianModel() {}

/* implementation of the Gaussian statistical model
 * @return @c true if no error occur and @c false otherwise.
 */
bool GaussianModel::run()
{
  // compute the mean
  compMean();
  // compute the covariance matrix
  compCovariance();
  // create p_law_ (will be deleted in base class)
  // update gaussian law (will be deleted in base class)
  if (!p_law_) p_law_ = new Law::MultivariateNormal<Point>(mean_, cov_);
  else         p_law_->update();
  // compute log likelihood of the gaussian law
  logLikelihood_ = static_cast<Law::MultivariateNormal<Point>* >(p_law_)->logLikelihood(*p_data_ );
  // everything ok
  return true;
}

/* implementation of the weighted Gaussian statistical model
 * @param weights the weights of the samples
 * @return @c true if no error occur and @c false otherwise.
 */
bool GaussianModel::run(Vector const& weights)
{
  // compute the mean
  compWeightedMean(weights);
  // compute the covariance matrix
  compWeightedCovariance(weights);
  // create p_law_ (will be deleted in base class)
  // update gaussian law (will be deleted in base class)
  if (!p_law_) p_law_ = new Law::MultivariateNormal<Point>(mean_, cov_);
  else         p_law_->update();
  // compute log likelihood of the gaussian law
  logLikelihood_ = static_cast<Law::MultivariateNormal<Point>* >(p_law_)->logLikelihood(*p_data_ );
  // everything ok
  return true;
}

/* compute the empirical covariance matrix. */
void GaussianModel::compCovariance()
{
  // resize mean
  cov_.resize(p_data_->rangeHo());
  // get dimensions
  const Integer first = p_data_->firstCol(), last = p_data_->lastCol();
  for (Integer i= first; i <= last; ++i)
  {
    cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), mean_[i]);
    for (Integer j= first; j < i; ++j)
    {
      cov_(i, j) = Stat::covarianceWithFixedMean<Vector>(p_data_->col(i), p_data_->col(j), mean_[i], mean_[j]);
      cov_(j, i) = cov_(i,j);
    }
  }
}

/* compute the empirical covariance matrix. */
void GaussianModel::compWeightedCovariance(Vector const& weights)
{
  // resize mean
  cov_.resize(p_data_->rangeHo());
  // get dimensions
  const Integer first = p_data_->firstCol(), last = p_data_->lastCol();
  for (Integer i= first; i <= last; ++i)
  {
    cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), weights, mean_[i]);
    for (Integer j= first; j < i; ++j)
    {
      cov_(i, j) = Stat::covarianceWithFixedMean<Vector>(p_data_->col(i), p_data_->col(j), weights, mean_[i], mean_[j]);
      cov_(j, i) = cov_(i,j);
    }
  }
}

} // namespace STK
