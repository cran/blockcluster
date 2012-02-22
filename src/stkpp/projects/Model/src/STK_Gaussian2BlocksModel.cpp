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

#include "../include/STK_Gaussian2BlocksModel.h"

#include "../../STatistiK/include/STK_Stat_UnivariateReal.h"
#include "../../STatistiK/include/STK_Stat_BivariateRealReal.h"
#include "../../STatistiK/include/STK_Law_MultivariateNormal.h"

namespace STK
{

Gaussian2BlocksModel::Gaussian2BlocksModel( const Matrix *p_data)
                                                  : GaussianModel(p_data)
                                                  , dim_(p_data_->sizeHo())
{ }


/* destructor */
Gaussian2BlocksModel::~Gaussian2BlocksModel()
{ }


/* compute the empirical covariance matrix. */
void Gaussian2BlocksModel::compCovariance()
{
  // resize mean
  cov_.resize(p_data_->rangeHo());
  cov_ = 0.;
  // get dimensions for the first block
  const Integer first1 = p_data_->firstCol(), last1 = min(p_data_->lastCol(), first1+dim_-1);
  for (Integer i= first1; i <= last1; ++i)
  {
    cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), mean_[i]);
    for (Integer j= first1; j < i; ++j)
    {
      cov_(i, j) = Stat::covarianceWithFixedMean<Vector>(p_data_->col(i), p_data_->col(j), mean_[i], mean_[j]);
      cov_(j, i) = cov_(i,j);
    }
  }
  // get dimensions for the second block
  const Integer first2 = last1+1, last2 = p_data_->lastCol(), size2 = last2 - last1 ;
  if (size2)
  {
    // compute variance of each column
    for (Integer i= first2; i <= last2; ++i)
    { cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), mean_[i]);}
    variance2_ = trace(MatrixSquare(cov_, Range(first2, last2)))/(Real)size2;
    // compute variance of each column
    for (Integer i= first2; i <= last2; ++i)
    { cov_(i, i) = variance2_;}
  }
}

/* compute the empirical covariance matrix. */
void Gaussian2BlocksModel::compWeightedCovariance(Vector const& weights)
{
  // resize mean
  cov_.resize(p_data_->rangeHo());
  // get dimensions for the first block
  const Integer first1 = p_data_->firstCol(), last1 = min(p_data_->lastCol(), dim_);
  for (Integer i= first1; i <= last1; ++i)
  {
    cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), weights, mean_[i]);
    for (Integer j= first1; j < i; ++j)
    {
      cov_(i, j) = 0.;
      cov_(j, i) = 0.;
    }
  }
  // get dimensions for the second block
  const Integer first2 = last1+1, last2 = p_data_->lastCol(), size2 = last2 - last1 ;
  if (size2)
  {
    // compute variance of each column
    for (Integer i= first2; i <= last2; ++i)
    { cov_(i, i) = Stat::varianceWithFixedMean<Vector>(p_data_->col(i), weights, mean_[i]);}
    variance2_ = trace(MatrixSquare(cov_, Range(first2, last2)))/(Real)size2;
    for (Integer i= first2; i <= last2; ++i)
    { cov_(i, i) = variance2_;}
  }
}

} // namespace STK
