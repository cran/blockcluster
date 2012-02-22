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
 * Purpose: implemen the interface class IGaussianModel .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IGaussianModel.cpp
 *  @brief In this file we implement the Interface class IGaussianModel.
 **/


#include "../include/STK_IGaussianModel.h"

#include "../../STatistiK/include/STK_Stat_UnivariateReal.h"

namespace STK
{

/* constructor */
IGaussianModel::IGaussianModel( Matrix const* p_data)
                              : RealModel(p_data)
                              , mean_(p_data_->rangeHo())
{ }

/* destructor */
IGaussianModel::~IGaussianModel() {}

/* compute the empirical mean */
void IGaussianModel::compMean()
{
  // resize mean
  mean_.resize(p_data_->rangeHo());
  // get dimensions
  const Integer first = p_data_->firstCol(), last = p_data_->lastCol();
  for (Integer j= first; j <= last; ++j)
  {
    mean_[j] = Stat::mean<Vector>(p_data_->col(j));
  }
}

/* compute the weighted empirical mean */
void IGaussianModel::compWeightedMean(Vector const& weights)
{
  // resize mean
  mean_.resize(p_data_->rangeHo());
  // get dimensions
  const Integer firstCol = p_data_->firstCol(), lastCol = p_data_->lastCol();
  for (Integer j= firstCol; j <= lastCol; ++j)
  {
    mean_[j] = Stat::mean<Vector>(p_data_->col(j), weights);
  }
}

} // namespace STK
