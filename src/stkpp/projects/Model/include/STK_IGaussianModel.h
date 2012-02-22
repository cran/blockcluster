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
 * Project:  stkpp::Model
 * created on: 13 août 2011
 * Purpose:  Create a gaussian statistical model.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IGaussianModel.h
 *  @brief In this file we define the Interface class IGaussianModel.
 **/

#ifndef STK_IGAUSSIANMODEL_H
#define STK_IGAUSSIANMODEL_H

#include "STK_ITModel.h"

#include "../../Arrays/include/STK_MatrixSquare.h"

namespace STK
{

/** @ingroup Model
 *  @brief Compute the the maximum likelihood estimates of a Gaussian statistical
 *  model.
 *
 *  This is an interface class that can be derived in order to impose various
 *  constraint on the covariance matrix.
 **/

class IGaussianModel : public ITModel<Real, Matrix>
{
  typedef ITModel<Real, Matrix> RealModel;

  protected:
    /** constructor.
     * @param p_data pointer on the data set
     **/
    IGaussianModel( Matrix const* p_data);

  public:
    /** destructor. */
    virtual ~IGaussianModel();

    /** get the empirical mean.
     * @return the empirical mean
     **/
    inline Point const& mean() const { return mean_;}

  protected:
    /** compute the empirical mean */
    virtual void compMean();
    /** compute the empirical weighted mean
     * @param weights the weights of the samples
     **/
    virtual void compWeightedMean(Vector const& weights);

    /** compute the empirical covariance matrix. */
    virtual void compCovariance() =0;
    /** compute the empirical weighted covariance matrix.
     * @param weights the weights of the samples
     **/
    virtual void compWeightedCovariance(Vector const& weights) =0;
    /** Vector of the empirical means */
    Point mean_;
};

} // namespace STK

#endif /* STK_IGAUSSIANMODEL_H */
