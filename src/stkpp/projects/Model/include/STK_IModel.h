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
 * created on: 22 juil. 2011
 * Purpose: define the Interface base class IModel (Statistical Model).
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IModel.h
 *  @brief In this file we define the Interface base class IModel.
 **/

#ifndef STK_IMODEL_H
#define STK_IMODEL_H

#include <cmath>

#include "../../Sdk/include/STK_IRunner.h"
#include "../../STKernel/include/STK_Real.h"
#include "../../STKernel/include/STK_Integer.h"

namespace STK
{

/** @ingroup Model
 *  @brief Interface base class for all statistical models.
 *  A statistical model is defined with two elements
 *  - A data set where the number of samples is the number of rows
 *    and the number of variables is the number of column
 *  - A parameterized probability density
 **/
class IModel
{
  protected:
    /** Constructor.
     * @param nbSample number of sample of the model
     * @param nbVar number of variable of the model
     **/
    IModel(Integer const& nbSample =0, Integer const& nbVar =0);

    /** virtual destructor */
    virtual ~IModel();

  public:
    /** get the total available observations
     *  @return the total available observations
     **/
    inline Integer const& nbSample() const { return nbSample_;}

    /** get the total available observations
     *  @return the total available observations
     **/
    inline Integer const& nbVar() const { return nbVar_;}

    /** get the log of the total available observations
     *  @return Log of the total available observations
     **/
    inline Real logNbSample() const
    { return (nbSample_ <= 0) ? Arithmetic<Real>::NA(): log((double)nbSample_);}

    /** @brief get the log-likelihood
     *  @return The log-likelihood
     **/
    inline Real logLikelihood() const { return logLikelihood_;}

    /** @brief get the likelihood
     *  @return The likelihood
     **/
    inline Real likelihood() const
    { return (Arithmetic<Real>::isFinite(logLikelihood_)) ? exp((double)logLikelihood_) : 0.;}

    /** get the number of free parameters
     * @return the number of free parameter
     */
    inline Integer const& nbFreeParameter() const { return nbFreeParameter_;}

  protected:
    /** total available samples */
    Integer nbSample_;
    /** total available variables */
    Integer nbVar_;
    /** likelihood of the sample */
    Real logLikelihood_;
    /** number of free parameter of the model */
    Integer nbFreeParameter_;
    /** set the default value of the parameters of the model
     *  @param nbSample number of sample of the model
     *  @param nbVar number of variable of the model
     * */
    void setDefault( Integer const& nbSample =0, Integer const& nbVar =0);
};

} // namespace STK

#endif /* STK_IMODEL_H */
