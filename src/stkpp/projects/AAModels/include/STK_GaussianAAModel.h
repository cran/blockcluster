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

/** @file STK_GaussianAAModel.h
 *  @brief In this file we declare the class GaussianAAModel for the
 *  auto-associative models.
 **/

#ifndef STK_GAUSSIANAAMODEL_H
#define STK_GAUSSIANAAMODEL_H

#include "STK_IAAModel.h"
#include "../../Model/include/STK_IModel.h"

namespace STK
{
/** @ingroup AAModels
 *  @brief Gaussian AutoAssociative models.
 *  A Gaussian Auto-Associative model is a p-dimensional vector \f$\mathbf{y}\f$
 *  with projection function \f$\mathbf{P}\f$, if it can be written
 *  \f[
 *    \mathbf{y} = \mathbf{Q}
 *    \left(
 *    \begin{pmatrix}
 *      x_1 \\ \vdots \\
 *      x_d \\
 *      \tilde{r}_{d+1}(\mathbf{x}) \\ \vdots \\
 *      \tilde{r}_p(\mathbf{x})
 *    \end{pmatrix}
 *     + \tilde{\varepsilon}
 *    \right)
 *    = R(\mathbf{x})+ \varepsilon,
 *  \f]
 *  where the \f$\tilde{r}_j(\mathbf{x})\f$, \f$ d+1 \leq j \leq p\f$, are arbitrary
 *  real functions from \f$\mathbb{R}^d\f$ to \f$\mathbb{R}\f$.
 *
 *  The vector \f$\mathbf{x}\f$ is a \f$d\f$-dimensional Gaussian random vector:
 *  \f[
 *    \mathbf{x} \sim \mathcal{N}(\mu_x, \Sigma_x)
 *  \f]
 *  with covariance matrix \f$\Sigma_x = \mathrm{Diag}(\sigma_1^2, \ldots, \sigma_d^2)\f$.
 *
 *  The Gaussian noise \f$\tilde{\varepsilon}\f$ is centered with the following covariance
 *  matrix \f$\Sigma_\varepsilon = \mathrm{diag}(0,\ldots,0,\sigma^2,\ldots,\sigma^2)\f$.
 *
 **/
class GaussianAAModel : public IAAModel
                      , public IModel
{
  public:
    /** Constructor.
     *  @param data a reference on the data set to process
     **/
    GaussianAAModel( Matrix& data);

    /** virtual destuctor. */
    virtual ~GaussianAAModel();

    /** Set a new working data set.
     *  @param workData the working data set to use
     **/
    virtual void setWorkData(Matrix& workData);

    /** compute the log-likelihood of the model */
    void computeLogLikelihood();
};

} // namespace STK

#endif //STK_GAUSSIANAAMODEL_H
