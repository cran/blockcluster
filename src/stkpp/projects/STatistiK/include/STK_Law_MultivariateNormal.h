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
 * Project:  stkpp::STatistiK::Law
 * created on: 29 juil. 2011
 * Purpose:  define the templated MultivariateNormal law.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Law_MultivariateNormal.h
 *  @brief In this file we define the templated MultivariateNormal law.
 **/

#ifndef STK_MUTIVARIATENORMAL_H
#define STK_MUTIVARIATENORMAL_H

#include "STK_Law_ITMultivariate.h"
#include "STK_Law_Normal.h"

#include "../../Algebra/include/STK_EigenvaluesSymmetric.h"

#include "../include/STK_Const_Math.h"

#include "../../Algebra/include/STK_LinAlgebra1D.h"
#include "../../Algebra/include/STK_LinAlgebra2D.h"

namespace STK
{

namespace Law
{

/** @ingroup Laws
 *  @brief  Class for the multivariate Normal distribution.
 *
 *  In probability theory and statistics, the "multivariate normal distribution"
 *  or "multivariate Gaussian distribution", is a generalization of the
 *  one-dimensional (univariate) @link Normal normal distribution @endlink to
 *  higher dimensions.  A random vector is said to be multivariate normally
 *  distributed if every linear combination of its components has a univariate
 *  normal distribution.  The multivariate normal distribution is often used to
 *  describe, at least approximately, any set of (possibly) correlated
 *  real-valued random variables each of which clusters around a mean value.
 *
 *  The multivariate normal distribution of a \f$ p\f$-dimensional random vector
 *  \f[
 *  \mathbf{X} = \left(
 *   X_1, X_2, \ldots, X_p
 *               \right)'
 *  \f]
 *  can be written in the following notation
 *  \f[
 *  \mathbf{X}\ \sim\ \mathcal{N}(\mu,\ \Sigma).
 *  \f]
 *  with \f$ p \f$-dimensional mean vector
 *  \f[
 *  \mu = \left(
 *        \mathrm{E}[X_1], \mathrm{E}[X_2], \ldots, \mathrm{E}[X_k]
 *        \right)'
 *  \f]
 *  and \f$ p \times p \f$ covariance matrix
 *  \f[
 *  \Sigma = [\mathrm{Cov}(X_i, X_j)]_{i=1,2,\ldots,p;\ j=1,2,\ldots,p}
 *  \f]
 */
template <class Container1D>
class MultivariateNormal: public ITMultivariate<Real, Container1D>
{
  public:
    /** Constructor.
     *  @param mu mean of the Normal distribution
     *  @param sigma covariance matrix of the Normal distribution
     **/
    MultivariateNormal( Container1D const& mu, MatrixSquare const& sigma)
                      : ITMultivariate<Real, Container1D>(_T("Normal Multivariate"))
                      , mu_(mu)
                      , sigma_(sigma)
                      , decomp_(&sigma_)
                      , invEigenvalues_(mu_.range())
                      , squareroot_(sigma_.range())
    { update();}

    /** destructor. */
    virtual ~MultivariateNormal() { ;}

    /** give the location parameter of the gaussian distribution.
     *  @@return the location parameter
     **/
    inline Container1D const& mu() { return mu_;}
    /** give the location parameter of the gaussian distribution.
     *  @@return the location parameter
     **/
    inline MatrixSquare const& sigma() { return sigma_;}
    /** give the eigenvalue decomposition of the covariance matrix sigma.
     *  @@return the eigenvalue decomposition
     **/
    inline EigenvaluesSymmetric const& decomp() { return decomp_;}
    /** update the parameters specific to the law. */
    virtual void update()
    {
      // check dimensions
      if (mu_.range() != sigma_.range())
      { throw runtime_error(_T("In MultivariateNormal::MultivariateNormal(mu, sigma) "
                               "mu_.range() != sigma_.range().\n"));
      }
      // decomposition of the covariance matrix
      if (!decomp_.run())
      { throw runtime_error(_T("In MultivariateNormal::MultivariateNormal(mu, sigma) "
                               "decomposition of sigma fail.\n"));
      }
      // compute the inverse of the eigenvalues of sigma_ and the squareroot_
      // matrix needed by rand
      invEigenvalues_.resize(mu_.range());
      squareroot_.resize(mu_.range());
      // get dimension
      Integer first = mu_.first(), rank = decomp_.rank(), last = first + rank -1;
      for (Integer j=first; j<= last; j++)
      {
        invEigenvalues_[j] = 1./decomp_.eigenvalues()[j];
        Real squareRoot = sqrt((double)decomp_.eigenvalues()[j]);
        for (Integer i= first; i <= last; ++i)
        {
          squareroot_(i,j) = decomp_.rotation()(i,j)*squareRoot;
        }
      }
      first = last+1; last = mu_.last();
      for (Integer j=first; j<= last; j++)
      {
        invEigenvalues_[j] = 0.;
        Real squareRoot = sqrt((double)decomp_.eigenvalues()[j]);
        for (Integer i= first; i <= last; ++i)
        {
          squareroot_(i,j) = decomp_.rotation()(i,j)*squareRoot;
        }
      }
    }
    /** @brief compute the probability distribution function (density) of the
     * multivariate normal law
     * \f[
     *  f(x) = \frac{1}{ (2\pi)^{k/2}|\Sigma|^{1/2} }
     *           \exp\!\left( {-\tfrac{1}{2}}(x-\mu)'\Sigma^{-1}(x-\mu) \right),
     * \f]
     *  Give the value of the pdf at the point x.
     *  @param x the multivariate value to compute the pdf.
     *  @return the value of the pdf
     **/
    virtual Real pdf( ITContainer1D<Real, Container1D> const& x) const
    {
      // check determinant is not 0
      if (decomp_.det() == 0.)
      { throw runtime_error("In MultivariateNormal::pdf(x) "
                                 "|sigma| == 0.\n"
                                 );
      }
      // check ranges
      if (x.range() != mu_.range() )
      {
        throw runtime_error("In MultivariateNormal::pdf(x)"
                                 "x.range() != mu_.range()\n"
                                );
      }
      return exp((double)lpdf(x));
    }

    /** @brief compute the log probability distribution function.
     *  Give the value of the log-pdf at the point x.
     *  @param x the multivariate value to compute the lpdf.
     *  @return the value of the log-pdf
     **/
    Real lpdf( ITContainer1D<Real, Container1D> const& x) const
    {
      // check ranges
      if (x.range() != mu_.range() )
      {
        throw runtime_error(_T("In MultivariateNormal::pdf(x)"
                               "x.range() != mu_.range()\n"));
      }
      // compute x - mu
      Vector xbar;
      diff(x, mu_, xbar);
      // compute P'(x-mu)
      Vector xrot;
      multLeftTranspose(decomp_.rotation(), xbar, xrot);
      // compute pdf using (x-mu)'PD^{-1}P'(x-mu)
      Real res = 0.5 * weightedNormTwo2<Vector, Vector>(xrot, invEigenvalues_)
               + invEigenvalues_.size() * Const::_LNSQRT2PI_
               + 0.5 * log((double)decomp_.det());
      return -res;
    }

    /** @brief compute the log likehood of a data set.
     *  sum the values of the log-pdf at the points stored in x.
     *  @param data the multivariate values to compute the lpdf.
     *  @return the value of the log-pdf
     **/
    template < class TContainer2D>
    Real logLikelihood( TContainer2D const& data) const
    {
      // check ranges
      if (data.rangeHo() != mu_.range() )
      {
        throw runtime_error(_T("Error in n MultivariateNormal::logLikehood(x)\nWhat: "
                               "data.rangeHo() != mu_.range()\n"));
      }

      // compute x - mu
      Vector xres(mu_.range()), xrot(mu_.range());
      // get dimensions of the samples and sum over all log-likelihood values
      const Integer first = data.firstRow(), last = data.lastRow();
      Real sum = 0.0;
      for (Integer i=first; i<= last; i++)
      {
        // compute residual
        xres = data(i) - mu_;
        // compute P'(x-mu)
        multLeftTranspose(decomp_.rotation(), xres, xrot);
        // compute lpdf using (x-mu)'PD^{-1}P'(x-mu)
        Real res = 0.5 * weightedNormTwo2(xrot, invEigenvalues_);
        sum += res;
      }
      sum += data.sizeVe()*( invEigenvalues_.size() * Const::_LNSQRT2PI_
                           + 0.5 * log((double)decomp_.det())
                           );
      return -sum;
    }

    /** @brief simulate a realization of the Multivariate Law and store the
     *  result in x.
     *  @param[out] x the simulated value.
     **/
    virtual void rand( ITContainer1D<Real, Container1D>& x) const
    {
      // create intermediary container
      Container1D x_iid(x.range());
      // fill it with iid N(0,1) variates
      Normal(0., 1.).rand1D<Container1D>(x_iid);
      // rotate with squareroot_
      mult(squareroot_, x_iid, x);
      // translate with mu_
      x.asLeaf() += mu_;
    }
    /** @brief simulate a realization of the Multivariate Law and store the
     *  result in x (using a reference vector). The class Container1D have
     *  to derive from IContainerRef.
     *  @param[out] x the simulated value.
     **/
    virtual void rand( ITContainer1D<Real, Container1D> const& x) const
    {
      // create intermediary container
      Container1D x_iid(x.range());
      // create intermediary reference container
      Container1D x_ref(x.asLeaf(), true);
      // fill it with iid N(0,1) variates
      Normal(0., 1.).rand1D<Container1D>(x_iid);
      // rotate with squareroot_
      mult(squareroot_, x_iid, x_ref);
      // translate with mu_
      x_ref += mu_;
    }

  protected:
    /** The position parameter. **/
    Container1D mu_;
    /** The covariance parameter. **/
    MatrixSquare sigma_;
    /** the decomposition in eigenvalues of the covariance matrix*/
    EigenvaluesSymmetric decomp_;
    /** inverse of the eigenvalues of sigma_ */
    Vector invEigenvalues_;
    /** The square root of the matrix @c Sigma_. **/
    MatrixSquare squareroot_;
};

} // namespace Law

} // namespace STK

#endif /* STK_MUTIVARIATENORMAL_H */
