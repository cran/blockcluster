/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  StatDesc
 * Purpose:  Compute multivariate elementary statistics for
 * a 2D container.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Stat_MultivariateReal.h
 *  @brief In this file we specialize the class Multivariate to real type.
 **/

#ifndef STK_STAT_MULTIVARIATEREAL_H
#define STK_STAT_MULTIVARIATEREAL_H

#include "../../STKernel/include/STK_Arithmetic.h"
#include "../../STKernel/include/STK_Misc.h"

#include "../../Arrays/include/STK_MatrixSquare.h"

#include "STK_Stat_UnivariateReal.h"
#include "STK_Stat_BivariateRealReal.h"
#include "STK_Stat_Multivariate.h"

namespace STK
{
namespace Stat
{

typedef Multivariate<Real, Matrix> MultivariateMatrix;

/** @ingroup StatDesc
 *  @brief Computation of the Multivariate Statistics of a 2D Container
 *  of Real.
 *
 *  This partial specialization if the class @c Multivariate is just a factory
 *  class for computing the mean, the variance and the covariance Matrix of a
 *  p_data set stored in a @c TContainer2D with n rows (the samples) and
 *  p columns (the variables).
 *
 *  The p_data set will be copied as a reference of the original p_data set.
 *
 *  The p_data set can be weighted with a vector of weight.
 **/
template < class TContainer2D >
class Multivariate<Real, TContainer2D > : public IRunnerPtr2D< Real, TContainer2D>
{
  typedef typename TContainer2D::TContainerHo TContainerHo;
  typedef typename TContainer2D::TContainerVe TContainerVe;
  /** type of runner */
  typedef IRunnerPtr2D< Real, TContainer2D> Runner2D;

  public:
    /** Constructor.
     *  Compute the Multivariate statistics of the Matrix p_data set.
     *  @param p_data the p_data set
     **/
    Multivariate( TContainer2D const* p_data)
                : Runner2D(p_data)
                , nbSamples_(0)
                , nbVar_(0)
                , mean_(0)
                , var_(0)
                , cov_(0)
    { update(); }

    /** virtual destructor.*/
    virtual ~Multivariate() { }

    /** get the number of variables
     * @return the number of variables in the p_data set (the number of columns)
     **/
    inline Integer const& nbVar() const {return nbVar_;}
    /** get the umber of samples
     *  @return the number of samples in the p_data set (the number of rows)
     **/
    inline Integer const& nbSamples() const {return nbSamples_;}
    /** get the Vector of the mean
     * @return the mean of the variables */
    inline TContainerHo const& mean() const { return mean_;}
    /** get the vector of the variance of the Variables
     * @return the variance of the variables
     **/
    inline TContainerHo const& variance() const { return var_;}
    /** Matrix of the covariance of the variables
     * @return the covariance of the variables
     * */
    inline MatrixSquare const& covariance() const { return cov_;}
    /** run the estimation of the Multivariate statistics. **/
    virtual bool run()
    {
      try
      {
        if (!this->p_data_)
          throw runtime_error("data have not be set.");
        // get dimensions
        const Integer first = this->p_data_->firstCol(), last  = this->p_data_->lastCol();
        // for each variables
        for (Integer j= first; j<= last; j++)
        {
          mean_[j] = Stat::mean(this->p_data_->col(j));
          // compute the variance
          var_[j]  = varianceWithFixedMean(this->p_data_->col(j), mean_[j]);
          cov_(j, j) = var_[j];
          // compute the covariances
          for (Integer i= first; i<j; i++)
          {
            cov_(i, j) = covarianceWithFixedMean(this->p_data_->col(i), this->p_data_->col(j), mean_[i], mean_[j]);
            cov_(j, i) = cov_(i, j);
          }
        }
      }
      catch (Exception error)
      {
        this->msg_error_  = _T("Error in MultivariateReal::run():\nWhat: ");
        this->msg_error_ += error.error();
        return false;
      }
      // no error
      return true;
    }

    /** run the estimation of the weighted multivariate statistics.
     * @param weights the weights of the samples
     **/
    virtual bool run( TContainerVe const& weights)
    {
      try
      {
        if (!this->p_data_)
          throw runtime_error("data have not be set.");
        // get dimensions
        const Integer first = this->p_data_->firstCol(), last  = this->p_data_->lastCol();
        // for each variables
        for (Integer j= first; j<= last; j++)
        {
          mean_[j] = Stat::mean(this->p_data_->col(j), weights);
          // compute the variance
          var_[j]  = varianceWithFixedMean(this->p_data_->col(j), weights, mean_[j]);
          cov_(j, j) = var_[j];
          // compute the covariances
          for (Integer i= first; i<j; i++)
          {
            cov_(i, j) = covarianceWithFixedMean(this->p_data_->col(i), this->p_data_->col(j), weights, mean_[i], mean_[j]);
            cov_(j, i) = cov_(i, j);
          }
        }
      }
      catch (Exception error)
      {
        this->msg_error_  = _T("Error in MultivariateReal::run(weights): ");
        this->msg_error_ += error.error();
        return false;
      }
      // no error
      return true;
    }

  protected:
    /** number of samples */
    Integer nbSamples_;
    /** Number of variables */
    Integer nbVar_;

    /** Vector of the mean of the Variables */
    TContainerHo mean_;
    /** Vector of the variance of the variables */
    TContainerHo var_;
    /** Matrix of the covariance of the variables */
    MatrixSquare cov_;

    /** udpating method in case we set a new data set */
    virtual void update()
    {
      // if there is no data there is nothing to update
      if (this->p_data_)
      {
        nbSamples_ = this->p_data_->sizeVe();
        nbVar_ = this->p_data_->sizeHo();
        mean_.resize(this->p_data_->rangeHo());
        var_.resize(this->p_data_->rangeHo());
        cov_.resize(this->p_data_->rangeHo());
      }
    }
};


/** @ingroup StatDesc
 *  Compute the mean of the container V and store the result in mean
 *  @param V the container with the Data
 *  @param mean the container with the result
 **/
template < class TContainerHo, class TContainer2D >
void mean( ITContainer2D< Real, TContainer2D> const&  V
         , ITContainer1D< Real, TContainerHo> & mean
         )
{
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  mean.resize(V.rangeHo());
  for (Integer j= firstVar; j<= lastVar; j++)
    mean[j] = Stat::mean(V.col(j));
}

/**  @ingroup StatDesc
 *  Compute the (weighted) mean of the variable V
 *  \f[ \hat{\mu} = \frac{1}{\sum_{i=1}^n W(i)}
 *                  \sum_{i=1}^n W(i) V(i)
 *  \f]
 * if there is weights.
 *  @param V the variable
 *  @param W the weights
 *  @param mean the computed weighted mean
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
void mean( ITContainer2D< Real, TContainer2D> const& V
         , ITContainer1D< Real, TContainerVe> const& W
         , ITContainer1D< Real, TContainerHo> & mean
         )
{
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  mean.resize(V.rangeHo());
  for (Integer j= firstVar; j<= lastVar; j++)
    mean[j] = Stat::mean(V.col(j), W);
}

/** @ingroup StatDesc
 *  Compute the mean of the container V.
 *  @param V the container with the Data
 **/
template < class TContainerHo, class TContainer2D >
TContainerHo mean( ITContainer2D< Real, TContainer2D> const&  V)
{
  TContainerHo mean(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
    mean[j] = Stat::mean(V.col(j));
  return mean;
}

/**  @ingroup StatDesc
 *  Compute the (weighted) mean of the variable V
 *  \f[ \hat{\mu} = \frac{1}{\sum_{i=1}^n W(i)}
 *                  \sum_{i=1}^n W(i) V(i)
 *  \f]
 * if there is weights.
 *  @param V the variable
 *  @param W the weights
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
TContainerHo mean( ITContainer2D< Real, TContainer2D> const& V
                 , ITContainer1D< Real, TContainerVe> const& W
                 )
{
  // create result
  TContainerHo mean(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
    mean[j] = Stat::mean(V.col(j), W);
  return mean;
}

/**  @ingroup StatDesc
 *  Compute the variance of the variable V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n-1}
 *                       \sum_{i=1}^n (V(i)-\hat{\mu})^2.
 *  \f]
 *  @param V variable
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainer2D >
typename TContainer2D::TContainerHo variance( ITContainer2D< Real, TContainer2D> const& V
                     , bool unbiased = false
                     )
{
  // create result
  typename TContainer2D::TContainerHo var(V.rangeHo());
  // get dimensions
  const Integer firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    var[j] = Stat::variance(V.col(j), unbiased);
  }
  // return variance
  return var;
}

/** @ingroup StatDesc
 *  Compute the weighted variance of the variable V.
 *  \f[ \hat{\sigma}^2 = \frac{\sum_{i=1}^n W(i)}
 *                          {\left( \sum_{i=1}^n W(i))\right)^2
 *                                - \sum_{i=1}^n W(i)^2}
 *                     \sum_{i=1}^n W(i) (V(i)-\hat{\mu})^2.
 *  \f]
 * If there is no weights, this definition reduces to the usual
 * definition of the variance with factor 1/(n-1).
 *  @param V variable
 *  @param W weights
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
typename TContainer2D::TContainerHo variance( ITContainer2D< Real, TContainer2D> const& V
                     , TContainerVe const& W
                     , bool unbiased = false
                     )
{
  // create result
    typename TContainer2D::TContainerHo var(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  // check if there exist a mean vector
  for (Integer j= firstVar; j<= lastVar; j++)
    var[j] = Stat::variance(V.col(j), W, unbiased);
  // return variance
  return var;
}

/**  @ingroup StatDesc
 *  Compute the variance of the variable V with fixed mean.
 *  \f[ \hat{\mu} = \frac{1}{n}
 *                  \sum_{i=1}^n (V(i) - \mu)^2.
 *  \f]
 *  @param V variable
 *  @param mu mean of the variable V
 * @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
TContainerHo varianceWithFixedMean( ITContainer2D< Real, TContainer2D> const& V
                                  , TContainerHo const& mu
                                  , bool unbiased = false
                                  )
{
  if (mu.range() != V.rangeHo())
    throw runtime_error("In varianceWithFixedMean(V, mu, unbiased) "
                             "mu.range() != V.rangeHo()");
  // create result
  TContainerHo var(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
    var[j] = Stat::varianceWithFixedMean(V.col(j), mu[j], unbiased);
  // return variance
  return var;
}

/**  @ingroup StatDesc
 *  Compute the weighted variance of the variables V with fixed mean.
 *  \f[ \hat{\mu} = \frac{1}{\sum_{i=1}^n W(i)}
 *                  \sum_{i=1}^n W(i) (V(i) - \mu)^2
 *  \f]
 *  @param V variable
 *  @param W weights
 *  @param mu weighted mean of V
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
TContainerHo varianceWithFixedMean( ITContainer2D< Real, TContainer2D> const& V
                                  , TContainerVe const& W
                                  , TContainerHo const& mu
                                  , bool unbiased = false
                                  )
{
#ifdef STK_DEBUG
  if (mu.range() != V.rangeHo())
    throw runtime_error(_T("In varianceWithFixedMean(V, W, mu) "
                           "mu.range() != V.rangeHo()"));
#endif
  // create result
  TContainerHo var(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
    var[j] = Stat::varianceWithFixedMean<TContainerVe>(V[j], W, mu[j], unbiased);
  // return variance
  return var;
}

/**  @ingroup StatDesc
 *  Compute the covariance of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n (V_i-\hat{\mu}) (V_i-\hat{\mu})'.
 *  \f]
 *  @param V variable
 *  @param cov the computed covariance
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainer2D >
void covariance( TContainer2D const& V
               , MatrixSquare & cov
               , bool unbiased = false
               )
{
  typename ContainerTraits<Real, TContainer2D>::TContainerHo mean;
  // compute the
  Stat::mean(V, mean);
  // get dimensions
  const Integer firstVar = V.firstCol(), lastVar = V.lastCol();
  cov.resize(V.rangeHo());
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    cov(j, j) = Stat::varianceWithFixedMean(V.col(j), mean[j], unbiased);
    for (Integer i= firstVar; i<j; i++)
    {
      cov(j,i) = ( cov(i, j) = Stat::covarianceWithFixedMean(V.col(i), V.col(j), mean[i], mean[j], unbiased));
    }
  }
}

/**  @ingroup StatDesc
 *  Compute the weighted covariance of the data set V.
 *  \f[ \hat{\sigma}^2 = \frac{1}{n}
 *                       \sum_{i=1}^n  w_{i} w_j (V_i-\hat{\mu}) (V_i-\hat{\mu})'.
 *  \f]
 *  @param V the variable
 *  @param W the weights
 *  @param cov the computed covariance
 *  @param unbiased @c true if we want an unbiased estimator of the variance,
 *  @c false otherwise (default is @c false)
 **/
template < class TContainerVe, class TContainer2D >
void covariance( ITContainer2D< Real, TContainer2D> const& V
               , ITContainer1D< Real, TContainerVe> const& W
               , MatrixSquare & cov
               , bool unbiased = false
               )
{
  typename ContainerTraits<Real, TContainer2D>::TContainerHo mean;
  // compute the
  Stat::mean(V, W, mean);
  // get dimensions
  const Integer firstVar = V.firstCol(), lastVar = V.lastCol();
  cov.resize(V.rangeHo(), V.rangeHo());
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    cov(j, j) = Stat::varianceWithFixedMean(V.col(j), W, unbiased);
    for (Integer i= firstVar; i<j; i++)
    {
      cov(j,i) = ( cov(i, j) = Stat::covarianceWithFixedMean(V.col(i), V.col(j), W, mean[i], mean[j], unbiased));
    }
  }
}


}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_MULTIVARIATEREAL_H*/
