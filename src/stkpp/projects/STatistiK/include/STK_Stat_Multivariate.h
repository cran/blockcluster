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
 * Project:  stkpp::StatistiK::StatDesc
 * Purpose:  Compute multivariate elementary statistics for a 2D container.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Stat_Multivariate.h
 *  @brief This file contain the declaration of the base class Multivariate.
 **/

#ifndef STK_STAT_MULTIVARIATE_H
#define STK_STAT_MULTIVARIATE_H

#include "../../Arrays/include/STK_Vector.h"
#include "../../Sdk/include/STK_IRunnerPtr2D.h"

namespace STK
{
namespace Stat
{
/** @ingroup StatDesc
 *  @brief Computation of the multivariate statistics of a Variable.
 *
 *  The template parameter @c TContainer2D is the type of container
 *  used for storing the data : It should derive from
 *  ITContainer2D and contain elements of type TYPE.
 *
 *  The template parameter TYPE is the Type of the data to analyze.
 **/
template < class TYPE, class TContainer2D >
class Multivariate : public IRunnerPtr2D< TYPE, TContainer2D>
{
  typedef IRunnerPtr2D< Real, TContainer2D> Runner2D;
  public:
    /** Constructor.
     *  Data set are initialized but no computation is done
     *  in this constructor. Statistics of the number of missing data
     *  and available data are delegated to derived classes.
     *  @param p_data the data set
     **/
    Multivariate( TContainer2D const* p_data)
                : Runner2D(p_data)
                , p_weights_(0)
                , nbSamples_(p_data->sizeVe())
                , nbVar_(p_data->sizeHo())
                , nMiss_(p_data->rangeHo(), 0)
//                , nbSamples_(p_data->rangeHo(), nbSamples_)
    { }

    /** virtual destructor.
     **/
    virtual ~Multivariate() { ;}
  
    /** Number of samples
     *  @return the number of samples in the data set (the number of rows)
     **/
    inline Integer const& nbSamples() const {return nbSamples_;}

    /** Number of variables
     * @return the number of variables in the data set (the number of columns)
     **/
    inline Integer const& nbVar() const {return nbVar_;}

    /** Number of missing values
     * @return An array with the number of missing values for each variables
     * of the data set
     **/
    inline Array1D<Integer> const& nbMissingSamples() const {return nMiss_;}

    /** Number of observed values
     * @return An array with the number of observed values for each variables
     * of the data set (not missing).
     **/
    inline Array1D<Integer> const& nbAvailableSamples() const {return nbAvailable_;}

    /** run the estimation of the Multivariate statistics. **/
    virtual bool run()
    {
      nbSamples_ = this->p_data_->sizeVe();
      nbVar_ = this->p_data_->sizeHo();
      nMiss_.resize(this->p_data_->rangeHo());
      nMiss_ = 0;
      nbAvailable_.resize(this->p_data_->rangeHo());
      // get dimensions
      const Integer first_ind = this->p_data_->firstRow();
      const Integer last_ind  = this->p_data_->lastRow();
      const Integer first_var = this->p_data_->firstCol();
      const Integer last_var  = this->p_data_->lastCol();
      // for each variables
      for (Integer j= first_var; j<= last_var; j++)
      {
        // number of not missing observations
        Integer nobs = nbSamples_;
        // compute the mean
        for (Integer i= first_ind; i<= last_ind; i++)
          if (!Arithmetic<TYPE>::isFinite((*this->p_data_)(i,j))) nobs--;
        nbAvailable_[j] = nobs;
        nMiss_[j] = nbSamples_ - nobs;
      }
      return true;
    }

    /** run the estimation of the weighted multivariate statistics.
     * @param weights the weights of the samples
     **/
    virtual void run( Vector const& weights)
    {
      p_weights_ = &weights;
      run();
    }

  protected:
    /** Vector of the weights */
    Vector const* p_weights_;
    /** number of samples */
    Integer nbSamples_;
    /** Number of variables */
    Integer nbVar_;
    /** number of missing data of each variables */
    Array1D<Integer> nMiss_;
    /** number of observed data of each variables */
    Array1D<Integer> nbAvailable_;
};

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_MULTIVARIATE_H*/
