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
 * Purpose: define the Interface base class ITModel (Statistical Model).
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ITModel.h
 *  @brief In this file we define the Interface base class ITModel.
 **/

#ifndef STK_ITMODEL_H
#define STK_ITMODEL_H

#include <cmath>

#include "STK_IModel.h"
#include "../../Sdk/include/STK_IRunnerPtr2D.h"
#include "../../STatistiK/include/STK_Law_ITMultivariate.h"

namespace STK
{

/** @ingroup Model
 *  @brief Interface base class for all statistical models.
 *  A statistical model is defined with various elements
 *  - A data set with the samples in rows and the variables in column
 *  - A parameterized probability density
 *
 *  We are making the assumption that the variables are all of the same
 *  type and stored in a class deriving from ITContainer2D.
 **/
template < class TYPE, class TContainer2D>
class ITModel : public IModel
              , public IRunnerPtr2D< TYPE, TContainer2D>
{
  /** Type of the container containing the data*/
  typedef ITContainer2D<TYPE, TContainer2D> Container2D;
  /** Type of the Runner*/
  typedef IRunnerPtr2D< TYPE, TContainer2D> Runner2D;
  /** Type of the column container */
  typedef typename TContainer2D::TContainerHo TContainerHo;
  /** Type of the law */
  typedef Law::ITMultivariate<TYPE, TContainerHo> MultiLaw;

  protected:
    /** Constructor */
    ITModel( Container2D const* p_data)
           : IModel()
           , Runner2D(p_data)
           , p_law_(0)
    {
      if (this->p_data_)
      {
        nbSample_ = this->p_data_->sizeVe();
        nbVar_ = this->p_data_->sizeHo();
      }
    }

  public:
    /** virtual destructor */
    virtual ~ITModel() { if (p_law_) delete p_law_;}

    /** give the probability law of the model
     * @return a constant pointer on the probability law.
     **/
    inline MultiLaw const* p_law() const { return p_law_;}

    /** Set the probability law of the model.
     *  @param p_law the probability law of the model
     */
    void setLaw( MultiLaw* p_law)
    {
      //
      if (p_law_) delete p_law_;
      p_law_ = p_law;
    }

  protected:
    /** a pointer on the probability law. */
    MultiLaw* p_law_;

    /** compute the log Likelihood of the statistical model. */
    virtual void compLogLikelihood()
    {
      // no data
      if (!this->p_data_) return;
      // check there exists a law
      if (!p_law_) throw runtime_error("In ITModel::compLogLikelihood() "
                                            "p_law_ is not initialized.");
      // get dimensions of the samples and sum over all log-likelihood values
      const Integer first = this->p_data_->firstRow(), last = this->p_data_->lastRow();
      Real sum = 0.0;
      for (Integer i=first; i<= last; i++)
      {
        sum += p_law_->lpdf(this->p_data_->row(i));
      }
      logLikelihood_ = sum;
    }
};

} // namespace STK

#endif /* STK_ITModel */
