/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2010  Serge Iovleff

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
 * Project:  stkpp::Regress
 * created on: 23 juin 2010
 * Purpose:  Interface base class for regression methods.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_IRegression.h
 *  @brief In this file we define the Interface base class IRegression.
 **/

#ifndef STK_IREGRESSION_H
#define STK_IREGRESSION_H

//#include "../../Sdk/include/STK_IRunner.h"

namespace STK
{

/** @ingroup Regress
 * @brief Interface base class for Regression methods.
 * The pure virtual function to implement are
 * @code
 *   virtual void regression();
 *   virtual void regression(WContainer const& weights);
 *   virtual void prediction();
 *   virtual Integer const& computeNbParameter() const;
 * @endcode
 * The virtual function
 * @code
 *   virtual void preRun();
 * @endcode
 * can be overloaded.
 */
template <class YContainer, class XContainer, class WContainer>
class IRegression
{
  protected:
    /** Constructor. Initialize the data members.
     * @param p_y container with the observed output
     * @param p_x container with the predictors (inputs) of the model
     */
    IRegression( YContainer const* p_y =0, XContainer const* p_x =0)
               : p_y_(p_y)
               , p_x_(p_x)
               , p_predicted_(0)
               , p_residuals_(0)
     {}

  public:
    /** virtual destructor. */
    virtual ~IRegression()
    { clear();}

    /** run the computations. */
    bool run()
    {
      // remove any existing storage
      preRun();
      // compute the regression
      regression();
      // Compute the number of parameter of the regression function.
      nbParameter_ = computeNbParameter();
      // predictions
      prediction();
      // compute residuals
      residuals();
      // return the result of the computations
      return true;
    }

    /** run the weighted computations.
     *  @param weights weights of the samples
     **/
    bool run( WContainer const& weights)
    {
      // perform any pre-operation needed befor the regression step
      preRun();
      // compute weighted regression
      regression(weights);
      // Compute the number of parameter of the regression function.
      nbParameter_ = computeNbParameter();
      // create container of the predicted value and compute prediction
      prediction();
      // create container of the residuals and compute them
      residuals();
      // return the result of the computations
      return true;
    }

    /** get the last error message.
     * @return the last error message
     **/
    inline String const& error() const { return msg_error_;}

    /** get the pointer of the container of the predicted values.
     * The container @c p_predicted_ will not be deleted by @c this.
     * @return the pointer on the predicted values
     **/
    inline YContainer* p_predicted() const
    { return p_predicted_;}
    /** get the pointer of the container of the residuals. The container
     *  @c p_residuals_ will not be deleted by @c this.
     *  @return the pointer on the residuals
     **/
    inline YContainer* p_residuals() const
    {  return p_residuals_;}

    /** Give the number of parameter of the regression function.
     *  @return the number of parameter of the regression function
     **/
    inline Integer const& nbParameter() const
    {  return nbParameter_;}

    /** Set the data set the regression method should use.
     * @param p_y data set to adjust
     * @param p_x data set of the predictors
     */
    void setData( YContainer const* p_y, XContainer const* p_x)
    { p_y_ = p_y; p_x_ = p_x;}

    /** Set the Y data set the regression method should use.
     * @param p_y data set to adjust
     */
    void setY( YContainer const* p_y)
    { p_y_ = p_y;}

    /** Set the X data set.
     * @param p_x data set of the predictors
     */
    void setX( XContainer const* p_x)
    { p_x_ = p_x;}

  protected:
    /** @brief perform any work needed before the call of the regression
     *  method.
     *  At this level do nothing.
     */
    virtual void preRun()
    { }

    /** @brief Compute the residuals of the model.
     * The residuals of the model are computed by computing the difference
     * between the observed outputs and the predicted outputs
     * of the model.
     */
    inline void residuals()
    {
      if (p_residuals_) delete p_residuals_;
       p_residuals_ =  p_y_->clone();
      *p_residuals_ = *p_y_ - *p_predicted_;
    }

    /** Container of the output to regress. */
    YContainer const* p_y_;
    /** Container of the regressors. */
    XContainer const* p_x_;
    /** Container of the predicted output. */
    YContainer* p_predicted_;
    /** Container of the residuals. */
    YContainer* p_residuals_;

  private:
    /** number of parameter of the regression method. */
    Integer nbParameter_;

    /** compute the regression function. */
    virtual void regression() =0;
    /** compute the weighted regression function.
     * @param weights the weights of the samples
     **/
    virtual void regression(WContainer const& weights) =0;
    /** Compute the predicted outputs by the regression function and store the
     * result in the p_predicted_ container. */
    virtual void prediction() =0;
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    virtual Integer computeNbParameter() const =0;
    /** delete allocated memory. */
    void clear()
    {
      if (p_predicted_) delete p_predicted_;
      if (p_residuals_) delete p_residuals_;
      p_predicted_ = 0;
      p_residuals_ = 0;
    }

  protected:
    /** String with the last error message. */
    String msg_error_;

};

}

#endif /* STK_IREGRESSION_H */
