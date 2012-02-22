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
 * created on: 31 juil. 2010
 * Purpose: definition of the BsplineRegression class.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_BSplineRegression.h
 *  @brief In this file we define the BsplineRegression class.
 **/

#ifndef STK_BSPLINEREGRESSION_H
#define STK_BSPLINEREGRESSION_H

#include "../../Arrays/include/STK_Matrix.h"

#include "STK_BSplineCoefficients.h"
#include "STK_IRegression.h"

namespace STK
{

/** @brief Compute a BSpline, multivalued, regression function using BSpline
 *  basis.
 *
 */
class BSplineRegression : public IRegression<Matrix, Vector, Vector>
{
  private:
    typedef BSplineCoefficients::KnotsPosition _Kposition;

  public:
    /** Constructor.
     * @param p_y d-dimensional array of output to fit
     * @param p_x uni-dimensional array of predictor
     * @param nbControlPoints number of control points of the spline
     * @param degree degree of the BSpline basis
     * @param position position of the knots to used
     **/
    BSplineRegression( Matrix const* p_y
                     , Vector const* p_x
                     , Integer const& nbControlPoints
                     , Integer const& degree = 3
                     , const _Kposition& position = BSplineCoefficients::uniform_
                     );

    /** virtual destructor. */
    virtual ~BSplineRegression();

    /** give the degree of the B-Spline curve.
     *  @return the degree of the B-Spline curve
     * */
    inline Integer const& degree() const
    { return degree_;}
    /** give the number of control points of the B-Spline curves.
     *  @return the number of control points of the B-Spline curve
     **/
    inline Integer const& nbControlPoints() const
    { return nbControlPoints_;}
    /** give the control points.
     *  @return the control points of the B-Spline curve
     **/
    inline Matrix const& controlPoints() const
    { return controlPoints_; }
    /** give the knots.
     *  @return the knots of the B-Spline curve
     **/
    inline Vector const& knots() const
    { return coefs_.knots(); }
    /** give the computed coefficients of the B-Spline curves.
     *  This is a matrix of size (p_x_->range(), 0:lastControlPoints).
     *  @return the coefficients of the B-Spline curve
     **/
    inline Matrix const& coefficients() const
    { return coefs_.coefficients();}

  protected:
    /** number of control points of the B-Spline curve. */
    Integer nbControlPoints_;
    /** degree of the B_Spline curve */
    Integer degree_;
    /** method of position of the knots of the B-Spline curve */
    _Kposition position_;
    /** Coefficients of the regression matrix */
    BSplineCoefficients coefs_;
    /** Estimated control points of the B-Spline curve */
    Matrix controlPoints_;

    /** compute the coefficients of the BSpline basis. This method will be
     *  called in the base class @c IRegression::run()
     **/
    inline virtual void preRun() {coefs_.run();}

    /** compute the regression function. This method will be
     *  called in the base class @c IRegression::run() after preRun()
     **/
    virtual void regression();
    /** compute the regression function. This method will be
     *  called in the base class @c IRegression::run(weights) after preRun()
     *  @param weights the weights of the samples
     **/
    virtual void regression(Vector const& weights);
    /** Compute the predicted outputs by the regression function. This method
     *  will be called in the base class @c IRegression::run() after preRun()
     **/
    virtual void prediction();
    /** Compute the number of parameter of the regression function.
     * @return the number of parameter of the regression function
     **/
    inline virtual Integer computeNbParameter() const
    { return controlPoints_.sizeHo() * controlPoints_.sizeVe(); }

};

} // namespace STK

#endif /* STK_BSPLINEREGRESSION_H */
