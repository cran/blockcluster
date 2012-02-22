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
 * Project:  stkpp::regress
 * created on: 31 juil. 2010
 * Purpose:  implementation of the BSplineRegression class.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_BSplineRegression.cpp
 *  @brief In this file we implement the BSplineRegression class.
 **/


#include "../../Algebra/include/STK_LinAlgebra2D.h"
#include "../../Algebra/include/STK_LinAlgebra3D.h"
#include "../../Algebra/include/STK_GinvSymmetric.h"

#include "../../Algebra/include/STK_TExpAlgebra.h"

#include "../include/STK_BSplineRegression.h"

namespace STK
{

BSplineRegression::BSplineRegression( Matrix const* p_y
                                    , Vector const* p_x
                                    , Integer const& nbControlPoints
                                    , Integer const& degree
                                    , const _Kposition& position
                                    )
                                    : IRegression<Matrix, Vector, Vector>(p_y, p_x)
                                    , nbControlPoints_(nbControlPoints)
                                    , degree_(degree)
                                    , position_(position)
                                    , coefs_(p_x, nbControlPoints_, degree_, position_)
                                    , controlPoints_()
{ }

BSplineRegression::~BSplineRegression()
{}

/* compute the regression function. */
void BSplineRegression::regression()
{
  // compute X'X
  MatrixSquare prod;
  multLeftTranspose(coefs_.coefficients(), prod);

  // compute (X'X)^{-1}
  GinvSymmetric inv;
  inv(prod);

  // compute X'Y
  Matrix temp;
  multLeftTranspose(coefs_.coefficients(), p_y_->asLeaf(), temp);

  // compute (X'X)^{-1}X'Y
  mult(prod, temp, controlPoints_);
}

/* compute the regression function. */
void BSplineRegression::regression(Vector const& weights)
{
  // compute X'X
  MatrixSquare* prod = weightedMultLeftTranspose(coefs_.coefficients(), weights);

  // compute (X'X)^{-1}
  GinvSymmetric inv;
  inv(prod);

  // compute X'Y
  Matrix* temp = weightedMultLeftTranspose(coefs_.coefficients(), p_y_->asLeaf(), weights);

  // compute (X'X)^{-1}X'Y
  mult(*prod, *temp, controlPoints_);
  // remove temporary storages
  delete temp;
  delete prod;
}

/* Compute the predicted outputs by the regression function. */
void BSplineRegression::prediction()
{
  // create predictions if it does not exist any
  if (!p_predicted_) p_predicted_ = new Matrix;
  // compute predictions
  mult(coefs_.coefficients(), controlPoints_, *p_predicted_);
}

}
