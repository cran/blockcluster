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
 * Project:  stkpp::
 * created on: 27 oct. 2010
 * Purpose: Definition of the class MultidimRegression .
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MultidimRegression.cpp
 *  @brief In this file we implement the class MultidimRegression.
 **/

#include "../../Algebra/include/STK_LinAlgebra2D.h"
#include "../../Algebra/include/STK_LinAlgebra3D.h"
#include "../../Algebra/include/STK_GinvSymmetric.h"

#include "../../Algebra/include/STK_TExpAlgebra.h"

#include "../include/STK_MultidimRegression.h"


namespace STK
{

MultidimRegression::MultidimRegression( Matrix const* y, Matrix const* x)
                                      : IRegression<Matrix, Matrix, Vector>(y, x)
                                      , coefs_()
{ }

MultidimRegression::~MultidimRegression()
{ }

/* compute the regression function. */
void MultidimRegression::regression()
{
  // compute X'X
  MatrixSquare prod;
  multLeftTranspose(p_x_->asLeaf(), prod);

  // compute (X'X)^{-1}
  GinvSymmetric inv;
  inv(prod);

  // compute X'Y
  Matrix temp;
  multLeftTranspose(p_x_->asLeaf(), p_y_->asLeaf(), temp);

  // compute (X'X)^{-1}X'Y
  mult(prod, temp, coefs_);
}

/* compute the regression function. */
void MultidimRegression::regression(Vector const& weights)
{
  // compute X'WX
  MatrixSquare* prod = weightedMultLeftTranspose(p_x_->asLeaf(), weights);

  // compute (X'WX)^{-1}
  GinvSymmetric inv;
  inv(prod);

  // compute X'WY
  Matrix* temp = weightedMultLeftTranspose(p_x_->asLeaf(), p_y_->asLeaf(), weights);

  // compute (X'WX)^{-1}X'WY
  mult(*prod, *temp, coefs_);
  // remove temporary storages
  delete temp;
  delete prod;
}

/* Compute the predicted outputs by the regression function. */
void MultidimRegression::prediction()
{
  // remove existing predictions if any (should not be the case)
  if (!p_predicted_) p_predicted_ = new Matrix;
  // compute predictions
  mult(p_x_->asLeaf(), coefs_, *p_predicted_);
}


}
