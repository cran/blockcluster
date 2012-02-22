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
 * Project:  stkpp::Algebra
 * Purpose:  Define level 3 Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_LinAlgebra2D.h
 *  @brief Define level 3 Linear Algebra methods for two dimensional
 *  Containers containing Real.
 **/
 
#ifndef STK_LINALGEBRA3D_H
#define STK_LINALGEBRA3D_H

// Point and Vector classes
#include "../../Arrays/include/STK_Point.h"
#include "../../Arrays/include/STK_Vector.h"
#include "../../Arrays/include/STK_Matrix.h"
#include "../../Arrays/include/STK_MatrixSquare.h"

#include "STK_LinAlgebra1D.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Matrix multiplication
 *
 *  Perform the matricial product of the Matrix A with the Matrix B
 * \f$ C = AB \f$.
 *
 *  @param[in] A Matrix to multiply
 *  @param[in] B Matrix which multiply
 *  @param[out] C Matrix result
 **/
template < class TContainer2D1, class TContainer2D2, class TContainer2D3>
void mult( ITContainer2D< Real, TContainer2D1 > const& A
         , ITContainer2D< Real, TContainer2D2 > const& B
         , ITContainer2D< Real, TContainer2D3 >& C
         )
{
  // resize C
  C.resize(A.rangeVe(), B.rangeHo());
  // indexes
  const Integer firstCol = B.firstCol(), lastCol = B.lastCol();
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  // for all cols and for all rows
  for (Integer j=firstCol; j<=lastCol; j++)
    for (Integer i=firstRow; i<=lastRow; i++)
      C(i, j) = dot(A.row(i), B.col(j));
}

/** @ingroup Algebra
 *  @brief Matrix multiplication
 *
 *  Perform the matricial product of the Matrix A with the Matrix B
 * \f$ C = AB \f$.
 *
 *  @param[in] A Matrix to multiply
 *  @param[in] B Matrix which multiply
 *  @param[out] C Matrix result
 **/
template < class TContainer2D1, class TContainer2D2>
void mult( MatrixSquare const& A
         , ITContainer2D< Real, TContainer2D1 > const& B
         , ITContainer2D< Real, TContainer2D2 >& C
         )
{
  // resize C
  C.resize(A.range(), B.rangeHo());
  // indexes
  const Integer first = A.first(), last = A.last();
  const Integer firstCol = B.firstCol(), lastCol = B.lastCol();
  // for all cols and for all rows
  for (Integer j=firstCol; j<=lastCol; j++)
    for (Integer i=first; i<=last; i++)
      C(i, j) = dot(A.row(i), B.col(j));
}

/** @ingroup Algebra
 *  @brief Matrix multiplication
 *
 *  Perform the matrix product \f$ C=A'B \f$.
 *
 *  @param[in] A Left operand
 *  @param[in] B Right operand
 *  @param[out] C Matrix result
 **/
template < class TContainer2D1, class TContainer2D2, class TContainer2D3>
void multLeftTranspose( ITContainer2D< Real, TContainer2D2 > const& A
                      , ITContainer2D< Real, TContainer2D3 > const& B
                      , ITContainer2D< Real, TContainer2D1 >& C
                      )
{
    // create result
    C.resize(A.rangeHo(), B.rangeHo());
    // indexes
    // indexes
    const Integer firstColB = B.firstCol(), lastColB = B.lastCol();
    const Integer firstColA = A.firstCol(), lastColA = A.lastCol();
    // for all cols and for all cols
    for (Integer j=firstColB; j<=lastColB; j++)
      for (Integer i=firstColA; i<=lastColA; i++)
        C(i, j) = dot(A.col(i), B.col(j));
}


/** @ingroup Algebra
 *  @brief Matrix multiplication by its transpose
 *
 *  Perform the matrix product \f$ A'A \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A the matrix to multiply by itself
 *  @param[out] R the result
 **/
void multLeftTranspose( Matrix const& A, MatrixSquare& R);

/** @ingroup Algebra
 *  @brief Matrix multiplication by its transpose
 *
 *  Perform the matrix product \f$ AA' \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A the matrix to multiply by itself
 *  @param[out] R the result
 **/
void multRightTranspose( Matrix const& A, MatrixSquare& R);

/** @ingroup Algebra
 *  @brief MatrixSquare multiplication [DEPRECATED]
 *
 *  Perform the matrix product of the square matrix A with the square
 *  matrix B. The matrices A and B must have the the correct range.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] B Right operand
 *  @return a pointer on the result as a MatrixSquare
 **/
MatrixSquare* mult( MatrixSquare const& A, MatrixSquare const& B);

/** @ingroup Algebra
 *  @brief Weighted matrix multiplication [DEPRECATED]
 *
 *  Perform the weighted matrix product \f$ A'WB \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] B Right operand
 *  @param[in] weights the weights to apply
 *  @return a pointer on the result as a Matrix
 **/
Matrix* weightedMultLeftTranspose( Matrix const& A, Matrix const& B, Vector const& weights);

/** @ingroup Algebra
 *  @brief Weighted matrix multiplication [DEPRECATED]
 *
 *  Perform the matrix product \f$ A'WA \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] weights the weights to apply
 *  @return a pointer on the result as a MatrixSquare
 **/
MatrixSquare* weightedMultLeftTranspose( Matrix const& A, Vector const& weights);

} // Namespace STK

#endif /* STK_LINALGEBRA3D_H */
