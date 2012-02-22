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
 * Project:  Algebra
 * Purpose:  Implement Level 3 Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_LinAlgebra2D.cpp
 *  @brief Implement Linear Algebra methods with Real for two dimensional
 *  containers.
 **/
 
#include "../include/STK_LinAlgebra3D.h"

namespace STK
{
/* @ingroup Algebra
 *  @brief Matrix multiplication by its transpose
 *
 *  Perform the matrix product \f$ A'A \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A the matrix to multiply by itself
 *  @param[out] R the result
 **/
void multLeftTranspose( Matrix const& A, MatrixSquare &R)
{
  R.resize(A.rangeHo());
  // indexes
  const Integer first = A.firstCol(), last = A.lastCol();
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  //
  for (Integer i=first; i<=last; i++)
  {
    // diagonal
    Real sum = 0.0;
    for (Integer k = firstRow; k<=lastRow; k++)
      sum += A(k, i) * A(k, i);
    R(i, i) = sum;
    // outside diagonal
    for (Integer j=first; j<i; j++)
    {
      Real sum = 0.0;
      for (Integer k = firstRow; k<=lastRow; k++)
        sum += A(k, i) * A(k, j);
      R(j, i) = (R(i, j) = sum);
    }
  }
}

/* @ingroup Algebra
 *  @brief Matrix multiplication [DEPRECATED]
 *
 *  Perform the matrix product \f$ AA' \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A the matrix to multiply by itself
 *  @param[out] R the result
 **/
void multRightTranspose( Matrix const& A, MatrixSquare& R)
{
  // create result
  R.resize(A.rangeVe());
  // get dimensions
  const Integer first = A.firstRow(), last = A.lastRow();
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
  // compute AA'
  for (Integer i=first; i<=last; i++)
  {
    // compute diagonal
    Real sum = 0.0;
    for (Integer k = firstCol; k<=lastCol; k++)
      sum += A(i, k) * A(i, k);
      R(i, i) = sum;
    // compute outside diagonal
    for (Integer j=first; j<i; j++)
    {
      Real sum = 0.0;
      for (Integer k = firstCol; k<=lastCol; k++)
        sum += A(i, k) * A(j, k);
      R(j, i) = (R(i, j) = sum);
    }
  }
}


/* @ingroup Algebra
 *  @brief Matrix multiplication
 *
 *  Perform the matrix product of the Matrix A with the Matrix B.
 *
 * The result is created on the stack. This
 * method should be specialized for upper and lower triangular
 * Matrix.
 *  @param[in] A Matrix to multiply
 *  @param[in] B Matrix which multiply
 *  @return a pointer on the result as a LEAF
 **/
MatrixSquare* mult( MatrixSquare const& A, MatrixSquare const& B)
{
#ifdef STK_DEBUG
  if (A.range() != B.range())
    throw runtime_error("In MatrixSquare* mult(A, B) "
                             "A.range() != B.range()");
#endif
  // create result
  MatrixSquare* res = new MatrixSquare(A.range());
  // indexes
  const Integer first = A.firstCol(), last = A.lastCol();
  for (Integer i=first; i<=last; i++)
  {
    for (Integer j=first; j<=last; j++)
    {
      Real sum = 0.0;
      for (Integer k = first; k<=last; k++)
        sum += A(i, k) * B(k, j);
      (*res)(i, j) = sum;
    }
  }
  return res;
}

/* @ingroup Algebra
 *  @brief Weighted matrix multiplication
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
Matrix* weightedMultLeftTranspose( Matrix const& A, Matrix const& B, Vector const& weights)
{
#ifdef STK_DEBUG
  if (A.rangeVe() != B.rangeVe())
    throw runtime_error("In Matrix* multTranspose(A, B) "
                             "A.rangeVe() != B.rangeVe()");
#endif
  // create result
  Matrix* res = new Matrix(A.rangeHo(), B.rangeHo());
  // indexes
  const Integer first = A.firstRow(), last = A.lastRow();
  const Integer firstRow = A.firstCol(), lastRow = A.lastCol();
  const Integer firstCol = B.firstCol(), lastCol = B.lastCol();
  //
  for (Integer i=firstRow; i<=lastRow; i++)
  {
    for (Integer j=firstCol; j<=lastCol; j++)
    {
      Real sum = 0.0;
      for (Integer k = first; k<=last; k++)
        sum += weights[k] * A(k, i) * B(k, j);
      (*res)(i, j) = sum;
    }
  }
  return res;
}

/* @ingroup Algebra
 *  @brief Weighted matrix multiplication
 *
 *  Perform the matrix product \f$ A'WA \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @return a pointer on the result as a MatrixSquare
 **/
MatrixSquare* weightedMultLeftTranspose( Matrix const& A, Vector const& weights)
{
  // create result
  MatrixSquare* res = new MatrixSquare(A.rangeHo());
  // indexes
  const Integer first = A.firstCol(), last = A.lastCol();
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  //
  for (Integer i=first; i<=last; i++)
  {
    // diagonal
    Real sum = 0.0;
    for (Integer k = firstRow; k<=lastRow; k++)
      sum += weights[k] * A(k, i) * A(k, i);
    (*res)(i, i) = sum;
    // outside diagonal
    for (Integer j=first; j<i; j++)
    {
      Real sum = 0.0;
      for (Integer k = firstRow; k<=lastRow; k++)
        sum += weights[k] * A(k, i) * A(k, j);
      (*res)(j, i) = ((*res)(i, j) = sum);
    }
  }
  return res;
}

} // Namespace STK

