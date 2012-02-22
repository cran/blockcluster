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
 * Purpose:  Implement 2D Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_LinAlgebra2D.cpp
 *  @brief Implement Linear Algebra methods with Real for two dimensional
 *  containers.
 **/
 
#include "../include/STK_LinAlgebra2D.h"

namespace STK
{
/* @ingroup Algebra
 *  @brief trace of a square matrix
 *
 *  Compute the trace of the matrix A.
 *
 *  @param[in] A the Matrix
 *  @return the sum of the diagonal elements
 **/
Real trace( MatrixSquare const& A)
{
  Real sum = 0.0;
  const Integer first = A.firstCol(), last = A.lastCol();
  for (Integer k = first; k<=last; k++)
    sum += A(k, k);
  return sum;
}

/* @ingroup Algebra
 *  @brief Matrix multiplication by a Vector
 *
 *  Perform the product of the matrix A by the Vector X
 *  The matrix A and X must have the correct range.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @return a pointer on the result as a Vector
 **/
Vector* multLeftTranspose( Matrix const& A, Vector const& X)
{
#ifdef STK_DEBUG
  if (A.rangeVe() != X.range())
    throw runtime_error("In Vector* multLeftTranspose(A, X) "
                             "A.rangeHo() != X.range()");
#endif
  // create result
  Vector* res = new Vector(A.rangeHo());
  // indexes
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();

  for (Integer i=firstCol; i<=lastCol; i++)
  {
      Real sum = 0.0;
      for (Integer k = firstRow; k<=lastRow; k++)
        sum += A(k, i) * X[k];
      (*res)[i] = sum;
  }
  return res;
}

} // Namespace STK

