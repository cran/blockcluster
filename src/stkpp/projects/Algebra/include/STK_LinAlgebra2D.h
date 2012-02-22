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
 * Purpose:  Define 2D Linear Algebra methods with Real Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_LinAlgebra2D.h
 *  @brief Define Linear Algebra methods for two dimensional
 *  Containers containing Real.
 **/
 
#ifndef STK_LINALGEBRA2D_H
#define STK_LINALGEBRA2D_H

// Point and Vector classes
#include "../../Arrays/include/STK_Point.h"
#include "../../Arrays/include/STK_Vector.h"
#include "../../Arrays/include/STK_Matrix.h"
#include "../../Arrays/include/STK_MatrixSquare.h"

#include "STK_LinAlgebra1D.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief trace of a square matrix
 *
 *  Compute the trace of the matrix A.
 *
 *  @param[in] A the Matrix
 *  @return the sum of the diagonal elements
 **/
/** @ingroup Algebra
 *  @brief Matrix multiplication by a Vector
 *
 *  Perform the product of the matrix A by the Vector X, Y= AX.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @param[out] Y the result
 *  @return a pointer on the result as a Vector
 **/
template < class TContainer2D , class TContainer1D>
void mult( ITContainer2D<Real, TContainer2D> const& A
         , ITContainer1D<Real, TContainer1D> const& X
         , ITContainer1D<Real, TContainer1D>& Y
         )
{
  // create result
  Y.resize(A.rangeVe());
  // indexes
  const Integer first = A.firstRow(), last = A.lastRow();
  for (Integer i=first; i<=last; i++)
  {
    Y[i] = dot(A.row(i), X);
  }
}

/** @ingroup Algebra
 *  @brief Square Matrix multiplication by a Vector.
 *
 *  Perform the product of the square matrix A by the Vector X, Y= AX.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @param[out] Y the result
 *  @return a pointer on the result as a Vector
 **/
template <class Container1D>
void mult( MatrixSquare const& A
         , ITContainer1D<Real, Container1D> const& X
         , ITContainer1D<Real, Container1D>& Y
         )
{
  // create result
  Y.resize(A.range());
  // indexes
  const Integer first = A.first(), last = A.last();
  for (Integer i=first; i<=last; i++)
  {
    Y[i] = dot<Point, Container1D>(A.row(i), X);
  }
}

/** @ingroup Algebra
 *  @brief Matrix multiplication by a Vector
 *
 *  Perform the product of the transposed matrix A by the Vector X, Y= A'X.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @param[out] Y the result
 *  @return a pointer on the result as a Vector
 **/
template < class TContainer2D , class TContainer1D>
void multLefTranspose( ITContainer2D<Real, TContainer2D> const& A
                     , ITContainer1D<Real, TContainer1D> const& X
                     , ITContainer1D<Real, TContainer1D>& Y
                     )
{
  // create result
  Y.resize(A.rangeHo());
  // indexes
  const Integer first = A.firstCol(), last = A.lastCol();
  for (Integer j=first; j<=last; j++)
  {
    Y[j] = dot(A.col(j), X);
  }
}

/** @ingroup Algebra
 *  @brief transposed square Matrix multiplication by a Vector.
 *
 *  Perform the product of the transposed square matrix A by the Vector X, Y= AX.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @param[out] Y the result
 *  @return a pointer on the result as a Vector
 **/
template <class Container1D>
void multLeftTranspose( MatrixSquare const& A
                      , ITContainer1D<Real, Container1D> const& X
                      , ITContainer1D<Real, Container1D>& Y
                      )
{
  // create result
  Y.resize(A.range());
  // indexes
  const Integer first = A.first(), last = A.last();
  for (Integer i=first; i<=last; i++)
  {
      Y[i] = dot<Vector, Container1D>(A[i], X);
  }
}

/** @ingroup Algebra
 *  @brief Compute the infinity norm of a 2D container
 * 
 *  Return the maximal absolute value of the 2D container x
 *  \f[ s= \max_{i=1}^n |x_i| \f]
 * 
 *  @param[in] A matrix to treat.
 *  @return the infinity norm of A
 **/
template < class TContainer2D >
Real normInf( ITContainer2D<Real, TContainer2D> const& A)
{
  // indexes
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  // find maximal absolute value
  Real scale = 0.0;
  for (Integer j=firstRow; j<=lastRow; j++)
    for (Integer i=firstCol; i<=lastCol; i++)
      scale = STK::max(scale, STK::abs(A(i,j)));
  return (scale);
}

/** @ingroup Algebra
 *  @brief transpose a matrix
 * 
 *  Transpose the Matrix A and give the result in At.
 * 
 *  @param[in] A the matrix to transpose
 *  @param[out] At the transposed matrix
 **/
template < class TContainer2D1, class TContainer2D2>
void transpose( ITContainer2D< Real, TContainer2D1> const& A
              , ITContainer2D< Real, TContainer2D2>& At
              )
{
  // Resize At.
  At.resize(A.rangeHo(), A.rangeVe());

  // indexes
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  // copy each column of A in each row of At
  for (Integer j=firstRow; j<=lastRow; j++)
    for (Integer i=firstCol; i<=lastCol; i++)
      At(i, j) = A(j, i);
}

/** @ingroup Algebra
 *  @brief Transpose a Matrix
 * 
 *  Transpose a matrix and overload it with the result
 * 
 *  @param[in, out] Q the matrix to transpose
 *  @return a reference on the matrix transposed
 **/
template < class TContainer2D >
ITContainer2D< Real, TContainer2D>& transpose(ITContainer2D< Real, TContainer2D>& Q)
{
  // check if we have to create an auxiliary container
  if (Q.rangeVe() != Q.rangeVe())
  {
    // create temporary Matrix
    Matrix R;
    transpose(R, Q);
    Q = R;
    return Q;
  }
  // Square Matrix
  const Integer first = Q.firstRow(), last = Q.lastRow();
  // swap elements
  for (Integer  i=first; i<=last; i++)
  {
    for ( Integer j= first; j< last; j++)
    {
      Real aux(Q(i,j));
      Q(i,j)  = Q(j,i);
      Q(j,i) = aux;
    }
  }
  return Q;
}

/** @ingroup Algebra
 *  @brief Compute the trace of a  square Matrix
 *
 *  @param[in] A the square matrix
 *  @return the trace of the matrix
 **/
Real trace( MatrixSquare const& A);


/** @ingroup Algebra
 *  @brief Matrix by Vector multiplication [DEPRECATED]
 *
 *  Perform the product \f$ AX \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @return a pointer on the result as a Vector
 **/
template <class Container1D>
Vector* mult( Matrix const& A, ITContainer1D<Real, Container1D> const& X)
{
#ifdef STK_DEBUG
  if (A.rangeHo() != X.range())
    throw runtime_error("In Vector* mult(A, X) "
                             "A.rangeHo() != X.range()");
#endif
  // create result
  Vector* res = new Vector(A.rangeVe());
  // indexes
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
  for (Integer i=firstRow; i<=lastRow; i++)
  {
      Real sum = 0.0;
      for (Integer k = firstCol; k<=lastCol; k++)
        sum += A(i, k) * X[k];
      (*res)[i] = sum;
  }
  return res;
}

/** @ingroup Algebra
 *  @brief Matrix by Vector multiplication [DEPRECATED]
 *
 *  Perform the product \f$ A'X \f$.
 *
 * The result is created on the stack.
 *
 *  @param[in] A Left operand
 *  @param[in] X Right operand
 *  @return a pointer on the result as a Vector
 **/
template <class Container1D>
Vector* multLeftTranspose( Matrix const& A, ITContainer1D<Real, Container1D> const& X)
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

#endif
// STK_LINALGEBRA2D_H
