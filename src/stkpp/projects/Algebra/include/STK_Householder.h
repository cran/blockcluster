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
 * Purpose:  Define Householder methods for Algebra classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Householder.h
 *  @brief In this file we define methods used by the Algebra classes.
 **/
 
#ifndef STK_HOUSEHOLDER_H
#define STK_HOUSEHOLDER_H

// Matrix class
#include "../../Arrays/include/STK_Matrix.h"

// for normInf
#include "STK_LinAlgebra1D.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute the Householder vector v of a vector x.
 * 
 *  Given a vector x, compute the vector v of the matrix of Householder
 *  \f$ P=I-2vv'/(v'v)  \f$ such that \f$ Px = v1 e_1 \f$. 
 *  The vector v is of the form : \f$ (1,x_2/s,...,x_n/s)' \f$ 
 *  and is stored in x. The value 1 is skipped and
 *  \f$ \beta = -2/(v'v) \f$ is stored in front of v.
 *  The method return the value v1.
 * 
 *  @param x the vector to rotate, it is overwritten by v
 **/
template <class TContainer1D>
Real house(ITContainer1D<Real, TContainer1D>& x)
{
  // compute L^{\infty} norm of X
  Real scale  = normInf(x);
  // first and last index fo the essential Householder vector
  Integer first = x.first()+1, last = x.last();
  // result and norm2 of X
  Real v1, norm2 = 0.0;
  // normalize the vector 
  if (scale)  // if not 0.0
  {
    for (Integer i=first; i<=last; i++)
    { x[i] /= scale; norm2 += x[i]*x[i];}
  }
  // check if the lower part is significative
  if (norm2 < Arithmetic<Real>::epsilon())
  {
    v1 = x.front(); x.front() = 0.0; // beta = 0.0
  }
  else
  {
    Real s, aux1 = x.front() / scale;
    // compute v1 = P_v X and beta of the Householder vector
    v1 =  (norm2 = sign(aux1, sqrt(aux1*aux1+norm2))) * scale;
    // compute and save beta
    x.front() = (s = aux1-norm2)/norm2;
    // comp v and save it
    for (Integer i=first; i<=last; i++) x[i] /= s;
  }
  return v1;
}

/** @ingroup Algebra
 *  @brief dot product with a Householder vector.
 * 
 *  Scalar product of a TContainer1D with a Householder vector
 *  d = < x,v>. The first composant of a Householder vector is 1.0
 *  @param x first vector
 *  @param v the Householder vector
 **/
template< class TContainer1D_1
        , class TContainer1D_2
        >
Real dotHouse( const ITContainer1D<Real, TContainer1D_1>& x
             , const ITContainer1D<Real, TContainer1D_2>& v
             )
{
  // first and last index fo the essential Householder vector
  const Integer first = v.first()+1, last = v.last();
  // compute the product
  Real sum = x[first-1] /* *1.0 */;
  for (Integer i=first; i<=last; i++) sum += x[i] * v[i];
  // return <x,v>
  return(sum);
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder vector.
 * 
 *  Perform a left multiplication of the matrix M with a Householder
 *  matrix \f$ H=I+beta vv' \f$. Overwrite M with HM.
 *  @param M the matrix to multiply (input/output)
 *  @param v the Householder vector (input)
 **/
template < class TContainer2D, class TContainer1D>
void leftHouseholder( ITContainer2D< Real, TContainer2D> const& M
                    , ITContainer1D< Real, TContainer1D> const& v
                    )
{
  // get beta
  Real beta = v.front();
  if (beta)
  {
    // get range of the Householder vector
    Range range_ve = v.range();
    // Multiplication of the cols by P=I+beta vv'
    for (Integer j=M.firstCol(); j<=M.lastCol(); j++)
    {
      // a ref on the jth column of M
      typename TContainer2D::TContainerVe Mj(M.asLeaf(), range_ve, j);
      // Computation of aux=beta* <v,M^j>
      Real aux =  dotHouse( Mj, v) * beta;
      // updating row X.first()
      Mj.front() += aux;
      // essential range of v
      const Integer first =  v.first()+1, last = v.last();
      // Computation of M^j + beta <v,M^j>  v = M^j + aux v
      for (Integer i=first; i<=last; i++)
        Mj[i] +=  v[i] * aux;
    }
  }
}

/** @ingroup Algebra
 *  @brief right multiplication by a Householder vector.
 * 
 *  Perform a right multiplication of the matrix M with a Householder
 *  matrix \f$ H=I+beta vv' \f$. Overwrite M with MH.
 * 
 *  @param M the matrix to multiply (input/output)
 *  @param v the Householder vector (input)
 **/
template < class TContainer2D, class TContainer1D>
void rightHouseholder( ITContainer2D< Real, TContainer2D> const& M
                     , ITContainer1D< Real, TContainer1D> const& v
                     )
{
  // get beta
  Real beta = v.front();
  if (beta)
  {
    // Multiplication of the cols by P=I+beta vv'
    for (Integer i=M.firstRow(); i<=M.lastRow(); i++)
    {
      // a ref on the ith row of M
      typename TContainer2D::TContainerHo Mi(M.asLeaf(), v.range(), i);
      // Computation of aux=beta* <v,M_i>
      Real aux =  dotHouse( Mi, v) * beta;
      // updating column X.first()
      Mi.front() += aux;
      // essential range of v
      const Integer first =  v.first()+1, last = v.last();
      // Computation of M_i + beta <v,M_i>  v = M_i + aux v'
      for (Integer i=first; i<=last; i++)
        Mi[i] +=  v[i] * aux;
    }
  }
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder Matrix.
 * 
 * Perform a left multiplication of the Matrix M with a Householder
 * Marix H. M <- HM with H = I + WZ'. The Householder vectors are
 * stored in the columns of H.
 * 
 * @param M the matrix to multiply
 * @param H the Householder Matrix
 **/
template < class TContainer2D>
void leftHouseholder( ITContainer2D< Real, TContainer2D> const& M
                    , Matrix const& H
                    )
{
  // compute the number of iterations
  Integer first = H.firstCol(), last = min( H.lastCol(), H.lastRow());
  // get range of the first Householder vector
  Range range_ve(last, H.lastRow());
  // iterations
  for (Integer j=last; j>=first; j--)
  {
    // apply left Householder vector to M
    leftHouseholder(M, H(range_ve, j));
    // decrease range of the Householder vector
    range_ve.decFirst();
  }
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder Matrix.
 * 
 * Perform a right multiplication of the matrix M with a Householder
 * Marix H. M <- MP with H = I + WZ'. The Householder vectors are
 * stored in the rows of H.
 * 
 * @param M the Matrix to multiply
 * @param H the Householder Matrix
 **/
template < class TContainer2D>
void rightHouseholder( ITContainer2D< Real, TContainer2D> const& M
                     , Matrix const& H
                     )
{
  // compute the number of iterations
  Integer first = H.firstCol(), last = min( H.lastCol(), H.lastRow());
  // get range of the first Householder vector
  Range range_ve(last, H.lastRow());
  // iterations
  for (Integer j=last; j>=first; j--)
  {
    // apply left Householder vector to M
    rightHouseholder(M, H(range_ve, j));
    // decrease range of the Householder vector
    range_ve.decFirst();
  }
}

} // namespace STK

#endif /*STK_HOUSEHOLDER_H*/
