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
 * Purpose:  Implement the Qr Class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Qr.cpp
 *  @brief In this file we implement the Qr Class (QR decomposition).
 **/
 
#include "../include/STK_Qr.h"

#include "../include/STK_Householder.h"
#include "../include/STK_Givens.h"

/* Templated Expression handling classes.                             */
#include "../include/STK_TOperator.h"
#include "../include/STK_TExpAlgebra.h"

namespace STK
{

/* Constructor */
Qr::Qr( const Matrix &A, bool ref)
      : Q_(A, ref)      // Creating Q
      , R_()            // Creating R
{ run();}

/* Computing the QR decomposition of the matrix Q_.                   */
void Qr::run()
{
  if (Q_.empty())     // if the container is empty
  {
    ncolr_ =0;
    ncolq_ =0;
    nrowq_ =0;
    compq_  = true;  // Q_ is computed

    return;
  }

  Q_.shift(1,1);      // translate the beg to 1
  ncolr_ = Q_.sizeHo(); // Number of cols of R
  ncolq_ = Q_.sizeHo(); // Number of column of Q
  nrowq_ = Q_.sizeVe(); // Number of rows of Q

  compq_  = false;    // Q_ is not computed

  // compute QR decomposition
  qr();
}


/* Computation of the QR decomposition                                */
void Qr::qr()
{
  Integer niter = min(nrowq_,ncolr_);   // number of iterations
  R_.resize(nrowq_, ncolr_);
  R_ = 0.0;                   // initialize to 0.0
  

  /* Main loop.                                                       */
  for (Integer iter=1, iter1=2; iter<=niter; iter++, iter1++)
  { 
    // A ref on the row of the matrix R_
    Point  Rrow0(R_, Range(iter, ncolr_), iter);
    // A ref on the row of the matrix R_
    Point  Qrow0(Q_, Range(iter, ncolq_), iter);
    // A ref on the column iter of the matrix Q_ : will contain the
    // Householder vector
    Vector u(Q_, Range(iter, nrowq_), iter);
    // compute the Householder vector of the current col
    Rrow0[iter] = house(u);
    // get beta
    Real beta = u.front();
    if (beta)
    {
      // ref on the essential part of the Householder vector
      Vector v(u, Range(iter1, nrowq_));
      // Apply Householder to next cols
      for (Integer j=iter1; j<=ncolr_; j++)
      {
        // Auxiliary data
        Real aux;
        // a ref on the jth column of Q_
        Vector Z(Q_, Range(iter1, nrowq_), j);
        // save the  current row of R_
        Rrow0[j] = Qrow0[j] + (aux = ( Qrow0[j] + dot(v, Z)) * beta);
        // update the next cols of Q_ 
        Z += v * aux;
      }
    }
    else
    {
      // just copy the row iter in R_
      for (Integer j=iter1; j<=ncolr_; j++)
      {
        Rrow0[j] = Qrow0[j];
      }
    }
  }
}


/* Computation of Q. */
void Qr::compQ()
{
  // if Q_ is computed yet
  if (compq_) return;
  // number of non zero cols of Q_  
  Integer ncol  = min(nrowq_, ncolq_);
  // Q_ is square
  if (nrowq_ < ncolq_)
  {
     Q_.popBackCols(ncolq_-nrowq_);
  }
  else
  { 
    Q_.pushBackCols(nrowq_-ncolq_);
    // initialization of the remaining cols of Q_ to 0.0
    Q_[Range(ncol+1,nrowq_)] = 0.0;
  }
  // the number of col_ is equal to the number of row
  ncolq_ = nrowq_;
  // Computation of Q_.
  // compute added col
  for (Integer iter=ncolq_; iter> ncol; --iter)
  { Q_(iter, iter) = 1.0;}
  // compute other cols
  for (Integer iter=ncol, iter1=ncol+1; iter>=1; iter--, iter1--)
  {
    // Get beta and test
    Real beta = Q_(iter,iter);
    if (beta)
    {
      // ref of the row iter
      Point P(Q_, Range(iter,ncolq_), iter);
      // ref of the column iter
      Vector X(Q_, Range(iter1,nrowq_), iter);
      // Update the cols from iter+1 to ncol
      for (Integer j=iter1; j<=ncolq_; j++)
      { Real aux;
        Vector Y(Q_, Range(iter1,nrowq_), j); // ref on the column j
        // Q_(iter, j) = beta * X'Y
        P[j] = (aux = dot( X, Y) * beta);
        // Q^j += aux * Q^iter
        Y += X * aux;
      }
      P[iter] = 1.0 + beta;
      // Q^iter *= beta
      X *= beta;
    }
    else // Q^iter = identity
    {
      Q_(iter, iter) = 1.0;
      Q_(Range(iter1,nrowq_), iter) = 0.0;
    }
    // update the column iter
    Q_(Range(1,iter-1), iter) = 0.0;
  }
  // Q_ is now computed
  compq_ = true;
}


/* Destructeur de la classe Qr                                        */
Qr::~Qr() { ;}

/* clear Q_ and R_.                                                   */    
void Qr::clear()
{
  Q_.clear();
  R_.clear();
}


/* Operator = : overwrite the Qr with S.                              */
Qr& Qr::operator=(const Qr& S)
{
  ncolr_ = S.ncolr_;       //< Number of cols of R actually computed
  ncolq_ = S.ncolq_;       //< Number of cols used by Q
  nrowq_ = S.nrowq_;       //< Number of rows used by Q
    
  compq_ = S.compq_;       //< Is Q computed ?

  Q_ = S.Q_;               //< Matrix V
  R_ = S.R_;               //< Singular values

  return *this;
}

/* Delete the jth column and update the QR decomposition : default
 * is the last col
 **/    
void Qr::popBackCols(Integer const& n)
{
  // delete n cols
  R_.popBackCols(n);
  ncolr_ -= n;
}

void Qr::eraseCol(Integer const& pos)
{
#ifdef STK_BOUNDS_CHECK
  if (pos < R_.firstCol())
  { throw out_of_range("Qr::eraseCol(pos) "
                       "pos < R_.firstCol()");
  }
  if (R_.lastCol() < pos)
  { throw out_of_range("Qr::eraseCol(pos) "
                       "R_.lastCol() < pos");
  }
#endif
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // compute the number of iteration for updating to zeroed
  Integer niter = R_.firstCol()-1+min(R_.sizeVe(), R_.sizeHo());
  // Zeroed the unwanted elements (z)
  for (Integer iter = pos+1; iter<=niter; iter++)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    R_(iter-1, iter) = compGivens( R_(iter-1, iter)
                                 , R_(iter, iter)
                                 , cosinus
                                 , sinus
                                 );
    R_(iter, iter)   = 0.0;
    // if necessary update R_ and Q_
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_[Range(iter+1, R_.lastCol())], true);
      // Update the next rows (iter1:ncolr_) of R_
      leftGivens (Rsub, iter-1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iter-1, iter, cosinus, sinus);
    }
  }
  // erase the column pos
  R_.eraseCols(pos);
  ncolr_--;
  
  // update the range of the remaining cols of the container
  R_.update(Range(pos, min(R_.lastRow(), R_.lastCol())));
}


/* Adding the last column and update the QR decomposition.               */    
void Qr::pushBackCol(Vector const& T)
{
  // check conditions
  if (T.range() != Q_.rangeVe())
  { throw runtime_error("Qr::pushBackCol(T) "
                        "T.range() != Q_.rangeVe()");
  }
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  R_.pushBackCols();
  ncolr_++;
  // Create an auxiliary container
  Vector Rncolr(Q_.rangeVe());
  // Multipliate T by Q'and put the result in Rncolr
  for (Integer j=Q_.firstCol(); j<=Q_.lastCol(); j++)
    Rncolr[j] = dot(Q_[j], T);
  // update Q_
  for (Integer iter = ncolq_-1, iter1 = ncolq_; iter>=ncolr_; iter--, iter1--)
  { 
    Real sinus, cosinus, y = Rncolr[iter], z = Rncolr[iter1] ;
    // compute the Givens rotition
    Rncolr[iter]  = compGivens(y, z, cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    {
      // Update the cols of Q_
      rightGivens(Q_, iter, iter1, cosinus, sinus);
    }
  }
  // update R_
  R_[ncolr_] = Rncolr[R_.compRangeVe(ncolr_)];
}


/* Adding the jth column and update the QR decomposition.               */
void Qr::insertCol(Vector const& T, Integer const& pos)
{
#ifdef STK_BOUNDS_CHECK
  if (pos < 1)
  { throw out_of_range("Qr::insertCol(T, pos) "
                            "j < 1");
  }
  if (ncolr_ < pos)
  { throw out_of_range("Qr::insertCol(T, pos) "
                            "ncolr_ < pos");
  }
#endif
#ifdef STK_DEBUG
  if (T.range() != Q_.rangeVe())
  { throw runtime_error("Qr::insertCol(T, pos) "
                             "T.range() != Q_.rangeVe()");
  }
#endif
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  R_.insertCols(pos);
  ncolr_++;
  // update the range of the remaining cols of R_
  R_.update(Range(pos+1, min(R_.lastRow(), R_.lastCol())));
  for (Integer i=pos+1; i<= min(R_.lastRow(), R_.lastCol()); ++i)
    R_(i,i) = 0.0;
  // A ref on the last column of R_
  Vector Rpos(Q_.rangeVe());
  // Multipliate T by Q'
  // we cannot use mult as we are using ncolq_ columns.
  for (Integer j=Q_.firstCol(); j<=Q_.lastCol(); j++)
    Rpos[j] = dot(Q_[j],T);
  // Zeroed the unwanted elements
  for (Integer iter= ncolq_-1, iter1= ncolq_; iter>=pos; iter--, iter1--)
  { 
    Real sinus, cosinus, y = Rpos[iter], z = Rpos[iter1] ;
    // compute the Givens rotation
    Rpos[iter]  = compGivens(y, z, cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_[Range(iter1, R_.lastCol())], true);
      // Update the next rows (iter1:ncolr_) of R_
      leftGivens( Rsub, iter, iter1, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iter, iter1, cosinus, sinus);
    }
  }
  // update R_
  R_[pos] = Rpos[R_.compRangeVe(pos)];
  R_.update(pos);
}

} // Namespace STK

