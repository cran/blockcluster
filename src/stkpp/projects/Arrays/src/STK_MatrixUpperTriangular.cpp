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
 * Purpose:  Define the MatrixUpperTriangular class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_MatrixUpperTriangular.cpp
 *  @brief In this file we implement the MatrixUpperTriangular class.
 **/

#include "../include/STK_MatrixUpperTriangular.h"

namespace STK
{
/* Default constructor
 * @param I range of the Rows
 * @param J range of the Cols
 **/
MatrixUpperTriangular::MatrixUpperTriangular( Range const& I, Range const& J)
                                            : _ITArray2DType(I, J)
                                            , defaultConst_(Real())
{
  // initialize vertically the container
  this->initializeCols(J);
}

/* constructor with rbeg, rend, cbeg and cend specified,
 *  initialization with a constant.
 *  @param I range of the Rows
 *  @param J range of the Cols
 *  @param v initial value of the container
 **/
MatrixUpperTriangular::MatrixUpperTriangular( Range const& I
                                            , Range const& J
                                            , Real const& v
                                            )
                                            : _ITArray2DType(I, J)
                                            , defaultConst_(Real())
{
  // initialize vertically the container
  this->initializeCols(J);
  // initialize with v
  for (Integer j = J.first(); j <= J.last(); j++)
  {
    Real* p(this->data(j));
    Integer beg(this->rangeCols_[j].first());
    Integer end(this->rangeCols_[j].last());

    for (Integer i = beg; i <= end; i++) p[i] = v;
  }
}

/* Copy constructor
 *  @param T the container to copy
 *  @param ref true if T is wrapped
 **/
MatrixUpperTriangular::MatrixUpperTriangular( const MatrixUpperTriangular &T
                                            , bool ref
                                            )
                                            : _ITArray2DType(T, ref)
                                            , defaultConst_(Real())
{
  if (!ref)
  {
    // initialize vertically the container
    this->initializeCols(T.rangeHo());
    // initialize with T
    for (Integer j=T.firstCol(); j<=T.lastCol(); j++)
    {
      Real* p(this->data(j));
      const Real* pt= T.data(j);

      for (Integer i=T.firstRow(); i<=T.lastRow(); i++) p[i]= pt[i];
    }
  }
}

/* constructor by reference, ref_=1.
 *  @param T the container to wrap
 *  @param I range of the Rows to wrap
 *  @param J range of the Cols to wrap
 **/
MatrixUpperTriangular::MatrixUpperTriangular( const _ITArray2DType& T
                                            , Range const& I
                                            , Range const& J
                                            )
                                            : _ITArray2DType(T, I, J)
                                            , defaultConst_(Real())
{ ;}

/* Wrapper constructor Contruct a reference container.
 *  @param q pointer on the data
 *  @param I range of the  Rows to wrap
 *  @param J range of the Cols to wrap
 **/
MatrixUpperTriangular::MatrixUpperTriangular( Real** q
                                            , Range const& I
                                            , Range const& J
                                            )
                                            : _ITArray2DType(q, I, J)
                                            , defaultConst_(Real())
{ ;}

/* virtual destructor : use destructor of Array2D.                            */
MatrixUpperTriangular::~MatrixUpperTriangular()
{ ;}

/* operator = : overwrite the MatrixUpperTriangular with T.
 *  We resize the object if this and T does not have the same size
 *  but if they have the same size, we don't modify the range
 *  of the object.
 *  @param T the container to copy
 **/
MatrixUpperTriangular& MatrixUpperTriangular::operator=(const MatrixUpperTriangular &T)
{
  // Resize if necessary.
  if ( (this->sizeVe() != T.sizeVe())||(this->sizeHo() != T.sizeHo()))
  {
    this->initialize(T.rangeVe(), T.rangeHo());
  }
  // coopy without overlapping
  if (T.firstCol()>=this->firstCol())
  {
    for ( Integer jt=T.firstCol(), j=this->firstCol()
        ; jt<=T.lastCol()
        ; j++, jt++)
    {
      Real* p(this->data(j));
      const Real* pt= T.data(jt);
      Integer beg(this->rangeCols_[j].first());
      Integer end(this->rangeCols_[j].last());
      Integer tbeg(T.rangeCols_[j].first());

      for (Integer i=beg, it=tbeg; i<=end; i++, it++)
        p[i]= pt[it];
    }
    return *this;
  }
  // T.firstCol()<this->firstCol()
  for ( Integer jt=T.lastCol(), j=this->lastCol()
      ; jt>=T.firstCol()
      ; j--, jt--)
  {
    Real* p(this->data(j));
    const Real* pt= T.data(jt);
    Integer beg(this->rangeCols_[j].first());
    Integer end(this->rangeCols_[j].last());
    Integer tbeg(T.rangeCols_[j].first());

    for (Integer i=beg, it=tbeg; i<=end; i++, it++) p[i]= pt[it];
  }
  return *this;
}

/** @ingroup Arrays
 *  ostream for MatrixUpperTriangular.
 *  @param s the output stream
 *  @param V the MatrixUpperTriangular to write
 **/
ostream& operator<<(ostream& s, MatrixUpperTriangular const& V)
{ return out2D(s,V);}

} // namespace STK
