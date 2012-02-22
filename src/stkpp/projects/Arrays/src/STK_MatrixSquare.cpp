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
 * Purpose:  Define the MatrixSquare class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_MatrixSquare.cpp
  * @brief In this file, we implement the MatrixSquare class.
 **/

#include "../include/STK_MatrixSquare.h"

namespace STK
{
/* Default constructor with rangeHo_=(1:0) and rangeVe_=(1:0).
 *  @param I range of the Rows and Cols
 **/
MatrixSquare::MatrixSquare( Range const& I)
                          : Matrix(I, I)
{ ;}

/* constructor with rangeHo_and rangeVe_ givens,
 *  initialization with a constant.
 *  @param I range of the Rows and Cols
 *  @param v initial value of the container
 **/
MatrixSquare::MatrixSquare( Range const& I, Real const& v)
                          : Matrix(I, I, v)
{ ;}

/* Copy constructor
 *  @param T the container to copy
 *  @param ref true if T is wrapped
 **/
MatrixSquare::MatrixSquare( const MatrixSquare &T, bool ref)
                          : Matrix(T, ref)
{ ;}

  /* constructor by reference, ref_=1.                                    */
MatrixSquare::MatrixSquare( const Matrix& T, Range const& I)
                          : Matrix(T, I, I)
{ ;}


/* destructor.                                                         */
MatrixSquare::~MatrixSquare() { ;}

} // namespace STK
