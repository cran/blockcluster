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
 * Purpose:  Define the GramScmidt Class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_GramSchmidt.h
 *  @brief In this file we define the GramSchidt method.
 *
 **/
 
#ifndef STK_GRAMSCHMIDT_H
#define STK_GRAMSCHMIDT_H

// Container classes
#include "../../Sdk/include/STK_ITContainer2D.h"
#include "STK_LinAlgebra1D.h"

namespace STK
{

/** @ingroup Algebra
 *  @brief The gramSchmidt method perform the Gram Schmidt orthonormalization
 *  of a ITContainer2D of Real.
 *  @param A the container to orthnormalize
 **/
template < class TContainer2D>
void gramSchmidt( ITContainer2D< Real, TContainer2D>& A)
{
  // get dimensions
  const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
  const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
  // orthonormalize
  for (Integer j= firstCol; j<= lastCol; j++)
  {
    for( Integer i= firstCol; i < j; i++)
    {
      const Real dotij = dot(A.col(i), A.col(j));
      for (Integer k= firstRow; k<= lastRow; k++)
      {
        A(k,j) -= A(k,i) * dotij;
      }
    }
    // normalize
    const Real norm = normTwo(A.col(j));
    if (norm)
    {
      for (Integer k= firstRow; k<= lastRow; k++)
      {
        A(k,j) /= norm;
      }
    }
  }
}

} // Namespace STK

#endif // STK_GRAMSCHMIDT_H

