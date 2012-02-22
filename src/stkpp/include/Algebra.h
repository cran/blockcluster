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
 * Purpose:  Primary include file for Algebra sub-project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file Algebra.h
 *  @brief This file include all the header files of the project Algebra.
 **/

/**
 *  @defgroup Algebra Algebra
 *  @brief The Algebra project provides structures, tools and methods of the
 *  usual algebra techniques.
 *
 * The main contribution of the Algebra project are
 * <ul>
 *   <li> standard Real containers like Point, Vector, Matrix, ...
 *        using containers specialization
 *   <li> templated expressions evaluation for scientific
 *        computing @see http://www.oonumerics.org/blitz/papers/
 *    <li> usual methods of linear algebra
 *    <li> QR and SVD decompositions for matrices
 *    <li> an Eigenvalue decomposition for symmetric matrices and a generalized
 *    inversion method for such matrices.
 *  </ul>
 **/


#ifndef Algebra_H
#define Algebra_H

#include "../projects/Algebra/include/STK_TExpAlgebra.h"

/* The Linear Algebra methods. */
#include "../projects/Algebra/include/STK_LinAlgebra1D.h"
#include "../projects/Algebra/include/STK_LinAlgebra2D.h"
#include "../projects/Algebra/include/STK_LinAlgebra3D.h"
#include "../projects/Algebra/include/STK_Givens.h"
#include "../projects/Algebra/include/STK_GramSchmidt.h"
#include "../projects/Algebra/include/STK_Householder.h"

/* The Algebra classes. */
#include "../projects/Algebra/include/STK_Qr.h"
#include "../projects/Algebra/include/STK_Svd.h"
#include "../projects/Algebra/include/STK_EigenvaluesSymmetric.h"

/* Functor: compute the generalized inverse of a symetric matrix. */
#include "../projects/Algebra/include/STK_GinvSymmetric.h"

#endif /*Algebra_H*/

