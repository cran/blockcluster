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
 * Project:  stkpp::Arrays
 * Purpose:  Primary include file for Base sub-project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file Arrays.h
 *  @brief This file include all the other header files of the
 *  project Arays.
 **/

/**
 * @defgroup Arrays concrete implementation of the TContainer Interfaces
 * @brief The Arrays project provides a concrete implementation of the
 * Interfaces ITContainer1D and ITContainer2D defined in the Sdk project.
 **/

#ifndef ARRAYS_H
#define ARRAYS_H

#include "../projects/Arrays/include/STK_Point.h"
#include "../projects/Arrays/include/STK_Vector.h"
#include "../projects/Arrays/include/STK_Matrix.h"
#include "../projects/Arrays/include/STK_MatrixSquare.h"
#include "../projects/Arrays/include/STK_MatrixUpperTriangular.h"
#include "../projects/Arrays/include/STK_MatrixLowerTriangular.h"

/* 2D Array. */
#include "../projects/Arrays/include/STK_Array2D.h"
#include "../projects/Arrays/include/STK_CArray2D.h"

#include "../projects/Arrays/include/STK_Array1D.h"
#include "../projects/Arrays/include/STK_ArrayHo.h"

/* Interface base class for all Arrays. */
#include "../projects/Arrays/include/STK_AllocatorBase.h"

#endif  /* ARRAYS_H */
