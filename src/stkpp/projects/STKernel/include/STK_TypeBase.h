/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004  Serge Iovleff

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
 * Project:  Base
 * Purpose:  Define the fundamental types.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_TypeBase.h
 *  @brief In this file we include all the headers of the fundamental types used
 *  internally in the STK applications. You can add any type matching your need.
 * 
 *   Actually (the 20070216) the fundamental types defined are
 * - Char (char)
 * - String     (std::string)
 * - Real      (double)
 * - Integer  (long)
 * - Binary    enum : {0, 1}
 * - Sign    enum : {-1, 1}
 * - NotAvailable enum : {'.'}
 * 
 *  If you add a type:
 *  -# Specialize the struct Arithmetic, otherwise
 *     the type will not have NA value defined.
 *  -# Specialize the struct IdTypeImpl,
 *     otherwise calling the returnType method in the class Variable
 *     will give the unknown type.
 *  -# Surdefine the operators >> and << for the input and output stream in
 *  order to handle the NA value in a transparent way.
 **/

#ifndef STK_TYPEBASE_H
#define STK_TYPEBASE_H

#include "STK_Arithmetic.h"
#include "STK_IdTypeImpl.h"

#include "STK_Char.h"
#include "STK_String.h"

#include "STK_Integer.h"
#include "STK_Real.h"

#include "STK_Binary.h"
#include "STK_Sign.h"

#include "STK_NotAvailable.h"

#endif /* STK_TYPEBASE_H */
