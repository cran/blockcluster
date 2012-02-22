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
 * Project:  stkpp::STKernel
 * Purpose:  Primary include file for the STKernel project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STKernel.h
 *  @brief This file include all the other header files of the
 *  project STKernel.
 *
 **/

/**
 * @defgroup STKernel Kernel tools
 * @brief The STKernel project is the low-level core library that forms the
 * basis of the project. It provides data class handling for C++.
 * It contains the sub-projects:
 * - @ref Base
 **/

/**
 *  @ingroup STKernel
 *  @defgroup Base Subproject STKernel::Base
 *  In the Base subproject, we define
 * <ul>
 * <li> the standard types for handling data:
 * <ul>
 *   <li> Integer for discrete data
 *   <li> Real    for quantitative data
 *   <li> Binary  for binary data in {0,1}
 *   <li> String  for string data
 *   <li> Sign    for signed data in {-1, 1}.
 * </ul>
 * For all these types, we define a not available special value that can be
 * displayed in a transparent way using proxy classes.
 * <li> standard input and output streams for all these types.
 * <li> and miscellaneous utilities functions.
 * </ul>
 **/

#ifndef STKERNEL_H
#define STKERNEL_H

/* Arithmetic classes for fundamental types. */
#include "../projects/STKernel/include/STK_Arithmetic.h"

/* RTTI class for fundamental types. */
#include "../projects/STKernel/include/STK_IdTypeImpl.h"

/* fundamental STK Char.  */
#include "../projects/STKernel/include/STK_Char.h"

/* STK streams parametrized with Char.  */
#include "../projects/STKernel/include/STK_Stream.h"

/* STK String parametrized with Char. */
#include "../projects/STKernel/include/STK_String.h"

/* Fundamental constant of STKpp. */
#include "../projects/STKernel/include/STK_String_Util.h"

/* Miscellaneous functions on Strings. */
#include "../projects/STKernel/include/STK_String_Util.h"

/* Proxy classes for the fundamental types of STKpp.  */
#include "../projects/STKernel/include/STK_Proxy.h"

/* Fundamental types of STKpp.  */
#include "../projects/STKernel/include/STK_TypeBase.h"

/* Chrono functions.  */
#include "../projects/STKernel/include/STK_Chrono.h"

/* Miscellaneous functions. */
#include "../projects/STKernel/include/STK_Misc.h"

/* Index range. */
#include "../projects/STKernel/include/STK_Range.h"

/* Standard exceptions */
#include "../projects/STKernel/include/STK_Exceptions.h"

#endif  /* STKERNEL_H */
