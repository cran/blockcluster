/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * created on: 26 nov. 2011
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_CArrayHo.h
 *  @brief In this file .
 **/

#ifndef STK_CARRAYHO_H
#define STK_CARRAYHO_H

#include "STK_AllocatorBase.h"

namespace STK
{

/** @brief Horizontal C-like container.
 *  The CArrayHo allow to get a view of a row of a CArray2D without copying the
 *  data stored in the two-dimensional container.
 */
template < class TYPE>
class CArrayHo : public ITContainer1D<TYPE, CArrayHo<TYPE> >
               , public AllocatorBase<TYPE>
{
  public:

    CArrayHo()
    { ;}
    virtual ~CArrayHo()
    { ;}
};

} // namespace STK

#endif /* STK_CARRAYHO_H */
