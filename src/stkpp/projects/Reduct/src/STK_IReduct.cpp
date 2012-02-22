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
 * Project:  stkpp::AAModels
 * created on: 17 avr. 2010
 * Purpose:  Abstract class for the computation of the Index in the
 * AAM models.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_IReduct.cpp In this file we implement the interface
 *  base class IReduct.
 **/

#include "../include/STK_IReduct.h"

namespace STK
{
/*
 * Constructor.
 * @param data the input data set
 */
IReduct::IReduct( Matrix const& data)
                : IRunnerConstRef<Matrix>(data)
                , dim_(0)
                , p_reduced_(0)
{;}

/*
 * Destructor
 */
IReduct::~IReduct()
{ if (p_reduced_) delete p_reduced_;}

} // namespace STK

