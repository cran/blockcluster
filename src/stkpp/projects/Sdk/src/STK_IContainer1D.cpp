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
 * Project:  stkpp::Sdk::TContainer
 * Purpose:  Define the Interface 1D Container class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IContainer1D.cpp
 *  @brief In this file we implement the IContainer1D class.
 **/

#include "../include/STK_IContainer1D.h"

namespace STK
{

/* Default constructor */
IContainer1D::IContainer1D() { }

/* Virtual destructor. */
IContainer1D::~IContainer1D() { }

/* STL compatibility. Resize the container.
 * - call @c shift
 * - call @c pushBack if there will be more elements
 * - call @c popBack if three will be less elements
 * @param I the range to set to the container
 **/
void IContainer1D::resize(Range const& I)
{
  // check if there is something to do
  if ( range() == I) return;
  // translate beg
  shift(I.first());
  // compute number of elements to delete or add
  const Integer inc = I.last() - last();
  // adjust size of the container
  if (inc > 0) pushBack(inc);  // more elements
  else         popBack(-inc);  // less elements
}

} // namespace STK
