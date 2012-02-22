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
 * Project:  stkpp::Sdk
 * created on: 14 déc. 2011
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ContainerTraits.h
 *  @brief In this file we define the main traits class we use for the
 *  stk++ Containers.
 **/


#ifndef STK_CONTAINERTRAITS_H
#define STK_CONTAINERTRAITS_H

namespace STK
{
/** The traits class ContainerTraits must be specialized for any 2D
 *  container in the interface base class ITContainer2D<TYPE, TContainer2D>.
 *  Each specialization must provide a TContainerHo typedef and a TContainerVe
 *  typedef. For example:
 *  @code
 *  template<class TYPE>
 *  struct ContainerTraits<TYPE, Array2D<TYPE> >
 *  {
 *    typedef ArrayHo<TYPE> TContainerHo;
 *    typedef Array1D<TYPE> TContainerVe;
 *  };
 *  @endcode
 */
template <class TYPE, typename TContainer>
struct ContainerTraits
{
};

} // namespace STK

#endif /* STK_CONTAINERTRAITS_H */
