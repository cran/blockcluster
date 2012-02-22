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
 * Project:  stkpp::STatistiK::Law
 * created on: 29 juil. 2011
 * Purpose:  .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Law_ITMultivariate.h
 *  @brief In this file we define the Interface base class ITMultivariate.
 **/

#ifndef STK_ITMULTIVARIATE_H
#define STK_ITMULTIVARIATE_H

#include "../../Sdk/include/STK_IRecursiveTemplate.h"
#include "../../Sdk/include/STK_ITContainer1D.h"

#include "../../STKernel/include/STK_Real.h"

#include "STK_Law_ILawBase.h"

namespace STK
{

namespace Law
{

/** @ingroup Laws
 *  @brief Interface base Class for the multivariate distributions.
 *
 *  @c TYPE if the type of data handled by the distribution.
 *  @c Container1D is the type of container containing the data set.
 */
template<class TYPE, class Container1D>
class ITMultivariate: public ILawBase
{
  protected:
    /** Constructor. */
    ITMultivariate( String const& name) : ILawBase(name) {}

  public:
    /** destructor. */
    virtual ~ITMultivariate() {}

    /** @brief compute the probability distribution function (density) of the
     * multivariate law.
     *  Give the value of the pdf at the point x.
     *  @param x the multivariate value to compute the pdf.
     *  @return the value of the pdf
     **/
    virtual Real pdf( ITContainer1D<TYPE, Container1D> const& x) const =0;

    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x the multivariate value to compute the lpdf.
     *  @return the value of the log-pdf
     **/
    virtual Real lpdf( ITContainer1D<TYPE, Container1D> const& x) const =0;

    /** @brief simulate a realization of the Multivariate Law and store the
     *  result in x.
     *  @param[out] x the simulated value.
     **/
    virtual void rand( ITContainer1D<TYPE, Container1D>& x) const =0;

    /** @brief simulate a realization of the Multivariate Law and store the
     *  result in x. In order to get a result, @c Container1D should be a
     *  @c IContainerRef container.
     *  @param x the simulated value.
     **/
    virtual void rand( ITContainer1D<TYPE, Container1D> const& x) const =0;
};

} // namespace Law

} // namespace STK

#endif /* STK_ITMULTIVARIATE_H */
