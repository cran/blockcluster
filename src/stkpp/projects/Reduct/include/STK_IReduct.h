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
 * Project:  stkpp::Reduct
 * created on: 17 avr. 2010
 * Purpose:  Abstract class for the computation of the Index in the SLAAM.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IReduct.h In this file we define the interface base
 *  class IReduct.
 **/

#ifndef STK_IREDUCT_H
#define STK_IREDUCT_H

#include "../../Sdk/include/STK_IRunner.h"
#include "../../Arrays/include/STK_Matrix.h"

namespace STK
{

/** @ingroup Reduct
 *  @brief A IReduct is an interface base class for dimension
 *  reduction techniques.
 *
 * The class receive a matrix in input of size (n,p).
 * - n is the number of samples,
 * - p is the number of variables.
 * The observations can be weighted.
 *
 * The derived class will compute a @em reduced data set of dimension (n,d).
 */
class IReduct : public IRunnerConstRef<Matrix>
{
  protected:
    /** Constructor.
     *  @param data The data set to reduce.
     * */
    IReduct( Matrix const& data);

  public:
    /** virtual destructor  */
    virtual ~IReduct();

    /** get the number of dimension.
     *  @return The number of dimension computed
     **/
    inline Integer const& dim() const { return dim_;}
    /** get a pointer on the reduced data set
     * @return a constant pointer on the data set
     **/
    inline Matrix* p_reduced() const { return p_reduced_; }
    /** set the number of dimension.
     *  @param dim the number of dimension to set
     **/
    inline void setDimension( Integer const& dim) { dim_ = dim;}

  protected:
    /** dimension of the reduced data set */
    Integer dim_;
    /** The reduced data set. */
    Matrix* p_reduced_;
};

} // namespace STK

#endif /* STK_ILINEARREDUCT_H */
