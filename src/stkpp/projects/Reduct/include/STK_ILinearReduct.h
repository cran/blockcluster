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
 * Purpose:  Abstract class for the computation of the Index in the SLAAM.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ILinearReduct.h In this file we define the interface base class
 *  ILinearReduct.
 **/

#ifndef STK_ILINEAREDUCT_H
#define STK_ILINEAREDUCT_H

#include "../../Arrays/include/STK_Vector.h"
#include "../../Arrays/include/STK_Matrix.h"

#include "STK_IReduct.h"

namespace STK
{
/** @ingroup AAModels
 *  @brief A ILinearReduct is an interface base class for criteria
 *  to maximize in order to find the axis of a linear reduction.
 *
 * The class receive a matrix in input of size (n,p).
 * - n is the number of samples,
 * - p is the number of variables.
 * The n samples can be weighted.
 *
 * The class compute the optimal axis (stored in the @c axis_ )
 * attribute and the @c projected data set (stored in the @c p_reduct_ atribute
 * of the base class IReduct) when the user use the virtual method @c run(nbAxis)
 * (not weighted observations) or @c run(weights, nbAxis) (weighted observations).
 *
 * The Matrix axis_ is computed by maximizing some criteria defined in
 * derived classes. It is constructed using the pure virtual functions:
 * @code
 *  virtual void maximizeIndex() =0;
 *  virtual void maximizeIndex() =0;
 * @endcode
 */
class ILinearReduct : public IReduct
{
  protected:
    /** Constructor.
     * @param data the data set to reduce.
     **/
    ILinearReduct( Matrix const& data);

  public:
    /** virtual destructor  */
    virtual ~ILinearReduct();

    /** Compute the Index.
     */
    bool run();
    /**
     * Compute the weighted index.
     * @param weights the weights to used
     */
    bool run(Vector const& weights);

    /** get the number of axis.
     *  @return The number of axis computed or to compute
     **/
    inline Integer const& nbAxis() const { return dim_;}
    /** get the values of the index of each axis
     * @return a constant reference Vector of the value of the Index
     **/
    inline Vector const& indexValues() const { return index_values_; }
    /** get the axis
     * @return a constant reference Matrix of the axis
     **/
    inline Matrix const& axis() const { return axis_; }

  protected:
    /** The values of the index for each axis. */
    Vector index_values_;
    /** The computed axis. */
    Matrix axis_;
    /** Find the axis by maximizing the Index. */
    virtual void maximizeIndex() =0;
    /** Find the axis by maximizing the weighted Index.
     *  @param weights the weights to used
     * */
    virtual void maximizeIndex( const Vector& weights) =0;
    /** Compute the projection of the data set on the Axis. */
    void projection();
};

} // namespace STK

#endif /* STK_ILINEARREDUCT_H */
