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
 * Project: stkpp::STatistiK::Statdesc
 * Purpose: Perform the usual transformation on data set.
 * Author:  Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Stat_Transform.h
 *  @brief In this file we implement the main transformation on data set.
 **/

#ifndef STK_STAT_TRANSFORM_H
#define STK_STAT_TRANSFORM_H

#include "../../STKernel/include/STK_Arithmetic.h"
#include "../../STKernel/include/STK_Misc.h"

#include "../../Sdk/include/STK_ITContainer2D.h"

#include "STK_Stat_UnivariateReal.h"

namespace STK
{
namespace Stat
{

/** Compute the mean of the variables in the container V and center V.
 *  @param V the container with the Data
 *  @param mean the Vector of the means
 **/
template < class TContainerHo, class TContainer2D >
void center( ITContainer2D< Real, TContainer2D>&  V
           , ITContainer1D< Real, TContainerHo>& mean
           )
{
  mean.resize(V.rangeHo());
  // get dimensions
  const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    // compute mean
    Real mu = Stat::mean(V[j]);
    // save current mean
    mean[j] = (mu = Arithmetic<Real>::isFinite(mu) ?  mu : 0.);
    // translate data
    for (Integer i= firstRow; i<= lastRow; i++)
    { V(i,j) -= mu;}
  }
}

/** Compute the mean of the variables in the container V and center V.
 *  @param V the container with the Data
 *  @param W the Vector of the weights
 *  @param mean the Vector of the means
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
void center( ITContainer2D< Real, TContainer2D>& V
           , ITContainer1D< Real, TContainerVe> const& W
           , ITContainer1D< Real, TContainerHo>& mean
           )
{
#ifdef STK_DEBUG
  if (!V.rangeVe().isIn(W.range()))
    throw runtime_error("In center(V, W, mean) "
                             "V.range() not include in w.range()");
#endif
  // create result
  mean.resize(V.rangeHo());
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    // compute mean
    Real mu = Stat::mean<TContainerVe>(V[j], W);
    // save current mean
    mean[j] = (mu = Arithmetic<Real>::isFinite(mu) ?  mu : 0.);
    // translate data
    for (Integer i= firstRow; i<= lastRow; i++)
    { V(i,j) -= mu;}
  }
}

/** Compute the mean and the variance of the variable V and standardize it.
 *  @param V the container with the Data
 *  @param mean the Vector of the means
 *  @param std the Vector of the standard deviation
 **/
template < class TContainerHo, class TContainer2D >
void standardize( ITContainer2D< Real, TContainer2D>& V
                , ITContainer1D< Real, TContainerHo>& mean
                , ITContainer1D< Real, TContainerHo>& std
               )
{
    // create result
    mean.resize(V.rangeHo());
    std.resize(V.rangeHo());
    // center
    center(V, mean);
    // get dimensions
    const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
    const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
    for (Integer j= firstVar; j<= lastVar; j++)
    {
      // compute standard deviation
      Real dev = Stat::variance(V[j]);
      // take square root and save result
      std[j] = (dev = Arithmetic<Real>::isFinite(dev) ?  sqrt((double)dev) : 0.);
      // standardize data
      if (dev)
      {
        for (Integer i= firstRow; i<= lastRow; i++)
        { V(i,j) /= dev;}
      }
    }
}

/** Compute the mean and the variance of the variable V and standardize it.
 *  @param V the container with the Data
 *  @param W the Vector of the weights
 *  @param mean the Vector of the means
 *  @param std the Vector of the standard deviation
 **/
template < class TContainerHo, class TContainerVe, class TContainer2D >
void standardize( ITContainer2D< Real, TContainer2D>& V
                , ITContainer1D< Real, TContainerVe> const& W
                , ITContainer1D< Real, TContainerHo>& mean
                , ITContainer1D< Real, TContainerHo>& std
               )
{
#ifdef STK_DEBUG
  if (W.range() != V.rangeVe())
    throw runtime_error("In standardize(V, W, mean) "
                             "W.range() != V.rangeVe()");
#endif
  // create result
  mean.resize(V.rangeHo());
  std.resize(V.rangeHo());
  // center
  center(V, W, mean);
  // get dimensions
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    // compute standard deviation
    Real dev = Stat::variance(V[j], W);
    // take square root and save result
    std[j] = (dev = Arithmetic<Real>::isFinite(dev) ?  sqrt((double)dev) : 0.);
    // standardize data
    if (dev)
    {
      for (Integer i= firstRow; i<= lastRow; i++)
      { V(i,j) /= dev;}
    }
  }
}

/** Add the mean of the centered variables in the container V.
 *  @param V the container with the data
 *  @param mean the Vector of the means
 **/
template < class TContainerHo, class TContainer2D >
void decenter( ITContainer2D< Real, TContainer2D>&  V
             , ITContainer1D< Real, TContainerHo> const& mean
             )
{
  if (V.rangeHo() != mean.range())
    throw runtime_error(_T("Error in Stat::decenter(V, mean): "
                             "ranges are not the same."));
  // get dimensions
  const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    // translate data
    Real mu = mean[j];
    for (Integer i= firstRow; i<= lastRow; i++)
    { V(i,j) += mu;}
  }
}

/** dilate the standardized variable V.
 *  @param V the container with the Data
 *  @param std the Vector of the standard deviation
 **/
template < class TContainerHo, class TContainer2D >
void destandardize( ITContainer2D< Real, TContainer2D>& V
                  , ITContainer1D< Real, TContainerHo> const& std
                  )
{
    if (V.rangeHo() != std.range())
      throw runtime_error(_T("Error in Stat::decenter(V, std): "
                               "ranges are not the same."));
  // get dimensions
  const Integer  firstRow = V.firstRow(), lastRow = V.lastRow();
  const Integer  firstVar = V.firstCol(), lastVar = V.lastCol();
  for (Integer j= firstVar; j<= lastVar; j++)
  {
    // dilate
    if (std[j])
    {
      Real dev = std[j];
      // get dimensions
      for (Integer i= firstRow; i<= lastRow; i++)
      { V(i,j) *= dev;}
    }
  }
}

/** dilate the standardized variable V and add the mean.
 *  @param V the container with the Data
 *  @param mean the Vector of the means
 *  @param std the Vector of the standard deviation
 **/
template < class TContainerHo, class TContainer2D >
void destandardize( ITContainer2D< Real, TContainer2D>& V
                  , ITContainer1D< Real, TContainerHo> const& mean
                  , ITContainer1D< Real, TContainerHo> const& std
               )
{
  // dilate
  destandardize(V, std);
  // decenter
  decenter(V, mean);
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_TRANSFORM_H*/
