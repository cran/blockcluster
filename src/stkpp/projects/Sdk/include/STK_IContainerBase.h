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
 * Project:  stkpp::Sdk::TContainer
 * Purpose:  Define the Interface base class for all containers classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_IContainerBase.h
 *  @brief In this file we define the interface base class IContainerBase.
 **/

#ifndef STK_ICONTAINERBASE_H
#define STK_ICONTAINERBASE_H

#include <vector>

#include "STK_IRecursiveTemplate.h"
#include "../../STKernel/include/STK_Range.h"

namespace STK
{
/** @ingroup TContainer
 *  @brief Interface base class for all containers implementing
 *  the curious recursive template paradigm.
 * 
 *  Implement the curious recursive template paradigm : the template
 *  parameter @c Container is the type of the leaf class that
 *  implements @c IContainerBase. A constant reference/pointer
 *  on the derived type can be obtained using the methods
 *  @c asLeaf and @c asPtrLeaf.
 *
 *  The default constructor of derived containers is of the form
 *  @code
 *    MyContainer1D c( Range const& r1);
 *    MyContainer2D c( Range const& r1, Range const& r2);
 *    MyContainer3D c( Range const& r1, Range const& r2, Range const& r3);
 *    ...
 *  @endcode
 *  up to the fourth dimension.
 *
 *  A container with the same type and the same dimensions can be constructed
 *  using the @c clone() pure virtual method defined in IRecursiveTemplate
 *  class.
 **/
template <class Container>
class IContainerBase : public IRecursiveTemplate<Container>
{
  public:
    /** Intrinsic dimension of the container : 1D, 2D, 3D or 4D. More
     *  dimensions are not allowed. */
    enum Dimension
    {
      _1D_ = 1,
      _2D_ = 2,
      _3D_ = 3,
      _4D_ = 4
    };

  protected:
    /** Default Constructor.
     *  The Leaf class have to give the dimension and to furnish a constructor
     *  with a vector of range as input parameter.
     *
     * @param dim the dimension of the container.
     **/
    IContainerBase( Dimension dim)
                  : IRecursiveTemplate<Container>()
                  , ranges_(dim)
                  , dim_(dim)
    { ;}

    /** Virtual Destructor. */
    virtual ~IContainerBase() { ;}

  public:
    /** Get dimension of the container.
     *  @return the dimension of the container
     **/
    inline Dimension const& dim() const { return dim_;}

  protected:
    /** Array of the ranges of the containers. */
    std::vector<Range const* > ranges_;

  private:
    /** Intrinsic dimension of the container :: 1D, 2D, 3D or 4D. More
     *  dimensions are not allowed. The dimension cannot be mdoified.
     **/
    const Dimension dim_;
};

} // namespace STK

#endif
// STK_ICONTAINERBASE_H
