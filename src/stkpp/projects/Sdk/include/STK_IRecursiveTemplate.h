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
 * created on: 7 août 2011
 * Purpose:  Create an interface base class for class using the Recursive
 * template paradigm.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IRecursiveTemplate.h
 *  @brief In this file we define the Interface base class IRecursiveTemplate.
 **/

#ifndef STK_IRECURSIVETEMPLATE_H
#define STK_IRECURSIVETEMPLATE_H

namespace STK
{

/** @ingroup Sdk
 *  @brief Interface base class for all classes implementing the curious
 *  recursive template paradigm.
 *
 *  Implement the curious recursive template paradigm : the template
 *  parameter @c Leaf is the type of the leaf class that
 *  implements @c IRecursiveTemplate. A constant reference/pointer
 *  on the derived type can be obtained using the methods
 *  @c asLeaf and @c asPtrLeaf.
 **/
template<class Leaf>
class IRecursiveTemplate
{
  public:
    /** constructor. */
    IRecursiveTemplate() {;}

    /** copy constructor. */
    IRecursiveTemplate(IRecursiveTemplate const& leaf) { ;}

    /** destructor. */
    virtual ~IRecursiveTemplate() {;}

    /** static cast : return a reference of this with a cast to the
     * derived class.
     * This allow to delegate public methods to derived classes.
     * @return a reference on this in the derived type
     **/
    inline Leaf& asLeaf()
    { return static_cast<Leaf&>(*this); }

    /** static cast : return a const reference of this with a cast to
     * the derived class.
     * This allow to delegate public method to derived classes.
     * @return a constant reference on this in the derived type
     **/
    inline Leaf const& asLeaf() const
    { return static_cast<Leaf const&>(*this); }

    /** static cast : return a ptr on a @c Leaf of this with a cast to the
     * derived class.
     * This allow to delegate public methods to derived classes.
     * @return a pointer on this in the derived type
     **/
    inline Leaf* asPtrLeaf()
    { return static_cast<Leaf*>(this); }

    /** static cast : return a ptr on a constant @c Leaf of this with a cast to
     * the derived class.
     * This allow to delegate public method to derived classes.
     * @return a constant pointer on this in the derived type
     **/
    inline Leaf const* asPtrLeaf() const
    { return static_cast<Leaf const*>(this); }

    /** create a leaf using the copy constructor of the Leaf class.
     *  @return A pointer on a cloned Leaf of @c this
     **/
    inline Leaf* clone() const
    { return new Leaf(this->asLeaf());}
};

} // namespace STK

#endif /* STK_IRECURSIVETEMPLATE_H */
