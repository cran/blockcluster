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

/** @file STK_IContainer1D.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array1D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ICONTAINER1D_H
#define STK_ICONTAINER1D_H

#include "../../STKernel/include/STK_Range.h"

namespace STK
{
/** @ingroup TContainer
 *  @brief Interface base class for homogeneous 1D containers.
 *
 * The IContainer1D class is an interface base class for all
 * one dimensional containers.
 *
 * The pure virtual function defined in this interface have the
 * following definition:
 * @code
 *   void insertElt(Integer const& pos, Integer const& n =1) = 0;
 *   void clear();
 *   void erase( Integer const& pos, Integer const& n=1);
 *   void shift(Integer const& beg =1) = 0;
 *   void pushBack(Integer const& n =1) = 0;
 *   void popBack( Integer const& n =1);
 * @endcode
 *
 * All these methods have to be implemented in the Container1D class.
 **/
class IContainer1D
{
  protected:
    /** Default constructor */
    IContainer1D();

    /** Virtual destructor. */
    ~IContainer1D();

  public:
    /**  get the index of the first element
     *   @return the index of the first element
     */
    virtual Integer const& first() const =0;

    /**  get the index of the last element
     *   @return the index of the last element
     */
    virtual Integer const& last() const =0;

    /**  get the number of elements
     *   @return the size of the container
     */
    virtual Integer const& size() const =0;

    /**  get the range of the container
     *   @return the range of the container
     */
    virtual Range const& range() const =0;

    /** Is there some data ?
     *  @return @c true if the container is empty, @c false otherwise
     **/
    inline bool empty() const { return range().empty();}

    /**  STL compatibility. Delete last element. */
    inline void pop_back() { erase(last());}

    /**  STL compatibility. Delete first element. */
    inline void pop_front() { erase(first());}

    /** STL compatibility. Resize the container.
     * - call @c shift
     * - call @c pushBack if there will be more elements
     * - call @c popBack if three will be less elements
     * @param I the range to set to the container
     **/
    void resize(Range const& I = Range());

//  protected:
    /** Insert n elements at the position pos of the container.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    virtual void insertElt(Integer const& pos, Integer const& n =1) = 0;

    /** clear Container from all elements and memory allocated. */
    virtual void clear() = 0;

    /** Delete n elements at the @c pos index from the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
     **/
    virtual void erase(Integer const& pos, Integer const& n=1) = 0;

    /** New first index for the object.
     *  @param beg the index of the first column to set
     **/
    virtual void shift(Integer const& beg =1) = 0;

    /** Add n elements at the end of the container.
     *  @param n number of elements to add
     **/
    virtual void pushBack(Integer const& n =1) = 0;

    /** Delete n last elements of the container.
     *  @param n number of elements to delete
     **/
    virtual void popBack(Integer const& n =1) = 0;
};

} // namespace STK

#endif
// STK_ICONTAINER1D_H
