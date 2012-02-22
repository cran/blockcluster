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
 * Purpose:  Define the Interface 1D templated Container class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ITContainer1D.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array1D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ITCONTAINER1D_H
#define STK_ITCONTAINER1D_H

#include "STK_IContainerBase.h"
#include "STK_IContainer1D.h"
#include "../../STKernel/include/STK_Exceptions.h"

namespace STK
{
/** @ingroup TContainer
 *  @brief Interface base class for homogeneous 1D containers.
 *
 * The ITContainer1D class is the templated base class for all
 * homogeneous one-dimensional containers containing element of type @c TYPE.
 *
 * Implement the curious recursive template paradigm : the template
 * parameter @c TContainer1D is the name of the class that
 * implements @c ITContainer1D. For example
 * <code>
 * template<class TYPE>
 * class TContainer1D : public ITContainer1D<TYPE, TContainer1D<TYPE> >
 * {...}
 * </code>
 *
 * The pseudo virtual function defined in this interface have the
 * following definition:
 * @code
 *   TYPE& elt(Integer const& pos);
 *   TYPE const& elt(Integer const& pos) const;
 *   TContainer1D elt(Range const& I) const;
 * @endcode
 *
 * All these methods have to be implemented in the TContainer1D class.
 **/
template <class TYPE, class TContainer1D>
class ITContainer1D : public IContainerBase<TContainer1D>
                    , virtual public IContainer1D
{
  protected:
    /** Default constructor
     *  @param I : the range of the container
     **/
    ITContainer1D( Range const& I = Range())
                 : IContainerBase<TContainer1D>(IContainerBase<TContainer1D>::_1D_)
                 , range_(I)
    { this->ranges_.push_back(&range_);}

    /** Copy constructor
     *  @param T the container to copy
     **/
    ITContainer1D( const ITContainer1D& T)
                 : IContainerBase<TContainer1D>(IContainerBase<TContainer1D>::_1D_)
                 , range_(T.range_)
    { this->ranges_.push_back(&range_);}

  public:
    /** Type of datra stored by the container. */
    typedef TYPE Type;
    /** Virtual destructor. */
    virtual ~ITContainer1D() { ;}

    /**  Index of the first element
     *   @return the index of the first element
     */
    inline Integer const& first() const { return range_.first();}

    /**  Index of the last element
     *   @return the index of the last element
     */
    inline Integer const& last() const { return range_.last();}

    /**  Number of elements
     *   @return the size of the the container
     */
    inline Integer const& size() const { return range_.size();}

    /**  range of the container
     *   @return the range of the container
     */
    inline Range const& range() const { return range_;}

    /** Access to one element.
     *  @param pos position of the element
     *  @return a reference on the pos-th element
     **/
    inline TYPE& elt(Integer const& pos)
    { return this->asLeaf().elt(pos);}

    /** Access to one element.
     *  @param pos position of the element
     *  @return a constant reference on the pos-th element
     **/
    inline TYPE const& elt(Integer const& pos) const
    { return this->asLeaf().elt(pos);}

    /** Access to many elements.
     *  @param I the range of the elements
     *  @return A container with the elements in the range I of this
     **/
    inline TContainer1D elt(Range const& I) const
    { return this->asLeaf().elt(I);}

    /** Operator [] : access to one element.
     * @param pos position of the element
     *  @return a reference on the pos-th element
     **/
    inline TYPE& operator[](Integer const& pos)
    { return this->asLeaf().elt(pos);}

    /** Operator [] : access to one element const.
     *  @param pos position of the element
     *  @return a constant reference on the pos-th element
     **/
    inline TYPE const& operator[](Integer const& pos) const
    { return this->asLeaf().elt(pos);}

    /** Operator [] : access to many elements.
     *  @param I the range of the elements
     *  @return a container with the elements in the range I of this
     **/
    inline TContainer1D operator[](Range const& I) const
    { return this->asLeaf().elt(I);}

    /** Swap two elements of the container using default constructor
     *  of the class @c TYPE.
     *  @param pos1 position of the first element
     *  @param pos2 position of the second element
     **/
    void swap(Integer const& pos1, Integer const& pos2)
    {
#ifdef STK_BOUNDS_CHECK
      // check indices
      if (pos1<this->first())
      { throw out_of_range("ITContainer1D::swap(pos1, pos2) "
                                "pos1<this->first()");
      }
      if (pos1>this->last())
      { throw out_of_range("ITContainer1D::swap(pos1, pos2) "
                                "pos1>this->last()");
      }
      if (pos2<this->first())
      { throw out_of_range("ITContainer1D::swap(pos1, pos2) "
                                "pos2<this->first()");
      }
      if (pos2>this->last())
      { throw out_of_range("ITContainer1D::swap(pos1, pos2) "
                                "pos2>this->last()");
      }
#endif
      // swap
      TYPE aux(this->asLeaf().elt(pos1));
      this->asLeaf().elt(pos1) = this->asLeaf().elt(pos2);
      this->asLeaf().elt(pos2) = aux;
    }

    /** STL compatibility : return the element number pos.
     *  @param pos position of the element
     * @return a reference on the pos-th element of the container
     **/
    inline TYPE& at(Integer const& pos)
    {
      if ((pos<this->first()))
      {
        throw out_of_range("ITContainer1D::at(pos) "
                                "pos<this->first()");
      }
      if ((pos>this->last()))
      { throw out_of_range("ITContainer1D::at(pos) "
                                "pos>this->last()");
      }
      return this->asLeaf().elt(pos);
    }

    /** STL compatibility : return the element number pos (const).
     *  @param pos position of the element
     * @return a constant reference on the pos-th element of the container
     **/
    inline TYPE const& at(Integer const& pos) const
    {
      if (pos<this->first())
      {
        throw out_of_range("ITContainer1D::at(pos) "
                                "pos<this->first()");
      }
      if (pos>this->last())
      { throw out_of_range("ITContainer1D::at(pos) "
                                "pos>this->last()");
      }
      return this->asLeaf().elt(pos);
    }

    /** STL compatibility : return the first element.
     * @return a reference on the first element of the container
     **/
    inline TYPE& front()
    { return this->asLeaf().elt(this->first()); }

    /** STL compatibility : return the first element (const).
     * @return a constant reference on the first element of the container
     **/
    inline TYPE const& front() const
    { return this->asLeaf().elt(this->first()); }

    /** STL compatibility : return the last element.
     * @return a reference on the last element of the container
     **/
    inline TYPE& back()
    { return this->asLeaf().elt(this->last()); }

    /** STL compatibility : return the last element (const).
     * @return a constant reference on the last element of the container
     **/
    inline TYPE const& back() const
    { return this->asLeaf().elt(this->last()); }

    /** STL compatibility : append an element v.
     *  @param v value to append
     **/
    inline void push_back(TYPE const& v)
    {
      this->asLeaf().pushBack();
      this->asLeaf().elt(this->last()) = v;
    }

    /** STL compatibility : push front an element.
     *  @param v value to append
     **/
//    inline void push_front(TYPE const& v)
//    { insert(Range(this->first(), this->first()), v);}

    /** STL compatibility : insert n elements v in the range I
     *  @param I range of the element to insert
     *  @param v value to insert
     **/
//    inline void insert(Range const& I, TYPE const& v)
//    {
//      this->insertElt(I.first(), I.size());
//      for (Integer i=I.first(); i<=I.last(); i++)
//        this->asLeaf().elt(i) = v;
//    }

  protected:
    /** Range of the container. */
    Range range_;

    /** Set range of the rows of the container.
     *  @param I the range to set (default empty)
     **/
    inline void setRange(Range const& I = Range())
    { range_ = I;}

    /** increment the range of the container (can be negative).
     *  @param inc increment to apply to the range
     **/
    inline void incRange(Integer const& inc =1)
    { range_.inc(inc);}

    /** increment the beginning of the container (can be negative).
     *  @param inc the increment to apply to the beginning of the range
     **/
    inline void incFirst(Integer const& inc = 1)
    { range_.incFirst(inc);}

    /** decrement the beginning of the container.
     *  @param inc the decrement to apply to the beginning of the range
     **/
    inline void decFirst(Integer const& inc = 1)
    { range_.decFirst(inc);}

    /** increment the end of the container (can be negative).
     *  @param inc the increment to apply to the end of the range
     **/
    inline void incLast(Integer const& inc =1)
    { range_.incLast(inc);}

    /** decrement the end of the container.
     *  @param inc the decrement to apply to the end of the range
     **/
    inline void decLast(Integer const& inc =1)
    { range_.decLast(inc);}
};

} // namespace STK

#endif
// STK_ITCONTAINER1D_H
