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
 * Project:  stkpp::Sdk::TContainer
 * created on: 26 nov. 2011
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ITContainer2D.h
 *  @brief In this file we define the Interface class ITContainer2D.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ITCONTAINER2D_H
#define STK_ITCONTAINER2D_H

#include "STK_IContainer2D.h"
#include "STK_ITContainer1D.h"
#include "../../STKernel/include/STK_Exceptions.h"

namespace STK
{
template <class TYPE, class TContainer2D>
class ContainerTraits;

/** @ingroup TContainer
 *  @brief Interface class for homogeneous 2D containers.
 *
 * Use the curious recursive template paradigm : the template
 * parameter @c TContainer2D is the name of the class that
 * implements the interface ITContainer2D.
 * For example
 * @code
 * template<class TYPE>
 * class TContainer2D : public ITContainer2D< TYPE
 *                                          , TContainerHo
 *                                          , TContainerVe
 *                                          , TContainer2D<TYPE>
 *                                          >
 * {...}
 * @endcode
 *
 * A @c TContainerVe is a class that allow to access to (or a part of) a
 * column of the @c TContainer2D class
 *
 * A @c TContainerHo is a class that allow to access to (or a part of) a
 * a row of the @c TContainer2D class.
 *
 * The pseudo virtual function defined in this interface and to implement
 * in the derived class have the following definitions:
 * @code
 *   TYPE& elt( Integer const& i, Integer const& j);
 *   TYPE const& elt( Integer const& i, Integer const& j) const;
 *   TContainerVe col( Integer const& j) const;
 *   TContainerVe col( Range const& I, Integer const& j) const;
 *   TContainerHo row( Integer const& i) const;
 *   TContainerHo row( Integer const& i, Range const& J) const;
 *   TContainer2D sub( Range const& I, Range const& J) const;
 * @endcode
 * while the pure virtual function (inherited from IContainer2D) to implement
 * are:
 * @code
 *    virtual void shift( Integer const& rbeg, Integer const& cbeg);
 *    virtual void eraseCols(Integer const& pos, Integer const& n = 1);
 *    virtual void pushBackCols( Integer const& n =1);
 *    virtual void popBackCols( Integer const& n =1);
 *    virtual void eraseRows(Integer const& pos, Integer const& n=1);
 *    virtual void pushBackRows( Integer const& n =1);
 *    virtual void popBackRows( Integer const& n =1);
 * @endcode
 *
 * ITContainer2D is a general interface that can be used as a contract in the
 * class using two dimensional arrays.
 **/
template <class TYPE, class TContainer2D>
class ITContainer2D : public IContainerBase<TContainer2D>, public IContainer2D
{
  protected:
    /** Default constructor. Default values are rangeHo=(1:0) and rangeVe=(1:0).
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    ITContainer2D( Range const& I = Range(), Range const& J = Range())
                 : IContainerBase<TContainer2D>(IContainerBase<TContainer2D>::_2D_)
                 , IContainer2D(I, J)
    { this->ranges_.push_back(&rangeVe_); this->ranges_.push_back(&rangeHo_);}

    /** Copy constructor.
     *  @param T the container to copy
     **/
    ITContainer2D( const ITContainer2D &T)
                 : IContainerBase<TContainer2D>(IContainerBase<TContainer2D>::_2D_)
                 , IContainer2D(T)
    { this->ranges_.push_back(&rangeVe_); this->ranges_.push_back(&rangeHo_);}

  public:
    /** Type of datra stored by the container. */
    typedef TYPE Type;
    /** type of the Vertical Container (column) of the two dimensional
     *  container. */
    typedef typename ContainerTraits<TYPE, TContainer2D>::TContainerVe TContainerVe;
    /** type of the Horizontal Container (row) of the two dimensional
     *  container. */
    typedef typename ContainerTraits<TYPE, TContainer2D>::TContainerHo TContainerHo;

    /** virtual destructor. */
    virtual ~ITContainer2D() { ;}

    /** access to one element.
     *  @param i index of the row
     *  @param j index of the column
     *  @return the element (i,j) of the 2D container.
     **/
    inline TYPE& elt(Integer const& i, Integer const& j)
    { return this->asLeaf().elt(i,j);}

    /** access to a constant element.
     *  @param i index of the row
     *  @param j index of the column
     *  @return a constant element (i,j) of the 2D container.
     **/
    inline TYPE const& elt(Integer const& i, Integer const& j) const
    { return this->asLeaf().elt(i,j);}

    /** Operator () : access to one element.
     *  @param i index of the row
     *  @param j index of the column
     *  @return the element (i,j) of the 2D container.
     **/
    inline TYPE& operator()(Integer const& i, Integer const& j)
    { return this->asLeaf().elt(i,j);}

    /** Operator () : access to one constant element.
     *  @param i index of the row
     *  @param j index of the column
     *  @return a constant element (i,j) of the 2D container.
     **/
    inline TYPE const& operator()(Integer const& i, Integer const& j) const
    { return this->asLeaf().elt(i,j);}

    /** Operator [] : access to one column.
     *  @param j index of the column
     *  @return a Vertical container containing the column @c j of the Container
     **/
    inline TContainerVe operator[](Integer const& j) const
    { return this->asLeaf().col(j);}

    /** access to one column.
     *  @param j index of the column
     *  @return a Vertical container containing the column @c j of the Container
     **/
    inline TContainerVe col(Integer const& j) const
    { return this->asLeaf().col(j);}

    /** Operator () : access to many elements of a column.
     *  @param I range of the index of the rows
     *  @param j index of the col
     *  @return a Vertical container containing the column @c j of the Container
     *  in the range @c I
     **/
    inline TContainerVe operator()(Range const& I, Integer const& j) const
    { return this->asLeaf().col(I, j);}

    /** Operator () : access to one row.
     *  @param i index of the row
     *  @return an Horizontal container containing the row @c i of the Container
     **/
    inline TContainerHo operator()(Integer const& i) const
    { return this->asLeaf().row(i);}

    /** Operator () : access to one row.
     *  @param i index of the row
     *  @return an Horizontal container containing the row @c i of the Container
     **/
    inline TContainerHo row(Integer const& i) const
    { return this->asLeaf().row(i);}

    /** Operator () : access to many elements of a row.
     *  @param i index of the row
     *  @param J index of the col
     *  @return an Horizontal container containing the row @c i of the Container
     *  in the range @c J
     **/
    inline TContainerHo operator()(Integer const& i, Range const& J) const
    { return this->asLeaf().row(i, J);}

    /** Operator () : access to a sub-array.
     *  @param I range of the index of the rows
     *  @param J range of the index of the cols
     *  @return a 2D container containing the Container in the range @c I, @c J
     **/
    inline TContainer2D operator()(Range const& I, Range const& J) const
    { return this->asLeaf().sub(I, J);}

    /** Operator () : access to many rows.
     *  @param I range of the index of the rows
     *  @return a 2D container containing the Container in the vertical range
     *  @c I
     **/
    inline TContainer2D operator()(Range const& I) const
    { return this->asLeaf().sub(I, this->rangeHo());}

    /** Operator [] : access to many columns.
     *  @param J range of the index of the cols
     *  @return a 2D container containing the Container in the Horizontal range
     *  @c J
     **/
    inline TContainer2D operator[](Range const& J) const
    { return this->asLeaf().sub(this->rangeVe(), J);}

    /** return safely the element (i, j).
     *  @param i index of the row
     *  @param j index of the col
     *  @return the element (i,j) of the 2D container.
     **/
    TYPE& at(Integer const& i, Integer const& j)
    {
      if (this->firstRow() > i)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->firstRow() > i");
      }
      if (this->lastRow() < i)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->lastRow() < i");
      }
      if (this->firstCol() > j)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->firstCol() > j");
      }
      if (this->lastCol() < j)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->lastCol() < j");
      }
      return this->asLeaf().elt(i, j);
    }

    /** return safely the constant element (i, j).
     *  @param i index of the row
     *  @param j index of the col
     *  @return the constant element (i,j) of the 2D container.
     **/
    TYPE const& at(Integer const& i, Integer const& j) const
    {
      // check bounds
      if (this->firstRow() > i)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->firstRow() > i");
      }
      if (this->lastRow() < i)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->lastRow() < i");
      }
      if (this->firstCol() > j)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->firstCol() > j");
      }
      if (this->lastCol() < j)
      { throw out_of_range("ITContainer2DBase::at(i, j) "
                                "this->lastCol() < j");
      }
      // return element
      return this->asLeaf().elt(i, j);
    }

    /** STL compatibility : return a sub-Array.
     *  @param I range of the index of the rows
     *  @param J range of the index of the cols
     **/
    inline TContainer2D at(Range const& I, Range const& J) const
    {
      if (this->firstRow() > I.first())
      { throw out_of_range("ITContainer2D::at(I, J) "
                                "this->firstRow() > I.first()");
      }
      if (this->lastRow() < I.last())
      { throw out_of_range("ITContainer2D::at(I, J) "
                                "this->lastRow() < I.last()");
      }
      if (this->firstCol() > J.first())
      { throw out_of_range("ITContainer2D::at(I, J) "
                                "this->firstCol() > J.first()");
      }
      if (this->lastCol() < J.last())
      { throw out_of_range("ITContainer2D::at(I, J) "
                                "this->lastCol() < J.last()");
      }
      return this->asLeaf().sub(I, J);
    }

    /** STL compatibility : return a part of the column j in the range I.
     *  @param I range of the index of the rows
     *  @param j index of the col
     **/
    inline TContainerVe at(Range const& I, Integer const& j) const
    {
      if (this->firstRow() > I.first())
      { throw out_of_range("TContainer2D::at(I, j) "
                           "this->firstRow() > I.first()");
      }
      if (this->lastRow() < I.last())
      { throw out_of_range("TContainer2D::at(I, j) "
                           "this->lastRow() < I.last()");
      }
      if (this->firstCol() > j)
      { throw out_of_range("TContainer2D::at(I, j) "
                           "this->firstCol() > j");
      }
      if (this->lastCol() < j)
      { throw out_of_range("TContainer2D::at(I, j) "
                           "this->lastCol() < j");
      }
      return this->asLeaf().col(I, j);
    }

    /** STL compatibility : return a part of the row i in the range J.
     *  @param i index of the row
     *  @param J range of the index of the cols
     **/
    inline TContainerHo at(Integer const& i, Range const& J) const
    {
      if (this->firstRow() > i)
      { throw out_of_range("TContainer2D::at(i, J) "
                                "this->firstRow() > i");
      }
      if (this->lastRow() < i)
      { throw out_of_range("TContainer2D::at(i, J) "
                                "this->lastRow() < i");
      }
      if (this->firstCol() > J.first())
      { throw out_of_range("TContainer2D::at(i, J) "
                                "this->firstCol() > J.first()");
      }
      if (this->lastCol() < J.last())
      { throw out_of_range("TContainer2D::at(i, J) "
                                "this->lastCol() < J.last()");
      }
      return this->asLeaf().row(J, i);
    }

    /** STL compatibility : return the column j.
     *  @param j index of the col
     **/
    inline TContainerVe atCol(Integer const& j) const
    {
      if (this->firstCol() > j)
      { throw out_of_range("TContainer2D::atCol(j) "
                                "this->firstCol() > j");
      }
      if (this->lastCol() < j)
      { throw out_of_range("TContainer2D::atCol(j) "
                                "this->lastCol() < j");
      }
      return this->asLeaf().col(this->rangeVe(), j);
    }

    /** STL compatibility : return the Container2D in column range J.
     *  @param J range of the index of the cols
    **/
    inline TContainer2D atCol(Range const& J) const
    {
      if (this->firstCol() > J.first())
      { throw out_of_range("TContainer2D::atCol(J) "
                                "this->firstCol() > J.first()");
      }
      if (this->lastCol() < J.last())
      { throw out_of_range("TContainer2D::atCol(J) "
                                "this->lastCol() < J.last()");
      }
      return this->asLeaf().sub(this->rangeVe(), J);
    }

    /** STL compatibility : return the row i.
     *  @param i the index of the row
     **/
    inline TContainerHo atRow(Integer const& i) const
    {
      if (this->firstRow() > i)
      { throw out_of_range("TContainer2D::atRow(i) "
                                "this->firstRow() > i");
      }
      if (this->lastRow() < i)
      { throw out_of_range("TContainer2D::atRow(i) "
                                "this->lastRow() < i");
      }
      return this->asLeaf().row(i, this->rangeHo());
    }

    /** STL compatibility : return the Container2D in the row range I.
     *  @param I range of the index of the rows
     **/
    inline TContainer2D atRow(Range const& I) const
    {
      if (this->firstRow() > I.first())
      { throw out_of_range("TContainer2D::atRow(I) "
                                "this->firstRow() > I.first()");
      }
      if (this->lastRow() < I.last())
      { throw out_of_range("TContainer2D::atRow(I) "
                                "this->lastRow() < I.last()");
      }
      return this->asLeaf().sub(I, this->rangeHo());
    }

    /** push back a column to the container with value V.
     *  @param V the values to add to the end of the container
     **/
    template<class Container1D>
    inline void pushBackCol(ITContainer1D<TYPE, Container1D> const& V)
    {
      const Integer firstRow = V.first(), lastRow = V.last();
      // check if the container is empty
      if (this->empty())
      {
        this->resize(V.range(), Range(1));
        for (Integer i=firstRow; i<=lastRow; i++)
          (*this)(i, 1) = V[i];
        return;
      }
#ifdef STK_BOUNDS_CHECK
      if (V.range() != this->rangeVe())
      { throw runtime_error("TContainer2D::pushBackCol(V) "
                                 "V.range() != this->rangeVe()");
      }
#endif
      // if the container is not empty we add a column and copy V inside
      this->asLeaf().pushBackCols();
      const Integer lastCol = this->lastCol();
      for (Integer i=firstRow; i<=lastRow; i++)
        (*this)(i, lastCol) = V[i];
    }
};

} // namespace STK

#endif
// STK_ITCONTAINER2D_H
