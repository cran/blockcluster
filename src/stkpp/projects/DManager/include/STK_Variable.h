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
 * Project:  stkpp::DManager
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Variable.h
 *  @brief Define a templated implementation of the interface class
 *  IVariable.
 **/

#ifndef STK_VARIABLE_H
#define STK_VARIABLE_H

#include "../../Arrays/include/STK_RecursiveArray1D.h"
#include "../../Arrays/include/STK_Display1D.h"

#include "STK_IVariable.h"

namespace STK
{

/** @ingroup DManager
  * @brief Variable is an implementation of the Base class IVariable
  * using The Array1D class for storing the data.
  * It implement all purely virtual methods defined in the IVariable base
  * class.
 **/
template< class TYPE>
class Variable : public IVariable
               , public RecursiveArray1D<TYPE, Variable<TYPE> >
{
  /** Type for the Base reference Class.                            */
  typedef AllocatorBase<TYPE*> _AllocatorBaseType_;
  /** Type for the IArray1DBase Class.                              */
  typedef RecursiveArray1D<TYPE, Variable<TYPE> > _RecArray1DType_;

  public:
    /** Default constructor
     *  @param I : range of the data
     *  @param name : name of the variable
     **/
    Variable( Range const& I = Range()
            , String const& name = Arithmetic<String>::NA()
            )
            : IVariable(IdTypeImpl<TYPE>::returnType(), name)
            , _RecArray1DType_(I)
    { ;}

    /** Misc constructor
     *  @param I : range of the data
     *  @param name : name of the variable
     *  @param v    : initial value
     **/
    Variable( Range const& I, TYPE const& v, String const& name = Arithmetic<String>::NA())
            : IVariable(IdTypeImpl<TYPE>::returnType(), name)
            , _RecArray1DType_(I, v)
     { ;}

    /** Copy ctor.
     *  @param V the Variable to copy
     *  @param ref true if we want to wrap V
     **/
    Variable( Variable const& V, bool ref = false)
            : IVariable(V)
            , _RecArray1DType_(V, ref)
     { ;}

    /** Reference constructor
     *  @param V the Variable to wrap
     *  @param I range of the data
     *  @param col the column of the data
     **/
    Variable( Variable const& V, Range const& I, Integer const& col)
            : IVariable(V)
            , _RecArray1DType_(V, I, col)
    { ;}

    /** Virtual destructor.
     **/
    virtual ~Variable()
    { ;}

    /** access to one element.
     *  @param pos index of the element
     **/
    inline TYPE& elt(Integer const& pos)
    { return this->data(pos);}

    /** access to one element const.
     *  @param pos index of the const element
     **/
    inline TYPE const& elt(Integer const& pos) const
    { return this->data(pos);}

    /** access to many elements.
     *  @param J the range of the elements
     **/
    inline Variable elt(Range const& J) const
    {
#ifdef STK_BOUNDS_CHECK
      if ((J.first()<this->first()))
      { throw out_of_range("Array1D::elt(J) "
                                "J.first()<this->first()");
      }
      if ((J.last()>this->last()))
      { throw out_of_range("Array1D::elt(J) "
                                "J.last()>this->last()");
      }
#endif
      return Variable(*this, J);
    }

    /** Operator = : overwrite the Variable with V.
     *  @param V the variable to copy
     * */
    inline Variable& operator=(Variable const& V)
    {
            IVariable *p3 = this;   // convert this to IVariable
      const IVariable *p4 = &V;     // convert V    to IVariable
      (*p3) = (*p4);                // use = of IVariable class

      // check size
      if (size()!=V.size()) resize(V.range());
      // copy without ovelapping.
      if (first() < V.first())
      { for (Integer i=first(), j=V.first(); j<=V.last(); i++, j++)
          this->setData(i, V.data(j));
      }
      else
      { for (Integer i=last(), j=V.last(); j>=V.first(); i--, j--)
          this->setData(i, V.data(j));
      }
      return *this;
    }

    /** Operator = : overwrite the Variable with the value v.
     *  @param v the value to set to the variable
     **/
    Variable& operator=(TYPE const& v)
    {
      for (Integer i=first(); i<=last(); i++)
        this->setData(i, v);
      return *this;
    }

    /** clone return a ptr on a copy of the Object.
     *  @param ref true if we want just a reference
     **/
    virtual Variable* clone( bool ref = false) const
    { return new Variable(*this, ref);}

    /** push back n NA values.
     *  @param n number of NA values to add
     **/
    inline void pushBackNAValues(Integer const& n=1)
    {
      this->insert(Range(last() +1, n), Arithmetic<TYPE>::NA());
    }

    /** overwrite the Variable by converting the strings
     *  contained in V into the Type. Give the number of success.
     *  @param V Variable of String
     *  @param f io flags
     *  @return number of successful conversion
     **/
    virtual Integer importString( Variable< String > const& V
                                , std::ios_base& (*f)(std::ios_base&) = std::dec
                                )
    {
      this->resize(V.range());
      this->setName(V.name());
      Integer nSuccess = V.size();
      for (Integer i=V.first(); i<=V.last(); i++)
        if (Arithmetic<String>::isNA(V[i])) // not Available
          this->at(i) = Arithmetic<TYPE>::NA();
        else
        if (!stringToType<TYPE>(this->at(i), V[i]))
          nSuccess--;
      return nSuccess;
    }

    /** Overwrite the variable V by converting the data into strings.
     *  @param V Variable of String
     *  @param f io flags
     **/
    virtual Variable const& exportString( Variable< String >& V
                                        , std::ios_base &(*f)(std::ios_base&)
                                          = std::dec
                                        ) const
    {
      V.resize(this->range());
      V.setName(this->name());
      for (Integer i=first(); i<=last(); i++)
      {
        V[i] = typeToString<TYPE>(this->at(i), f);
      }
      return *this;
    }

    /** Operator << : overwrite the IVariable by converting the strings
     *  contained in V into the Type.
     *  @param V the Variable of string to import
     **/
    virtual Variable& operator<<( Variable< String > const& V)
    {
      this->resize(V.range());
      this->setName(V.name());
      for (Integer i=V.first(); i<=V.last(); i++)
        stringToType<TYPE>(this->at(i), V[i]);
      return *this;
    }

    /** Operator >> : convert the Variable V into strings.
     *  @param V Variable of String
     **/
    virtual Variable const& operator>>(Variable< String >& V) const
    {
      V.resize(this->range());
      V.setName(this->name());
      for (Integer i=first(); i<=last(); i++)
        V[i] = typeToString<TYPE>(this->at(i));

      return *this;
    }
    /** ostream for Variable.
     **/
    friend ostream& operator<<(ostream& s, Variable<TYPE> const& V)
    {
      s << V.name() << STRING_NL;
      return out1D(s,V);
    }
};

} // Namespace STK

#endif //STK_VARIABLE_H
