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
 * Purpose:  Define the Abstract Variable class
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IVariable.h
 *  @brief Define the Interface base class of all types of Variables.
 **/

#ifndef STK_IVARIABLE_H
#define STK_IVARIABLE_H

#include "STK_DManager_Util.h"

#include "../../STKernel/include/STK_String.h"
#include "../../STKernel/include/STK_String_Util.h"
#include "../../Sdk/include/STK_IContainer1D.h"

namespace STK
{
// forward definition
template <class TYPE>
class Variable;

/** @ingroup DManager
  * @brief IVariable is an Interface class for all Variables classes.
  *
  * A Variable have a name and a type, the implementation have to choose
  * some one-dimensional container class derived from ITContainer1D
  * for storing the data.
  *
  * The pure virtual function defined in this class are the one needed by
  * the DataFrame container.
 **/
class IVariable : virtual public IContainer1D
{
  protected:
    /** Id Type of the variable. */
    const IdType type_;

    /** Name of the variable. */
    String name_;

  protected:
    /** Default constructor                                                 */
    IVariable( const IdType &type, String const& name)
             : type_(type)
             , name_(name)
    { ;}

    /** Copy  constructor                                                   */
    IVariable( const IVariable& V)
             : type_(V.type_)
             , name_(V.name_)
    { ;}

  public:
    /** type Variable of String (used for exporting the variable).    */
    typedef Variable<String> _VarStringType;

    /** destructor.                                                         */
    virtual ~IVariable() { ;}

    /** Operator = : overwrite the IVariable with T.                  */
    IVariable& operator=(const IVariable &V)
    {
      name_ = V.name_;  // copy name_ of the variable
      return *this;
    }

    /** Return a default name for variable of the form : prefix + num.
     **/
    static inline String giveName( Integer const& num    = 0
                                 , const String  &prefix = STRING_VAR
                                 )
   { return (prefix+typeToString<Integer> (num));}

    /** Get the type of the variable.                                 */
    inline const IdType& getType() const
    { return type_;}

    /** Get the name of the variable.                                 */
    inline String const& name() const
    { return name_;}

    /** Set a default name for variable of the form : prefix + num.   */
    inline void setName( Integer const& num     = 0
                       , const String  &prefix  = STRING_VAR
                )
    { name_ = (prefix+typeToString<Integer> (num));}

    /** Set the name of the variable.                                 */
    inline void setName( String const& name)
    { name_ = name;}

    /** push back n NA values.
     *  @param n number of NA values to add
     **/
    virtual void pushBackNAValues(Integer const& n=1) = 0;

    /** overwrite the IVariable by converting the strings
     *  contained in V into the Type. Give the number of success.
     *  @param V Variable of String
     *  @param f io flags
     *  @return number of successful conversion
     **/
    virtual Integer importString( const _VarStringType& V
                                , std::ios_base& (*f)(std::ios_base&) = std::dec
                                ) = 0;

    /** Overwrite the variable V by converting the data into strings.
     *  @param V Variable of String
     *  @param f io flags
     **/
    virtual const IVariable& exportString( _VarStringType& V
                                         , std::ios_base& (*f)(std::ios_base&)
                                           = std::dec
                                         ) const =0;

    /** Operator << : overwrite the IVariable by converting the strings
     *  contained in V into the Type.
     *  @param V the Variable of string to import
     **/
    virtual IVariable& operator<<( const _VarStringType& V) =0;

    /** Operator >> : convert the Variable V into strings.
     *  @param V Variable of String
     **/
    virtual const IVariable& operator>>( _VarStringType& V) const =0;

    /** clone return a ptr on a copy of the object.
     *  @param ref true if we don't want to copy the data
     **/
    virtual IVariable* clone( bool ref = false) const = 0;
};

} // Namespace STK

#endif //STK_IVARIABLE_H
