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
 * Project:  stkpp::STKernel::Base
 * Purpose:  Define miscenaleous utility functions for Strings.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_String_Util.h
 *  @brief In this file we define miscellaneous utilities
 *  functions for Strings and the global String used for handling the NA values.
 **/

#ifndef STK_STRING_UTIL_H
#define STK_STRING_UTIL_H

#include "STK_String.h"
#include "STK_Stream.h"

#include "STK_Proxy.h"

namespace STK
{
/** @ingroup Base
  * @brief Representation of a Not Available value.
  *
  * We represent a Not Available value of any type as a "." (like in
  * (SAS(R))) for the end-user.
  **/
static const String STRING_NA      = String(_T("."));

/** @ingroup Base
  * @brief Size (in number of Char) of a Not Available value.
  *
  * We represent a Not Available value of any type as a "." (like in
  * (SAS(R))) for the end-user.
  **/
static const int STRING_NA_SIZE  = 1;

/** @ingroup Base
 *  @brief convert a String to TYPE
 *
 *  This method return true if the String s could be converted into
 *  a correct TYPE t.
 *  http://www.codeguru.com/forum/showpost.php?p=678440&postcount=1
 *  http://c.developpez.com/faq/cpp/?page=strings#STRINGS_is_type
 *
 *  The conversion is successful if it does not remain
 *  Char inside the String.
 *  The operator >> have been overloaded for each base type in order to
 *  return a NA value if fail.
 *
 *  @param t The value to get from the String
 *  @param s the String to convert
 *  @param f flags
 *  @return @c true if the conversion succeed, @c false otherwise
 **/
template <class TYPE>
bool stringToType( TYPE  &t
                 , String const& s
                 , std::ios_base& (*f)(std::ios_base&) = std::dec
                 )
{
  istringstream iss(s);
  Proxy<TYPE> wrapper_t(t);
  bool flag1 = (iss >> f >> wrapper_t).fail();
  iss >> std::ws;
  // ok if the conversion success and the String is exhausted
  return ( !flag1 && iss.eof() );
}

/** @ingroup Base
 *  @brief convert a TYPE to String
 *
 *  This method return the TYPE t into a String s.
 *  @see http://www.codeguru.com/forum/showpost.php?p=678440&postcount=1
 *  @see http://c.developpez.com/faq/cpp/?page=strings#STRINGS_convertform
 *
 *  @param t The value to convert to String
 *  @param f flags
 **/
template <class TYPE>
String typeToString( TYPE const& t
                   , std::ios_base& (*f)(std::ios_base&) = std::dec
                   )
{
  ostringstream oss;
  oss << f << ConstProxy<TYPE>(t);
  return oss.str();
}

/** @ingroup Base
 *  @brief Overloading of the operator << for the type TYPE using a
 *  constant Proxy. All output stream should use a ConstProxy in
 *  a STK application. For the enumerated types, we have also to define
 *  the standard output.
 * @param os the output stream
 * @param output the value to send to the stream
 **/
template <class TYPE>
ostream& operator << (ostream& os, const ConstProxy<TYPE>& output)
{ if (Arithmetic<TYPE>::isNA(output))
    return (os <<  STRING_NA);
  else
   return (os << static_cast<TYPE const &>(output));
}


/** @ingroup Base
 *  @brief Overloading of the operator >> for the type TYPE using a
 *  Proxy. All input stream should use a Proxy in a STK application.
 *  For the enumerated and String types, we have to overload the method.
 *  Due to the instruction
 *  @code
 *   is >> buff
 *  @endcode
 *  this operator will only work for the fundamental C/C++ types. For the other
 *  types, the operator
 *  @code
 *  operator >> (istream& is, TYPE& input);
 *  @endcode
 *  have to be implemented.
 *  @param is the input stream
 *  @param input the value to get from the stream
 **/
template <class TYPE>
istream& operator >> (istream& is, Proxy<TYPE>& input)
{
  TYPE buff;
  // get current position in the stream
  typename std::ios::pos_type pos = is.tellg();

  // If the standard Conversion failed
  if ((is >> buff).fail())
  {
    // clear failbit state
    is.clear(is.rdstate() & ~std::ios::failbit);
    // clear eofbit state if necessary and rewind position
    if (is.eof())
    {
      is.seekg(pos);
      is.clear(is.rdstate() & ~std::ios::eofbit);
    }
    // in all case input is a NA object
    input = Arithmetic<TYPE>::NA();

    // get current position in the stream
    pos = is.tellg();

    // Create a String buffer
    Char Lbuff[STRING_NA_SIZE+1];

    // try to read a NA String
    is.getline(Lbuff, STRING_NA_SIZE+1);

    // if we don't get a NA String, rewind stream
    if (!(STRING_NA.compare(Lbuff) == 0))
    { is.seekg(pos);}
  }
  else
  { input = buff;}
  // return the stream
  return is;
}

/** @ingroup Base
 *  @brief Overloading of the istream >> for the type String.
 *  @param is the input stream
 *  @param input the value to get from the stream
 **/
istream& operator >> (istream& is, Proxy<String>& input);

/** @ingroup Base
 *  @brief convert the characters of the String to upper case
 *
 *  @param s The String to convert
 *  @return the string converted to upper case
 **/
String& toUpperString( String& s);

/** @ingroup Base
 *  @brief convert the characters of the String to upper case
 *
 *  @param s The String to convert
 *  @return a String in upper case
 **/
String toUpperString( String const& s);

} // namespace STK

#endif /* STK_STRING_UTIL_H */
