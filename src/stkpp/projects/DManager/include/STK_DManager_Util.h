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
 * Project:  stkpp::DManager
 * created on: 12 juil. 2011
 * Purpose:  useful methods for processing dta, strings and i/o streams.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_DManager_Util.h
 *  @brief In this file we define useful constant and method specific to the
 *  DManager project.
 **/

#ifndef STK_DMANAGER_UTIL_H
#define STK_DMANAGER_UTIL_H

#include <list>

#include "../../STKernel/include/STK_Char.h"
#include "../../STKernel/include/STK_String_Util.h"
#include "../../STKernel/include/STK_Stream.h"

namespace STK
{

/** @ingroup DManager
  * @brief Representation of Not Available String value.
  * Empty String : the empty string is also the NA value of the class String.
  **/
static const String STRING_EMPTY   = String();

/** @ingroup DManager
 * @brief Representation of a New Line String. */
static const String STRING_NL      = _T("\n");

/** @ingroup DManager
 * @brief  Representation of a blank value. */
static const String STRING_BLANK  = _T(" ");

/** @ingroup DManager
 * @brief  default prefix of a variable name. */
static const String STRING_VAR  = _T("Var");

/** @ingroup DManager
 * @brief  The char indicating the beginning of a comment in an
 *  option file.
 **/
static const Char CHAR_COMMENT = _T('#');

/** @ingroup DManager
 * @brief  The char indicating an equality in an option file. */
static const Char CHAR_EQUAL = _T('=');

/** @ingroup DManager
 * @brief  The blank space char. */
static const Char CHAR_BLANK = _T(' ');

/** @ingroup DManager
 * @brief  The tab char. */
static const Char CHAR_TAB = _T('\t');

/** @ingroup DManager
 * @brief  The default separator char in list of option. */
static const Char CHAR_SEP = _T(',');

/** @ingroup DManager
 * @brief  The open bracket char. */
static const Char CHAR_OPENBRACKET = _T('[');

/** @ingroup DManager
 * @brief  The close bracket char. */
static const Char CHAR_CLOSEBRACKET = _T(']');

namespace DManager
{
/**  @ingroup DManager
 *   @brief check if a string represent a boolean.
 *
 *  A String is a boolean if it is written "TRUE" or "FALSE". There is no need
 *  to use upper case.
 *  @param str the string to check
 *  @return @c true if the String i a boolean, @c false otherwise.
 **/
bool checkStringToBoolean( String const& str);

/** @ingroup DManager
 *  @brief convert a string to a boolean.
 *
 *  A String is a boolean if it is written "TRUE" or "FALSE". There is no need
 *  to use upper case.
 *  @param str the string to convert
 *  @return @c true if the String is "TRUE, @c false otherwise.
 **/
bool StringToBoolean( String const& str);

/** @ingroup DManager
 *  @brief remove all occurrences of the char @c c at the beginning and the end
 *  of the string @c str.
 *  @param str the string to treat
 *  @param c the character to remove before and after
 */
void removeCharBeforeAndAfter( String & str, Char c );

/**  @ingroup DManager
 *  @brief Get the current field from the input stream.
 *
 *  A field is between the current position and a delimiter or an end of line in
 *  the stream. All blank spaces and tabulations before and after the field are
 *  removed.
 *  @param is the stream to treat
 *  @param delimiter the delimiter of the current field
 *  @return the field extracted from the input stream
 **/
String getField( istream& is, Char delimiter);

/** @ingroup DManager
 * Read a list of value of type @c TYPE stored in a line
 * @param strBuffer the string with the list of value
 * @param lst the resulting list
 * @param sep the separator character
 */
template<class TYPE>
void readList(String const& strBuffer, std::list<TYPE>& lst, Char sep = CHAR_SEP)
{
  // Declare an input string stream
  istringstream instream;
  // Use strBuffer as source of input.
  instream.str(strBuffer);
  // read the line
  do
  {
    // get field
    String strbuff = getField(instream, sep);
    // check if it is a blank field
    if (strbuff.empty()) {  break;}
    // append Data to the list
    TYPE value;
    if (stringToType(value, strbuff))
      lst.push_back(value);
    // TODO: else emit warning or Exception
  }
  while(1);
}

/**  @ingroup DManager
 * @brief Write a list of value of type @c TYPE stored in a line.
 * @param os the output stream
 * @param lst the list to write
 * @param sep the separator character
 */
template<class TYPE>
void writeList( ostream& os, std::list<TYPE> const& lst, Char sep = CHAR_SEP)
{
  if (lst.empty()) return;
  typename std::list<TYPE>::const_iterator it = lst.begin();
  if (it == lst.end()) return;
  os << *it;
  it++;
  for ( ; it != lst.end(); it++)
  {
    os << sep << STRING_BLANK << *it;
  }
}


} // namespace DManager

} // namespace STK

#endif /* STK_DMANAGER_UTIL_H */
