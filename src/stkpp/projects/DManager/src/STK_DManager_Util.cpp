/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2010  Serge Iovleff

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
 * created on: 12 oct. 2010
 * Purpose:  useful methods for processing strings and i/o streams..
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_DManager_Util.cpp
 *  @brief In this file we implement various methods useful for processing
 *  strings and i/o streams in the DManager project.
 **/
#include "../include/STK_DManager_Util.h"

namespace STK
{

namespace DManager
{

/*  @ingroup DManager
 *  @brief check if a string represent a boolean.
 *
 *  A String is a boolean if it is written "TRUE" or "FALSE". There is no need
 *  to use upper case.
 *  @param str the string to check
 *  @return @c true if the String i a boolean, @c false otherwise.
 **/
bool checkStringToBoolean( String const& str)
{
  // is it TRUE ?
  if (toUpperString(str) == _T("TRUE")) { return true;}
  // is it FALSE ?
  if (toUpperString(str) == _T("FALSE")) { return true;}
  // not a bolean string
  return false;
}

/* @ingroup DManager
 *  @brief check if a string represent a boolean.
 *
 *  A String is a boolean if it is written "TRUE" or "FALSE". There is no need
 *  to use upper case.
 *  @param str the string to check
 *  @return @c true if the String is a boolean, @c false otherwise.
 **/
bool StringToBoolean( String const& str)
{
  // is it TRUE ?
  if (toUpperString(str) == _T("TRUE")) { return true;}
  // if it's not true, it's false
  return false;
}


/* remove all occurrences of the char @c c at the beginning and the end
 *  of the string.
 */
void removeCharBeforeAndAfter( String & str, Char c )
{
  // erase first whitespaces
  str.erase(0, str.find_first_not_of(c));
  // erase remaining whitespaces
  size_t found =str.find_last_not_of(c);
  if (found != str.npos)
    str.erase(found+1);
  else
    str.clear(); // str is all whitespace
}

/* Get the current field from the stream.
 *  @param is the stream to treat
 *  @param delimiter the delimiter of the current field
 **/
String getField( istream& is, Char delimiter)
{
  String strbuff;
  std::getline( is, strbuff, delimiter);
  removeCharBeforeAndAfter(strbuff, CHAR_BLANK);
  removeCharBeforeAndAfter(strbuff, CHAR_TAB);
  return strbuff;
}


} // namespace DManager

} // namespace STK


