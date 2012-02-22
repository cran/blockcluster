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
 * Project:  Base
 * Purpose:  Define the fundamental type NotAvailable.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_NotAvailable.cpp
 *  @brief In this file we implement the fundamental type
 *  @c NotAvailable.
 **/

#include "../include/STK_NotAvailable.h"
#include "../include/STK_String_Util.h"

namespace STK
{
/*  Overloading of the ostream << for the type NotAvailable.
 **/
ostream& operator << (ostream& os, const NotAvailable& output)
{ return (os <<  STRING_NA);}

/*  Overloading of the istream >> for the type NotAvailable.
 **/
istream& operator >> (istream& is, Proxy<NotAvailable>& input)
{
  // in all case input is a NA object
  input = Arithmetic<NotAvailable>::NA();
  // get current file position
  std::ios::pos_type pos = is.tellg();

  // Create a String buffer
  Char Lbuff[STRING_NA_SIZE+1];

  // try to read a NA String
  is.getline(Lbuff, STRING_NA_SIZE+1);

  // if we don't get a NA String, rewind stream
  if (!(STRING_NA.compare(Lbuff) == 0))
  { is.seekg(pos); }
  return is;
}

} // namespace STK
