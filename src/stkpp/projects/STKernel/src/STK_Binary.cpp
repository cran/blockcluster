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
 * Project:  stkpp::STKernel::Base
 * Purpose:  Implement the fundamental type Binary.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Binary.cpp
 *  @brief In this file we implement the fundamental type Binary.
 **/

#include "../include/STK_Binary.h"
#include "../include/STK_String_Util.h"

namespace STK
{

/*  Overloading of the ostream << for the type Binary.
 **/
ostream& operator << (ostream& os, const Binary& output)
{ if (Arithmetic<Binary>::isNA(output))
    return (os <<  STRING_NA);
  else
    return (os << static_cast<int> (output));
}

/*  Overloading of the istream >> for the type Binary.
 **/
istream& operator >> (istream& is, Proxy<Binary>& input)
{
  // get current file position
  std::ios::pos_type pos = is.tellg();
  // try to read a discrete value
  int buff;
  // failed to read a discrete value
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
    input = Arithmetic<Binary>::NA();
    // get current position in the stream
    pos = is.tellg();

    // Create a String buffer
    Char Lbuff[STRING_NA_SIZE+1];

    // try to read a NA String
    is.getline(Lbuff, STRING_NA_SIZE+1);

    // if we don't get a NA String, rewind stream
    if (!(STRING_NA.compare(Lbuff) == 0))
    { is.seekg(pos); }
  }
  else
    input = static_cast<Binary>(buff);
  return is;
}



} // namespace STK
