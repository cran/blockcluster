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
 * Project:  IOaccess
 * Purpose:  Implementation of the classes ImportExportToCsv.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 */

/** @file STK_ExportToCsv.cpp
 *  @brief In this file we implement the classes ImportCsv and ExportToCsv.
 **/

#include "../include/STK_ExportToCsv.h"

namespace STK
{

ExportToCsv::ExportToCsv() : p_data_(0)
                           , isColNamed_(true)
{
  // create an empty ReadWriteCsv
  p_data_ = new ReadWriteCsv();
}

ExportToCsv::ExportToCsv( const DataFrame& df)
                        : p_data_(0)
                        , isColNamed_(true)
{
  // create an empty ReadWriteCsv
  p_data_ = new ReadWriteCsv();

  // for each field Try a String conversion
  for(Integer iVar = df.firstCol(); iVar<=df.lastCol(); iVar++)
  {
    // add an empty string variable
    p_data_->push_back();
    // use virtual method exportString() for getting Strings into Csv
    if (df.elt(iVar)) df.elt(iVar)->exportString(p_data_->back());
  }
}

ExportToCsv::~ExportToCsv()
{ if (p_data_) delete p_data_;}

/* Set a name to each column of the ReadWriteCsv using the form
 *  prefix + number.
 *  @param prefix the prefix to use for the names of the columns
 **/
void ExportToCsv::setColumnsNames(String const& prefix)
{
  // get first and last column indexes
  const Integer first = p_data_->first(), last = p_data_->last();
  // set names of the variables
  for(Integer i = first; i<=last; i++)
  {
    p_data_->setName(i, prefix + typeToString(i)) ;
  }
  isColNamed_ = true;
}

} // namespace STK
