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
 * Purpose:  Implementation of the class ImportFromCsv.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 */

/** @file STK_ImportFromCsv.cpp
 *  @brief In this file we implement the classes ImportFromCsv.
 **/

#include "../include/STK_ImportFromCsv.h"
#include "../../STKernel/include/STK_Exceptions.h"

namespace STK
{

// Default constructor for ImportCsv
ImportFromCsv::ImportFromCsv( ReadWriteCsv const& import)
                            : p_dataFrame_(new DataFrame())
                            , import_(import)
{ }


// virtual dtor
ImportFromCsv::~ImportFromCsv()
{ }

/* launch the conversion from the ReadWriteCsv to a Frame. */
bool ImportFromCsv::run( Import::TypeImport type)
{
  switch (type)
  {
    case Import::numeric_:
      return runNumeric();
      break;
    case Import::only_numeric_:
      return runOnlyNumeric();
      break;
    case Import::string_:
      return runString();
      break;
    default:
      return false;
      break;
  };
}

/* launch the conversion from the ReadWriteCsv to a DataFrame. */
bool ImportFromCsv::runNumeric()
{
  try
  {
    // for each field Try a numeric conversion
    for (Integer j =import_.first(); j <=import_.last(); j++)
    {
      Variable<Real>* pvReal = new Variable<Real>();
      // test number of successful conversion
      if (convertToTYPE(j, *pvReal))
      {  // if no failure add variable to the dataframe
         p_dataFrame_->pushBackVariable(pvReal);
      }
      else
      {
        delete pvReal; // delete varReal
        Variable<String>* pvString = new Variable<String>(import_[j]);
        p_dataFrame_->pushBackVariable(pvString);
      }
    }
  }
  catch(const Exception& error)
  {
    msg_error_  = error.error();
    msg_error_ += _T("\nIn ImportCsv::runNumeric()");
    return false;
  }
  return true;
}

/** launch the conversion from the ReadWriteCsv to a Frame. */
bool ImportFromCsv::runOnlyNumeric()
{
  try
  {
    // for each field Try a numeric conversion
    for (Integer j =import_.first(); j <=import_.last(); j++)
    {
      Variable<Real>* pvReal = new Variable<Real>();
      // test number of successful conversion
      if (convertToTYPE(j, *pvReal))
      {  // if no failure add variable to the dataframe
         p_dataFrame_->pushBackVariable(pvReal);
      }
      else delete pvReal;
    }
  }
  catch(const Exception& error)
  {
    msg_error_  = error.error();
    msg_error_ += _T("\nIn ImportCsv::runOnlyNumeric()");
    return false;
  }
  return true;
}

/* launch the conversion from the ReadWriteCsv to a DataFrame. */
bool ImportFromCsv::runString()
{
  try
  {
    // for each field Try a numeric conversion
    for (Integer j =import_.first(); j <=import_.last(); j++)
    {
      Variable<String>* pvString = new Variable<String>(import_[j]);
      p_dataFrame_->pushBackVariable(pvString);
    }
  }
  catch(const Exception& error)
  {
    msg_error_ = error.error();
    msg_error_ += _T("\nIn ImportCsv::runString()");
    return false;
  }
  return true;
}


} // namespace STK
