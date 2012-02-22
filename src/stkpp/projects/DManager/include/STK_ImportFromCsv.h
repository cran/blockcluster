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
 * Project:  stkp::DManager
 * Purpose:  Declaration of the class ImportFromCsv.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 */

/** @file STK_ExportToCsv.h
 *  @brief In this file we define the class ImportFromCsv.
 *
 *  These classes allow to import or export data from/to csv files.
 **/

#ifndef STK_IMPORTFROMCSV_H
#define STK_IMPORTFROMCSV_H

#include "STK_Import_Util.h"
#include "STK_ReadWriteCsv.h"
#include "STK_DataFrame.h"

namespace STK
{

/** @ingroup DManager
 *  @brief import data from a Csv File.
 *
 *  A ImportFromCsv object create a DataFram from a given ReadWriteCsv object.
 *  It will try to convert the given ReadWriteCsv to the predefined
 *  type given by the user.
 **/
class ImportFromCsv
{
  public:
    /** Constructor. Instantiates an instance of ImportFromCvs with the
     *  readWriteCsv to import.
     *  @param import the ReadWriteCsv to import
     **/
    ImportFromCsv( ReadWriteCsv const& import);

    /** destructor. */
    virtual ~ImportFromCsv();

    /** launch the conversion from the ReadWriteCsv to a Frame. */
    bool run( Import::TypeImport type = Import::numeric_);

    /** delete the DataFrame. The DataFrame imported is not deleted when
     * the object ImportFromCsv is deleted. The user can call this method
     * if needed in order to release the memory.
     **/
    inline void eraseReadWriteCsv()
    {
      if (p_dataFrame_) delete p_dataFrame_;
      p_dataFrame_ = 0;
    }

    /** delete the DataFrame. */
    inline String const& error() const
    {return msg_error_;}

    /** Return a ptr on the the data read. */
    inline DataFrame const* dataFrame() const
    { return p_dataFrame_;}

  protected:
    /** A ptr on the resulting DataFrame. */
    DataFrame* p_dataFrame_;

    /// Contain the last error message
    mutable String msg_error_;

  private:
    /** a constant reference on the the original ReadWriteCsv. */
    ReadWriteCsv const& import_;

    /** convert a column of the csv in a real vector.
     *  @param iCol the column to convert
     *  @param col the result stored in a vector
     *  @return @c true if the conversion is successful, @c false otherwise
     **/
    template<class TYPE>
    bool convertToTYPE( Integer const& iCol, Variable<TYPE>& col)
    {
      // get dimensions of the variable
      // try to ConvertType strings to TYPE
      Integer nSuccess = col.importString(import_[iCol]);
      return (nSuccess == col.size());
    }

    /** launch the conversion from the ReadWriteCsv to a DataFrame
     * with a numeric conversion.
     **/
    bool runNumeric();
    /** launch the conversion from the ReadWriteCsv to a DataFrame
     *   with only the successful numeric conversion.
     **/
    bool runOnlyNumeric();
    /** launch the conversion from the ReadWriteCsv to a DataFrame
     *   with only the successful numeric conversion.
     **/
    bool runString();
};

} // namespace STK

#endif /*STK_IMPORTFROMCSV_H*/
