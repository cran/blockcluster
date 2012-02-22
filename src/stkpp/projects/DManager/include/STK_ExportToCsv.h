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
 * Project:  DManager
 * Purpose:  Declaration of the class ImportExportToCsv.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 */

/** @file STK_ExportToCsv.h
 *  @brief In this file we define the class ExportToCsv.
 *
 *  These classes allow to import or export data from/to csv files.
 **/

#ifndef STK_EXPORTCSV_H
#define STK_EXPORTCSV_H

#include "../../STKernel/include/STK_String_Util.h"

#include "STK_ReadWriteCsv.h"
#include "STK_DataFrame.h"

#include "../../Sdk/include/STK_ITContainer2D.h"

namespace STK
{

namespace Csv
{
  /** @ingroup DManager
   *  Defines the default prefix to used when naming the columns
   *  of the ReadWriteCsv.
   **/
  static const Char* DEFAULT_COLUMN_PREFIX  =  _T("Var");
}


/** @ingroup DManager
 *  @brief Export data to a Csv stream.
 *
 * A ExportToCsv object create a @c ReadWriteCsv from a container of data
 * like a DataFrame, a Vector or a Matrix. The data are stored in a String
 * format in the @c ReadWriteCsv struct.
 **/
class ExportToCsv
{
  public:
    /** Default constructor. Create an instance of ExportToCvs.
     **/
    ExportToCsv();
    /** Constructor : create an instance of ExportToCvs with a DataFrame.
     *  @param df the DataFrame to export
     **/
    ExportToCsv(const DataFrame& df);

    /** Instantiates an instance of ExportToCvs with a ITContainer1D.
     *  @param A the ITContainer1D to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE, class TContainer1D>
    ExportToCsv( ITContainer1D<TYPE, TContainer1D> const& A, String const& prefix=STRING_NA
               )
               : p_data_(0)
               , isColNamed_(false)
    {
      // create an empty ReadWriteCsv with no variable name
      p_data_ = new ReadWriteCsv(false);

      // for each field Try a String conversion
      const Integer first = A.first(), last = A.last();

      // add an empty string variable
      p_data_->push_back();
      // add strings to the String variable
      for(Integer i = first; i<=last; i++)
      {
        p_data_->back().push_back(typeToString<TYPE>(A.at(i)));
      }
      if (prefix != STRING_NA)
        p_data_->back().setName(prefix) ;
    }

    /** Instantiates an instance of ExportToCvs with a ITArray2D.
     *  @param A the IArray2d to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE, class TContainer2D>
    ExportToCsv( ITContainer2D<TYPE, TContainer2D> const& A, String const& prefix=STRING_NA)
               : p_data_(0)
               , isColNamed_(false)
    {
      // create an empty ReadWriteCsv with no name
      p_data_ = new ReadWriteCsv(false);

      // for each field try a String conversion
      const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
      const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
      for(Integer iVar = firstCol; iVar<=lastCol; iVar++)
      {
        // add an empty string variable (an empty column)
        p_data_->push_back();
        for (Integer irow=firstRow; irow<=lastRow; irow++)
        {
          p_data_->back().push_back(typeToString<TYPE>(A.at(irow,iVar)));
        }
        if (prefix != STRING_NA)
          p_data_->back().setName(prefix+typeToString<Integer>(iVar)) ;
      }
    }

    /** Append a 2D container.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE, class TContainer2D>
    void append( TContainer2D const& A, String const& prefix=STRING_NA)
    {
      // for each field try a String conversion
      const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
      const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
      for(Integer iVar = firstCol; iVar<=lastCol; iVar++)
      {
        // add an empty string variable (an empty column)
        p_data_->push_back();
        for (Integer irow=firstRow; irow<=lastRow; irow++)
        {
          p_data_->back().push_back(typeToString<TYPE>(A.at(irow,iVar)));
        }
        if (prefix != STRING_NA)
          p_data_->back().setName(prefix+typeToString<Integer>(iVar)) ;
      }
    }

    /** Append a 1D container.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE, class TContainer1D>
    void append1D( TContainer1D const& A, String const& prefix=STRING_NA)
    {

        // for each field Try a String conversion
        const Integer first = A.first(), last = A.last();

        // add an empty string variable
        p_data_->push_back();
        // add strings to the String variable
        for(Integer i = first; i<=last; i++)
        {
          p_data_->back().push_back(typeToString<TYPE>(A.at(i)));
        }
        if (prefix != STRING_NA)
          p_data_->back().setName(prefix) ;
    }

    /** Append a 1D container.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE>
    void appendData( TYPE const& A, String const& prefix=STRING_NA)
    {
        // add an empty string variable
        p_data_->push_back();
        // add strings to the String variable
        p_data_->back().push_back(typeToString<TYPE>(A));
        if (prefix != STRING_NA)
          p_data_->back().setName(prefix) ;
    }

    /** Set a name to each column of the ReadWriteCsv using the form
     *  prefix + number.
     *  @param prefix the prefix to use for the names of the columns
     **/
    void setColumnsNames(String const& prefix = Csv::DEFAULT_COLUMN_PREFIX);

    /** destructor.
     *  The protected field p_data_ will be liberated.
     **/
    virtual ~ExportToCsv();

    /** Accesor. Return a ptr on the the ReadWriteCsv. */
    inline ReadWriteCsv* const p_readWriteCsv() const
    { return p_data_;}

    /** release the ReadWriteCsv. It will be freed by the user. */
    inline void release() { p_data_ =0;}

  protected:
    /** ptr on the ReadWriteCsv containing the data. */
    ReadWriteCsv* p_data_;
    /** @c true if the columns have a name. @c false otherwise*/
    bool isColNamed_;
};

} // namespace STK

#endif /*STK_EXPORTCSV_H*/
