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
 * Purpose:  Declaration of the class ReaWriteCsv.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ReadWriteCsv.h
 *  @brief In this file we define the class ReadWriteCsv.
 *
 *  This class allow to handle csv files in various ways.
 **/

#ifndef READWRITECSV_H
#define READWRITECSV_H

#include "STK_Variable.h"
#include "../../Arrays/include/STK_Array1D.h"

namespace STK
{

/** The STK::Csv namespace contains the constants used in the DManager
 *  project when using csv files.
 **/
namespace Csv
{
  /** @ingroup DManager
   *  Set Array1D<string>::reserve() with this value*/
  static const Integer   DEFAULT_RESERVE    = Integer(0x0FFFF);

  /** @ingroup DManager
   *  Defines the default field delimiter in a file
   **/
  static const Char* DEFAULT_DELIMITER  =     _T(",");

  /** @ingroup DManager
   *  enumeration of read flags
   *  These flags determine the behavior of the reading methods.
   **/
  enum readflags
  {
    RF_APPEND_DATA    = 0x00000004,
    RF_REPLACE_DATA   = 0x00000008
  };
}


/** @ingroup DManager
 *  @brief the ReadWriteCsv class : allow to write and/or to read a csv
 *  file.
 *
 *  It is possible to merge two csv files and to extract subregion of
 *  the file too. All data are stored in a String format.
 **/
class ReadWriteCsv
{
  public:
    /** The default constructor. Instantiates an instance of ReadWriteCsv
     *  and initialize data members to default values.
     *  @param read_names true if we want ot read the names of the
     *  variables at the first line of the file
     **/
    ReadWriteCsv(bool read_names = true);

    /** Misc. constructor.  Instantiates an instance of ReadWriteCsv with the
     *  specified read flags.
     *  @param file_name name of the file to read/write
     *  @param read_names true if we want ot read the names of the
     *  variables at the first line of the file
     *  @param delimiter a String of delimiters
     *  @param reserve the place to reserve for the data
     **/
    ReadWriteCsv( std::string const& file_name
                , bool read_names = true
                , String const& delimiter = Csv::DEFAULT_DELIMITER
                , Integer const&  reserve   = Csv::DEFAULT_RESERVE
                );

    /** Copy constructor. Instantiates an instance of ReadWriteCsv with
     *  the contents of another ReadWriteCsv.
     *  @param df the ReadWriteCsv to copy
     **/
    ReadWriteCsv(ReadWriteCsv const& df);

    /** destructor : Calls clear() before the ReadWriteCsv is destroyed.
     **/
    ~ReadWriteCsv();

    /** Clears the data contained in a ReadWriteCsv
     *  and reclaims any allocated memory.
     **/
    void clear();

    /** Returns the beginning of the container (should be 1).
     *  @return the index of the first variable.
     **/
    inline Integer first() const
    { return str_data_.first(); }

    /** Returns the end of the container (should be NumberOfCols).
     * @return the index of the last variable
     **/
    inline Integer last() const
    { return str_data_.last(); }

    /** Returns the number of variables currently in the ReadWriteCsv.
     * @return The current number of variables of the ReadWriteCsv
     **/
    inline Integer size() const
    { return str_data_.size(); }

    /** Returns the first index of the samples of the variable @c icol.
     *  @param icol index of the variable
     *  @return the first index in the column @c icol
     **/
     inline Integer firstVe( Integer const& icol) const
     { return str_data_.at(icol).first();}

    /** Returns the last index of the samples of the variable @c icol.
     *  @param icol index of the column
     *  @return the last index in the column @c icol
     **/
     inline Integer lastVe( Integer const& icol) const
     { return str_data_.at(icol).last();}

   /** Returns the number of samples currently in the Variable icol.
    *  @param icol index of the column
    *  @return the number of element in the column @c icol
    * */
    inline Integer sizeVe( Integer const& icol) const
    { return str_data_.at(icol).size();}

    /** Returns the lower number of first index of samples. */
    Integer firstVe() const;

    /** Returns the largest number of end index of samples. */
    Integer lastVe() const;

    /** Returns the last error encountered by the class.
     * @return the last error encountered
     **/
    inline String const& error() const
    { return msg_error_; }

    /** return the element with the index icol.
     *  @param icol index of the col
     **/
    inline Variable<String>& at(Integer icol)
    { return str_data_.at(icol); }

    /** return the element with the index icol (const).
     *  @param icol index of the col
     **/
    inline const Variable<String>& at(Integer icol) const
    { return str_data_.at(icol); }

    /** Returns a reference to the value at the specified location.
     *  @param icol index of the col
     **/
    inline Variable<String>& operator[](Integer const& icol)
    { return str_data_.at(icol); }

    /** Returns a const reference to the value at the specified location.
     *  @param icol index of the col
     **/
    inline const Variable<String>& operator[](Integer const& icol) const
    { return str_data_.at(icol); }

   /** return the first element. */
    inline Variable<String>& front()
    { return str_data_.at(first());}

    /** return the first element (const). */
    inline const Variable<String>& front() const
    { return str_data_.at(first());}

    /** return the last element. */
    inline Variable<String>& back()
    { return str_data_.at(last());}

    /** return the last element (const). */
    inline const Variable<String>& back() const
    { return str_data_.at(last());}

    /** Reads the default file with the specified read flags.
     *  @return  @c true if successful, @c false if an error is encountered.
     **/
    bool read();

    /** Attempts to add a column with the values contained in data.
     *  @param data the column to push back
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool push_back(const Variable<String>& data = Variable<String>());

    /** Attempts to add a column with the values contained in data.
     *  @param data the column to push front
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool push_front(const Variable<String>& data = Variable<String>());

    /** Looks up the index of a variable given its name, starting at the
     *  specified index.
     *  @param variable_name Name of the variable to search
     *  @param iStartingIndex first index
     *  @return the index if successful, -1 if an error is encountered.
     **/
    Integer colIndex( const String  &variable_name
                    , Integer const& iStartingIndex = 0) const;

    /** Looks up the index of a variable given its name and original
     *  source file, starting at the specified index.
     *  @param variable_name Name of the variable to search
     *  @param sourceFilename Name of the file to search
     *  @param iStartingIndex first index
     *  @return the index if successful, -1 if an error is encountered.
     **/
    Integer colIndex( String const& variable_name
                    , std::string const& sourceFilename
                    , Integer const& iStartingIndex = 0) const;

    /** Assigns the variable name at the specified index to rStr.
     *  @param icol index of the variable
     *  @param rStr result name of the variable
     *  @return the new length of rStr if successful, -1 if an error is
     *  encountered.
     **/
    Integer name( Integer const& icol, String& rStr) const;

    /** Set the variable name @c name at the specified index.
     *  @param icol index of the variable
     *  @param name name of the variable to set
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool setName( Integer const& icol, String const& name);

    /** Assigns lpStr with the data at the target variable.
     *  @param icol index of the col
     *  @param irow index of the row
     *  @param lpStr returned value
     *  @return the new length of lpStr. Returns -1 if an error is encountered.
     **/
    Integer data( Integer const& icol, Integer const& irow, String &lpStr) const;

    /** Assigns lpStr with the data at the target variable.
     *  @param variable_name name of the variable
     *  @param irow index of the row
     *  @param lpStr returned value
     *  @return the new length of lpStr. Returns -1 if an error is encountered.
     **/
    Integer data( String const& variable_name, Integer const& irow, String &lpStr) const;

    /** Assigns rVector with the data at the target variable.
     *  @param icol index of the col
     *  @param rVector vector of the returned data
     *  @return the new length of lpStr. Returns -1 if an error is encountered.
     **/
    Integer data( Integer const& icol, Variable<String>& rVector) const;

    /** Assigns rVector with the data at the target variable.
     *  @param variable_name name of the variable
     *  @param rVector vector of the returned data
     *  @return the new length of lpStr. Returns -1 if an error is encountered.
     **/
    Integer data( String const& variable_name, Variable<String>& rVector) const;

    /** Sets the delimiters to use for parsing data
     *  (delimiters_ is mutable).
     *  @param delimiters delimiters to use
     **/
    inline void setDelimiters( String const& delimiters) const
    { delimiter_ = delimiters; }

    /** Sets the with_names_ value for reading/writting variables names
     *  (with_names_ is mutable).
     *  @param with_names true if we want to read the names of the
     *  variables
     **/
    inline void setWithNames( bool with_names = true) const
    { with_names_ = with_names; }

    /** Sets the reserve value for data storage (reserve_ is mutable).
     *  @param reserve number of place to reserve
     **/
    inline void setReserve(Integer const& reserve) const
    { reserve_ = reserve; }

    /** Attempts to set the specified value to the element (icol, irow).
     *  @param icol index of the col
     *  @param irow index of the row
     *  @param value the value to set
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool setData( Integer const& icol, Integer const& irow, String const& value);

    /** Attempts to append a data to the variable specified by icol.
     *  @param icol index of the col
     *  @param value value to set
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool appendData( Integer const& icol, String const& value);

    /** Attempts to append values from data to the variable specified by icol.
     *  @param icol index of the col
     *  @param data values to set
     *  @return @c true if successful, @c false if an error is encountered.
     **/
    bool appendData( Integer const& icol, const Variable<String>& data);

    /** Deletes the variable whose index is icol from a ReadWriteCsv.
     *  @param icol index of the column to erase
     *  @return  @c true if successful, @c false if an error is encountered.
     **/
    bool eraseColumn(Integer const& icol);

    /** Reads the specified file with the specified read flags.
     *  @param file_name name of the file to read
     *  @return  @c true if successful, @c false if an error is encountered.
     **/
    bool read(std::string const& file_name);

    /** Assigns a ReadWriteCsv equal to another ReadWriteCsv.
     *  @param df the ReadWriteCsv to copy
     **/
    ReadWriteCsv& operator=( ReadWriteCsv const& df);

    /** Appends a ReadWriteCsv to another ReadWriteCsv.
     *  @param df the ReadWriteCsv to append
     **/
    ReadWriteCsv& operator+=( ReadWriteCsv const& df);

    /** Combines ReadWriteCsv(s)
     *  @param df the ReadWriteCsv to add
     **/
    ReadWriteCsv operator+( ReadWriteCsv const& df) const;

    /** Attempts to write the ReadWriteCsv to the location specified by
     *  file_name using the delimiters specified by delimiter_.
     *  @param file_name name of the file to write
     *  @return  @c true if successful, @c false if an error is encountered.
     **/
    bool write( std::string const& file_name) const;

    /** Write the csv to an output stream.
     *  @param os the output stream
     **/
    void write( ostream& os) const;

    /** Write to output stream a selection based on the coordinates
     *  passed (Think of it as highlighting cells in Excel).
     *  @param os the output stream
     *  @param top the top index
     *  @param bottom th bottom index
     *  @param left the left index
     *  @param right the right index
     **/
    void writeSelection( ostream& os
                       , Integer const& top
                       , Integer const& bottom
                       , Integer const& left
                       , Integer const& right
                       ) const;

    /** Returns a reference of the value specified by the given
     *  coordinates.
     *  @param icol index of the col
     *  @param irow index of thhe row
     **/
    inline String& operator()( Integer const& icol, Integer const& irow)
    { return str_data_.at(icol).at(irow);}

    /** Returns a const reference of the value specified by the given
     *  coordinates.
     *  @param icol index of the col
     *  @param irow index of the row
     **/
    inline String const& operator()( Integer const& icol
                                   , Integer const& irow) const
    { return str_data_.at(icol).at(irow);}

  protected:
    /// Name of the Current file read
    mutable std::string  file_name_;
    /// Read and Write names of the variables
    mutable bool with_names_;
    /// Delimiter(s)
    mutable String  delimiter_;
    /// Size of the buffer
    mutable Integer reserve_;
    /// Contain the last error message
    mutable String msg_error_;

    /** Array for the source file_names. */
    Array1D<std::string> source_file_names_;
    /** Array of array for the data.*/
    Array1D< Variable<String> > str_data_;

    /** Protected member function for internal bookeeping.
     *  @param name name of the variable to find
     *  @param offset offset to apply
     *  @return -1 if the variable is not found
     **/
    Integer lookupVariableIndex( String const& name
                               , Integer const& offset=0) const;

    /** Returns the largest number of samples. */
    Integer largestNumberOfRows() const;

    // friend operators
    friend istream& operator>>( istream& is, ReadWriteCsv& df);
    friend ostream& operator<<( ostream& os, ReadWriteCsv const& df);
};

/** @ingroup DManager
 * Read the data from the stream and returns the stream when done.
 *  @param is input stream
 *  @param df the ReadWriteCsv to read
 **/
istream& operator>>( istream& is, ReadWriteCsv& df);

/** @ingroup DManager
 * write the data into the stream and returns the stream when done.
 *  @param os output stream
 *  @param df the ReadWriteCsv to write
 **/
ostream& operator<<( ostream& os, ReadWriteCsv const& df);

} // Namespace STK

#endif // READWRITECSV_H
