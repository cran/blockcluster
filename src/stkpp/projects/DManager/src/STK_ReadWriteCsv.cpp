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
 * Purpose:  Implementation of the class ReadWriteCsv.h.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ReadWriteCsv.cpp
 *  @brief In this file we implement the class ReadWriteCsv.
 **/

// required STL headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "../include/STK_ReadWriteCsv.h"
#include "../include/STK_DManager_Util.h"

namespace STK
{

/** @ingroup DManager
  * @brief some errors messages.
 **/
static const std::string ERRORCODES[] =
{
  "In ReadWriteCsv UNKNOWN ERROR",
  "In ReadWriteCsv::variable() An unknown error occurred!",
  "In ReadWriteCsv Variable name not found!",
  "In ReadWriteCsv Filename name not found!",
  "In ReadWriteCsv File not found!",
  "In ReadWriteCsv The Number of Headers is different"
                   " from the Number of Data Columns!"
 };

// using declarations
using namespace Csv;

/**  @ingroup DManager
 *  Compute the maximal Length of the strings in V.
 *  @param V the String to treat.
 **/
static Integer maxLength(const Variable<String>& V)
{
  if (V.empty()) return 0;

  // initialize
  Integer maxlength = -Arithmetic<Integer>::max();
  // loop over the values
  for (Integer i=V.first(); i<=V.last(); i++)
  {
    if (!Arithmetic<String>::isNA(V[i])) // update
      maxlength = max(maxlength, (Integer(V[i].size())));
  }
  return maxlength;
}


/** @ingroup DManager
 *  Counts the number of columns in a line stored in a String
 *  and return the position of the delimiters and its types.
 *  @param line The String to parse
 *  @param delimiters The String with the delimiters
 *  @param typeDelimiter The String with all the delimiters encountered
 *  @return the number of field in the line
 **/
static Integer CountCols( String const& line
                        , String const& delimiters
                        , Variable<String>& typeDelimiter
                        )
{
  // number of fields in the file
  Integer numField = 0;
  // no delimiters
  typeDelimiter.resize();
  // Find first delimiter
  String::size_type idx = line.find_first_of(delimiters);
  // if the position of the delimiter char is found (there is a position)
  while (idx != line.npos)
  {
    // add a column
    ++numField;
    // save delimiter
    typeDelimiter.push_back(line.substr(idx,1));
    // find next delimiter
    idx = line.find_first_of(delimiters, idx+1);
  }
  // last delimiter is the end of line
  typeDelimiter.push_back(_T("\n"));
  // return the number of fields
  return ++numField;
}

// The default constructor.  Instantiates an instance of CDataFile and 
// initialize data members to default values.
// Instantiates an instance of CDataFile with the
// specified read flags.
ReadWriteCsv::ReadWriteCsv( bool read_names)
                          : file_name_("")
                          , with_names_(read_names)
                          , delimiter_(DEFAULT_DELIMITER)
                          , reserve_(DEFAULT_RESERVE)
                          , msg_error_("")
{ }
  
// Misc. constructor.  Instantiates an instance of ReadWriteCsv
// with specified read flags.
ReadWriteCsv::ReadWriteCsv( std::string const& file_name
                          , bool read_names
                          , String const& delimiter
                          , Integer const& reserve
                          )
                          : file_name_(file_name)
                          , with_names_(read_names)
                          , delimiter_(delimiter)
                          , reserve_(reserve)
                          , msg_error_("")
{ }

// Copy constructor.  Instantiates an instance of ReadWriteCsv with the
// contents of another ReadWriteCsv.
ReadWriteCsv::ReadWriteCsv(ReadWriteCsv const& df)
{ *this = df;}

// destructor
ReadWriteCsv::~ReadWriteCsv()
{ clear();}

// sets one ReadWriteCsv equal to another ReadWriteCsv.
ReadWriteCsv& ReadWriteCsv::operator=(ReadWriteCsv const& df)
{
  delimiter_         = df.delimiter_;
  reserve_           = df.reserve_;
  with_names_        = df.with_names_;
  file_name_         = df.file_name_;
  msg_error_         = df.msg_error_;
  source_file_names_ = df.source_file_names_;
  str_data_          = df.str_data_;

  return *this;
}

// adds one ReadWriteCsv to another ReadWriteCsv
ReadWriteCsv& ReadWriteCsv::operator+=(ReadWriteCsv const& df)
{
  for ( Integer i=df.source_file_names_.first()
      ; i<=df.source_file_names_.last()
      ; i++)
    source_file_names_.push_back(df.source_file_names_[i]);

  for ( Integer i=df.str_data_.first()
      ; i<=df.str_data_.last()
      ; i++)
    str_data_.push_back(df.str_data_[i]);

  return *this;
}

// adds ReadWriteCsv(s) together
ReadWriteCsv ReadWriteCsv::operator+(ReadWriteCsv const& df) const
{ // copy this, add df and return the result
  return ReadWriteCsv((*this)) += df;
}

// try to set a value at the (icol, irow) place
bool ReadWriteCsv::setData( Integer const& icol
                          , Integer const& irow
                          , String const& value)
{
  try
  {
    str_data_.at(icol).at(irow) = value;
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

bool ReadWriteCsv::appendData(Integer const& icol, String const& value)
{
  try
  {
    if (Arithmetic<String>::isNA(value))
      str_data_[icol].push_back(Arithmetic<String>::NA());
    else
      str_data_[icol].push_back(value);
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

// Deletes the variable whose index is icol from a ReadWriteCsv.
// Returns true if successful, false if an error is encountered.
bool ReadWriteCsv::eraseColumn(Integer const& icol)
{
  try
  {
    // delete the variable from source_file_names_
    source_file_names_.erase(icol),

    // delete the variable from str_data_
    str_data_.erase(icol);

    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

Integer ReadWriteCsv::name( Integer const& icol, String& rStr) const
{
  try
  {
    rStr = str_data_.at(icol).name();
    return static_cast<Integer> (rStr.length());
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...)  { msg_error_ = ERRORCODES[1]; }
  return -1;
}

/* Set the variable name @c name at the specified index.
 *  @param icol index of the variable
 *  @param name name of the variable to set
 *  @return @c true if successful, @c false if an error is encountered.
 **/
bool ReadWriteCsv::setName( Integer const& icol, String const& name)
{
  try
  {
    str_data_.at(icol).setName(name);
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...)  { msg_error_ = ERRORCODES[1]; }
  return false;
}

Integer ReadWriteCsv::largestNumberOfRows() const
{
  Integer retVal = 0;

  for (Integer i=str_data_.first(); i<=str_data_.last(); i++)
    retVal = max(retVal, str_data_[i].size());

  return retVal;
}

/* Returns the largest number of end index of samples.
**/
Integer ReadWriteCsv::lastVe() const
{
  Integer retVal = Arithmetic<Integer> ::min();

  const Integer first = str_data_.first(), last = str_data_.last();
  for (Integer i=first; i<=last; i++)
  {
    //stk_cout << str_data_.at(i);
    retVal = max(retVal, lastVe(i));
  }
  // return the maximal number of row
  return retVal;
}

/* Returns the lower number of end index of samples.
**/
Integer ReadWriteCsv::firstVe() const
{
  Integer retVal = Arithmetic<Integer> ::max();

  const Integer first = str_data_.first(), last = str_data_.last();
  for (Integer i=first; i<=last; i++)
    retVal = min(retVal, firstVe(i));

  return retVal;
}

// Clears all data in the ReadWriteCsv
void ReadWriteCsv::clear()
{
  msg_error_ = "";
  source_file_names_.clear();
  str_data_.clear();
}

// Returns the length of the String if successful. 
// Returns -1 if an error is encountered.
Integer ReadWriteCsv::data( Integer const& icol
                                   , Integer const& irow
                                   , String &lpStr) const
{
  Integer retVal = 0;

  try
  {
    lpStr  = str_data_.at(icol).at(irow);
    retVal = static_cast<Integer> (lpStr.length());
  }
  catch(const Exception& e) 
  {  
    msg_error_ = e.error(); 
    retVal = -1;
  }
  catch(...) // other Exceptions
  { 
    msg_error_ = ERRORCODES[1]; 
    retVal = -1;
  }
  return retVal;
}

// Returns the length of the String if successful. 
// Returns -1 if an error is encountered.
Integer ReadWriteCsv::data( String const& variable_name
                                  , Integer const& irow
                                  , String& rStr) const
{
  Integer retVal = 0;
  Integer iVar   = lookupVariableIndex(variable_name); // find col
  
  if(iVar != -1)
    retVal = data(iVar, irow, rStr);
  else
  {
    msg_error_ = ERRORCODES[5]; 
    retVal = -1;
  }
  return retVal;
}


// Returns the new size of rVector if successful.
// Returns -1 if an error is encountered.
Integer ReadWriteCsv::data( Integer const& icol, Variable<String>& rVector) const
{
  Integer retVal = 0;

  try
  {
    rVector = str_data_.at(icol);
    retVal  = rVector.size();
  }
  catch(const Exception& e) 
  {  
    msg_error_ = e.error();
    retVal = -1;
  }
  catch(...) 
  { 
    msg_error_ = ERRORCODES[2]; 
    retVal = -1;
  }
  return retVal;
}

// Returns the new size of rVector if successful.
// Returns -1 if an error is encountered.
Integer ReadWriteCsv::data( String const& variable_name
                              , Variable<String>& rVector) const
{
  Integer index = lookupVariableIndex(variable_name);
  
  if(index != -1)
    return data(index, rVector);
  else
    msg_error_ = ERRORCODES[5]; 
  return -1;
}

// Returns the index of the first instance of the specified variable
// found AFTER iStartingIndex.
// Returns -1 if the variable is not found.
Integer ReadWriteCsv::colIndex( String const& variable_name
                              , Integer const& iStartingIndex
                              ) const
{
  return lookupVariableIndex(variable_name, iStartingIndex);
}

// Returns the index of the specified variable.
// Returns -1 if the variable is not found.
Integer ReadWriteCsv::colIndex( String const& variable_name
                              , std::string const& sourceFilename
                              , Integer const& iStartingIndex
                              ) const
{
  Integer it = source_file_names_.first();
  while (it <= source_file_names_.last())
  { if (source_file_names_[it] == sourceFilename) break;
    it++;
  }

  if(it == source_file_names_.last()+1) // sourceFilename was not found
  {
    msg_error_ = ERRORCODES[6];
    return -1;
  }

  Integer offset = source_file_names_.first();
  while(it != source_file_names_.first())
  {
    offset++;
    it--;
  }
  
  return lookupVariableIndex(variable_name, iStartingIndex+offset);
}

// Returns the index of the first variable name that matches szName.
// Returns -1 if szName is not found.
Integer ReadWriteCsv::lookupVariableIndex( String const& variable_name
                                     , Integer const& offset /*=0*/) const
{
  Integer it = str_data_.first() + offset;
  while (it <= str_data_.last())
  { if (str_data_[it].name() == variable_name) break;
    it++;
  }

  if(it == str_data_.last()+1) // variable Name was not found
    return -1;

  Integer retVal = str_data_.first();
  while(it != str_data_.first())
  {
    retVal++;
    it--;
  }
  return retVal;
}

// Attempts to add a variable with the name specified by data
// and the values contained in data. Returns true if successful,
// false if an error is encountered.
bool ReadWriteCsv::push_back( const Variable<String>& data)
{
  try
  {
    source_file_names_.push_back(file_name_);
    str_data_.push_back(data);
    str_data_.back().reserve(reserve_);
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

// Attempts to add a variable with the values contained in data.
// Returns true if successful,
// false if an error is encountered.
bool ReadWriteCsv::push_front( const Variable<String>& data)
{
  try
  {
    source_file_names_.push_front(file_name_);
    str_data_.push_front(data);
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

// Reads the default file.
// Returns true if successful, false if an error occurred.
bool ReadWriteCsv::read()
{
  return read(file_name_);
}
// Reads the specified file.
// Returns true if successful, false if an error occurred.
bool ReadWriteCsv::read(std::string const& file_name)
{
  try
  {
    // update file_name
    file_name_ = file_name;
    // input file stream
    ifstream inFile;
    // open file
    inFile.open(file_name.c_str());
    // check error
    if (inFile.rdstate() & std::ios::failbit)
    {
      inFile.close();
      msg_error_ = ERRORCODES[4];
      msg_error_ += "\nFile: " + file_name;
      return false;
    }
    // read file
    inFile >> *this;
    // close file
    inFile.close();
    return true;
  }
  catch(const Exception& e)
  {
    msg_error_ = e.error();
    msg_error_ += "\nIn ReadWriteCsv::read(" + file_name + ")";
  }
  catch(...)
  {
    msg_error_ = ERRORCODES[0];
    msg_error_ += "\nIn ReadWriteCsv::read(" + file_name + ")";
  }

  return false;
}

// write the ReadWriteCsv in a file
bool ReadWriteCsv::write( const std::string &file_name) const
{
  file_name_  = file_name;
  try
  {
    ofstream os(file_name.c_str());
    writeSelection( os
                  , firstVe()
                  , lastVe()
                  , str_data_.first()
                  , str_data_.last()
                  );
    os.close();
    return true;
  }
  catch(const Exception& e) { msg_error_ = e.error(); }
  catch(...) { msg_error_ = ERRORCODES[0]; }
  return false;
}

// write a selection
void ReadWriteCsv::write( ostream& os) const
{
  writeSelection(os, firstVe(), lastVe(), first(), last());
}

// write a selection
void ReadWriteCsv::writeSelection(       ostream& os
                                 , Integer const& top
                                 , Integer const& bottom
                                 , Integer const& left
                                 , Integer const& right) const
{
  // create a vector for the format of the output
  Array1D<Integer>  format(Range(left, right), 0);
  // for each var, find the largest size
  for(Integer iVar=left; iVar<=right; iVar++)
  {
    format.at(iVar)   = maxLength(str_data_.at(iVar));
    if (with_names_)
      format.at(iVar) = max( format.at(iVar)
                           , (Integer )str_data_.at(iVar).name().size());
  }
  // write if needed names variables
  if (with_names_)
    for(Integer iVar=left; iVar<=right; iVar++)
    {
      os << std::setw(format[iVar]) << std::right
         << ConstProxy<String>(str_data_.at(iVar).name())
         << ((iVar==right) ? _T('\n') : delimiter_.at(0));
    }

  // write data
  for(Integer irow = top; irow<=bottom; irow++)
    for(Integer iVar = left; iVar<=right; iVar++)
    {
      try
      {
        os << std::setw(format[iVar]) << std::right
           << ConstProxy<String>(str_data_.at(iVar).at(irow));
      }
      catch(...)
      {
        // if an error occur, we put NA value
        os << std::setw(format[iVar]) << std::right << STRING_NA;
      }
      os << ((iVar==right) ? _T('\n') : delimiter_.at(0));
    }
}

// reads the data from the stream and returns the stream when done.
istream& operator>>(istream& is, ReadWriteCsv& df)
{
  try
  {
    // clear previous ReadCvs if we don't want to append
    if (!(Csv::RF_APPEND_DATA)) df.clear();
    // compute number of existing variables
    Integer colOffset = df.size();
    // initialize the initial number of variables to 0
    Integer nbVars = 0;

    // aux variable for handling delimiters
    Variable<String>    typeDelimiter;
    // set filname    
    df.source_file_names_.push_back(df.file_name_);
    // load file in memory
    stringstream inBuffer;
    inBuffer << is.rdbuf();

    // If the names are at the top line
    if (df.with_names_)
    {
      // get current line
      String lineBuffer;
      // Count the number of names of the first line
      Integer numField;
      do
      {
        // get current line in strBuffer
        std::getline(inBuffer, lineBuffer);
        DManager::removeCharBeforeAndAfter(lineBuffer, CHAR_BLANK);
        (lineBuffer.size() == 0) ?
          numField = 0
        : numField = CountCols( lineBuffer, df.delimiter_, typeDelimiter);
      }
      while ((numField == 0)&&(!inBuffer.eof()));
      // break if we get the end of file
      if (inBuffer.eof()) return is;

      // Declare an input string stream
      istringstream instream;
      // Reset from possible previous errors.
      instream.clear();
      // Use strBuffer as source of input.
      instream.str(lineBuffer);
      // Loop over the columns
      for(Integer icol=1; icol<=numField; icol++)
      {
        // Append a Col
        df.push_back(Variable<String>());
        df.str_data_.at(icol+colOffset).setName( DManager::getField( instream
                                               , typeDelimiter.elt(icol).at(0))
                                               );
      }
      // Update the number of var
      nbVars = STK::max(nbVars, numField);
    }

    // Read data : loop for all rows
    for (Integer irow=1; !inBuffer.eof(); irow++)
    {
      // current line
      String lineBuffer;
      // number of fields in the line
      Integer numField;
      // loop until we encounter a line
      do
      {
        // get current line in strBuffer
        std::getline(inBuffer, lineBuffer);
        DManager::removeCharBeforeAndAfter(lineBuffer, CHAR_BLANK);
      }
      while ((lineBuffer.size() == 0)&&(!inBuffer.eof()));
      // pass empty line
      if (lineBuffer.size() == 0) continue;
      numField = CountCols( lineBuffer, df.delimiter_, typeDelimiter);
      // Declare an input string stream
      istringstream instream;
      // Reset from possible previous errors.
      instream.clear();
      // Use lineBuffer as source of input.
      instream.str(lineBuffer);
      String fieldValue;
      // first loop on the exisiting cols with data
      Integer attCols = min(numField, nbVars);
      for (Integer icol=1; icol<=attCols; icol++)
      {
        // Read field
        fieldValue = DManager::getField( instream, typeDelimiter.elt(icol).at(0));
        // append Data to the columnn
        df.appendData(icol+colOffset, fieldValue);
      }

      // second loop on the existing cols without data
      // will be passed if (numField < nbVars)
      for (Integer icol=attCols+1; icol<=nbVars; icol++)
        df.str_data_.at(icol+colOffset).push_back(Arithmetic<String>::NA());

      // loop on the non-exisiting cols
      // will be passed if (numField > nbVars)
      for (Integer icol=nbVars+1; icol<=numField; icol++)
      {
        // we create each column with NA values
        df.push_back(Variable<String>( irow-1
                                     , Arithmetic<String>::NA()
                                     , IVariable::giveName(icol+colOffset)
                                     )
                    );
        // Read field
        fieldValue = DManager::getField( instream, typeDelimiter.elt(icol).at(0));
        df.appendData(icol+colOffset, fieldValue);
      }
      // Update the number of variables
      nbVars = STK::max(nbVars, numField);
      // break if we get the end of the file
      if (inBuffer.eof()) break; 
    } // irow loop
  }
  catch(const Exception& e) 
  { 
    df.msg_error_ = e.error();
    throw e; 
  }
  catch(...) 
  { 
    df.msg_error_ = ERRORCODES[0];
    throw Exception();
  }

  return is;
}

ostream& operator<<(ostream& os, ReadWriteCsv const& df)
{
  try
  {
    df.writeSelection( os
                     , df.firstVe()
                     , df.lastVe()
                     , df.str_data_.first()
                     , df.str_data_.last()
                     );
  }
  // catch and re-throw any Exceptions
  catch(const Exception& e) { throw e; }
  catch(...) { throw Exception(); }
  return os;
}

} // Namespace STK
