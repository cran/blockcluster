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
 * Project:  stkpp::Dmanager
 * created on: 18 oct. 2010
 * Purpose:  .
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Option.cpp
 *  @brief In this file we implement the Option class.
 **/

#include "../include/STK_Option.h"
#include "../include/STK_IPage.h"

#include "../../STKernel/include/STK_String_Util.h"
#include "../../STKernel/include/STK_Exceptions.h"

namespace STK
{

/* defaut constructor
 * @param name the name of the Option
 **/
Option::Option( String const& name, TypeOption type, bool isOptional)
              : name_(name)
              , sep_(CHAR_SEP)
              , type_(type)
              , isOptional_(isOptional)
              , isValued_(false)
              , p_String_(0)
{
  toUpperString(name_);
  setDefaultValue();
}

/* Special constructor. This will construct an option as a sub-page. The
 * name of the option will be the name of the page.
 * @param page the page to set as option
 **/
Option::Option( IPage const& page)
              : name_(page.name())
              , sep_(CHAR_SEP)
              , type_(page_)
              , isOptional_(page.isOptional())
              , isValued_(true)
              , p_Page_(page.clone())
{
  toUpperString(name_);
}

/* Copy constructor.
 * @param opt the Option to copy
 **/
Option::Option( Option const& opt)
              : name_(opt.name_)
              , sep_(opt.sep_)
              , type_(opt.type_)
              , isOptional_(opt.isOptional_)
              , isValued_(opt.isOptional_)
              , p_String_(0)
{
  // copy rhs value
  switch (type_)
  {
    case string_:
      if (opt.p_String_) set(*(opt.p_String_));
      else p_String_ = 0;
      break;
    case real_:
      if (opt.p_Real_) set(*(opt.p_Real_));
      else p_Real_ = 0;
      break;
    case integer_:
      if (opt.p_Integer_) set(*(opt.p_Integer_));
      else p_Integer_ = 0;
      break;
    case range_:
      if (opt.p_Range_) set(*(opt.p_Range_));
      else p_Range_ = 0;
      break;
    case lstring_:
      if (opt.p_lString_) set(*(opt.p_lString_));
      else p_lString_ = 0;
      break;
    case lreal_:
      if (opt.p_lReal_) set(*(opt.p_lReal_));
      else p_lReal_ = 0;
      break;
    case linteger_:
      if (opt.p_lInteger_) set(*(opt.p_lInteger_));
      else p_lInteger_ = 0;
      break;
    case lrange_:
      if (opt.p_lRange_) set(*(opt.p_lRange_));
      else p_lRange_ = 0;
      break;
    case page_:
      if (opt.p_Page_) setPage(*(opt.p_Page_));
      else p_Page_ = 0;
      break;
    case unknown_:
      runtime_error("Error in Option(opt). Unknown type option.\n");
  };
}

/*
 * destructor.
 */
Option::~Option()
{ deleteValue();}

/*
 * copy
 */
Option& Option::operator=( Option const& opt)
{
  // Do the assignment operation of the members
  name_ = opt.name_;
  sep_ = opt.sep_;
  type_ = opt.type_;
  isOptional_ = opt.isOptional_;
  isValued_ = opt.isValued_;
  // copy rhs value
  switch (type_)
  {
    case string_:
      if (opt.p_String_) set(*(opt.p_String_));
      else p_String_ = 0;
      break;
    case real_:
      if (opt.p_Real_) set(*(opt.p_Real_));
      else p_Real_ = 0;
      break;
    case integer_:
      if (opt.p_Integer_) set(*(opt.p_Integer_));
      else p_Integer_ = 0;
      break;
    case range_:
      if (opt.p_Range_) set(*(opt.p_Range_));
      else p_Range_ = 0;
      break;
    case lstring_:
      if (opt.p_lString_) set(*(opt.p_lString_));
      else p_lString_ = 0;
      break;
    case lreal_:
      if (opt.p_lReal_) set(*(opt.p_lReal_));
      else p_lReal_ = 0;
      break;
    case linteger_:
      if (opt.p_lInteger_) set(*(opt.p_lInteger_));
      else p_lInteger_ = 0;
      break;
    case lrange_:
      if (opt.p_lRange_) set(*(opt.p_lRange_));
      else p_lRange_ = 0;
      break;
    case page_:
      if (opt.p_Page_) setPage(*(opt.p_Page_));
      else p_Page_ = 0;
      break;
    case unknown_:
      runtime_error("Error in Option::operator=(opt). Unknown type option.\n");
  };
  // return this
  return *this;  // Return a reference to myself.
}
/* Convert a string in a value
 * @param str the string to convert
 * @return @c true if the conversion success, @c false otherwise
 */
bool Option::setValue( String const& str )
{
  // list of values
  Range rangeValue;
  std::list<String> lStringValue;
  std::list<Real> lRealValue;
  std::list<Integer> lIntegerValue;
  std::list<Range> lRangeValue;
  // choose type
  switch (type_)
  {
    case string_:
      set(str);
      break;
    case real_:
      Real realValue;
      if (stringToType(realValue, str))
        set(realValue);
      else
        runtime_error("Error in Option::setValue(value). Conversion failed.\n");
      break;
    case integer_:
      Integer integerValue;
      if (stringToType(integerValue, str))
        set(integerValue);
      else
        runtime_error("Error in Option::setValue(value). Conversion failed.\n");
      break;
    case range_:
      if (stringToType(rangeValue, str))
      { set(rangeValue);}
      else
        runtime_error("Error in Option::setValue(value). Conversion failed.\n");
      break;
    case lstring_:
      DManager::readList(str, lStringValue);
      set(lStringValue);
      break;
    case lreal_:
      DManager::readList(str, lRealValue);
      set(lRealValue);
      break;
    case linteger_:
      DManager::readList(str, lIntegerValue);
      set(lIntegerValue);
      break;
    case lrange_:
      DManager::readList(str, lRangeValue);
      set(lRangeValue);
      break;
    case page_:
      runtime_error("Error in Option::setValue(value). Page option.\n");
      break;
    case unknown_:
      runtime_error("Error in Option::setValue(value). Unknown type option.\n");
  };
  // error if an error occur in readList ?
  return true;
}

/*  set a value from a Page.
 *  @param value the Page to set
 **/
void Option::setPage( IPage const& value )
{
  deleteValue();
  p_Page_ = value.clone();
  isValued_ = true;
  type_ = page_;
}


/* @brief write out the option in the output stream
 *  @param os output stream
 */
void Option::write( ostream& os) const
{
  // write option name if it's not a page
  if (type_ != page_)
  {  // write name and " = "
    os << name_ << STRING_BLANK << CHAR_EQUAL << STRING_BLANK;
  }
  // write option value
  switch (type_)
  {
    case string_:
      if (p_String_)
        os << *p_String_;
      break;
    case real_:
      if (p_Real_)
        os << *p_Real_;
      break;
    case integer_:
      if (p_Integer_)
        os << *p_Integer_;
      break;
    case range_:
      if (p_Range_)
        os << *p_Range_;
      break;
    case lstring_:
      if (p_lString_)
        DManager::writeList(os, *p_lString_, sep_);
      break;
    case lreal_:
      if (p_lReal_)
      DManager::writeList(os, *p_lReal_, sep_);
      break;
    case linteger_:
      if (p_lInteger_)
        DManager::writeList(os, *p_lInteger_, sep_);
      break;
    case lrange_:
      if (p_lRange_)
        DManager::writeList(os, *p_lRange_, sep_);
      break;
    case page_:
      if (p_Page_)
        p_Page_->write(os);
      break;
    case unknown_:
      runtime_error("Error in Option::setValue(value). Unknown type option.\n");

  };
}

/* @brief read out the option from the input stream
 *  @param is input stream
 */
void Option::read( istream& is)
{
  // read option value
  if (type_ == page_)
  {
      if (p_Page_)
      {
        p_Page_->read(is);
        if (!p_Page_->validate())
        {
          throw runtime_error("Error in Option::read(is) <A sub-page is not validated>.\n");
        }
      }
  };
}

/*  set a value from string.
 *  @param value the string value to set
 **/
void Option::set( String const& value )
{
  deleteValue();
  p_String_ = new String(value);
  isValued_ = true;
}

/*  set a value from a Real.
 *  @param value the real value to set
 **/
void Option::set( Real const& value )
{
  deleteValue();
  p_Real_ = new Real(value);
  isValued_ = true;
}

/*  set a value from an Integer.
 *  @param value the integer value to set
 **/
void Option::set( Integer const& value )
{
  deleteValue();
  p_Integer_ = new Integer(value);
  isValued_ = true;
}

/*  set a value from a Range.
 *  @param value the Range value to set
 **/
void Option::set( Range const& value )
{
  deleteValue();
  p_Range_ = new Range(value);
  isValued_ = true;
}

/*  set a value from a list of Real.
 *  @param value the list real values to set
 **/
void Option::set( std::list<String> const& value )
{
  deleteValue();
  p_lString_ = new std::list<String>(value);
  isValued_ = true;
}

/*  set a value from a list of Real.
 *  @param value the list real values to set
 **/
void Option::set( std::list<Real> const& value )
{
  deleteValue();
  p_lReal_ = new std::list<Real>(value);
  isValued_ = true;
}

/*  set a value from a list of Integer.
 *  @param value the list of integer value to set
 **/
void Option::set( std::list<Integer> const& value )
{
  deleteValue();
  p_lInteger_ = new std::list<Integer>(value);
  isValued_ = true;
}
/*  set a value from a list of Range.
 *  @param value the list of Range to set
 **/
void Option::set( std::list<Range> const& value )
{
  deleteValue();
  p_lRange_ = new std::list<Range>(value);
  isValued_ = true;
}

// delete allocated value
void Option::deleteValue()
{
  switch (type_)
  {
    case string_:
      if (p_String_) delete p_String_;
      p_String_ = 0;
      break;
    case real_:
      if (p_Real_) delete p_Real_;
      p_Real_ = 0;
      break;
    case integer_:
      if (p_Integer_) delete p_Integer_;
      p_Integer_ = 0;
      break;
    case range_:
      if (p_Range_) delete p_Range_;
      p_Range_ = 0;
      break;
    case lstring_:
      if (p_lString_) delete p_lString_;
      p_lString_ = 0;
      break;
    case lreal_:
      if (p_lReal_) delete p_lReal_;
      p_lReal_ = 0;
      break;
    case linteger_:
      if (p_lInteger_) delete p_lInteger_;
      p_lInteger_ = 0;
      break;
    case lrange_:
      if (p_lRange_) delete p_lRange_;
      p_lRange_ = 0;
      break;
    case page_:
      if (p_Page_) delete p_Page_;
      p_Page_ = 0;
      break;
    default:
      break;
  };
  isValued_ = false;
}

/*set A default value to the option.
 */
void Option::setDefaultValue()
{
  //deleteValue();
  switch (type_)
  {
    case string_:
      p_String_ = new String;
      *p_String_ = STRING_NA;
      break;
    case real_:
      p_Real_ = new Real;
      *p_Real_ = Arithmetic<Real>::NA();
      break;
    case integer_:
      p_Integer_ = new Integer;
      *p_Integer_ = Arithmetic<Integer>::NA();
      break;
    case range_:
      p_Range_ = new Range;
      break;
    case lstring_:
      p_lString_ = new std::list<String>;
      break;
    case lreal_:
      p_lReal_ = new std::list<Real>;
      break;
    case linteger_:
      p_lInteger_ = new std::list<Integer>;
      break;
    case lrange_:
      p_lRange_ = new std::list<Range>;
      break;
    case page_:
      p_Page_ = 0;
      break;
    default:
      break;
  };
}

} // namespace STK
