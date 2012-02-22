/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * Project:  stkpp::STKernel
 * created on: 6 janv. 2012
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Exceptions.h
 *  @brief In this file we define the standard Exceptions of the stk++ library.
 **/

#ifndef STK_EXCEPTIONS_H
#define STK_EXCEPTIONS_H

#include <exception>
#include "STK_String.h"

namespace STK
{
/** @ingroup Base
 *  @brief Base class for all library Exceptions.
 *
 *  This is the base class for all Exceptions thrown by the standard
 *  library, and by certain language expressions.  You are free to derive
 *  your own %Exception classes, or use a different hierarchy, or to
 *  throw non-class data (e.g., fundamental types).
 */
class Exception : public std::exception
{
  public:
    Exception() throw();
    virtual ~Exception() throw();

    /** Returns a C-style character string describing the general cause
     *  of the current error.  */
    virtual const Char* error() const throw();
};


/**  @ingroup Base
 *  Logic errors represent problems in the internal logic of a program;
 *  in theory, these are preventable, and even detectable before the
 *  program runs (e.g., violations of class invariants).
 *  @brief One of two subclasses of Exception.
 */
class logic_error : public Exception
{
    String _M_msg;

  public:
    /** Takes a character String describing the error.  */
    explicit logic_error(const String& __arg);
    virtual ~logic_error() throw();

    /** Returns a C-style character String describing the general cause of
     *  the current error (the same String passed to the ctor).  */
    virtual const Char* error() const throw();
};

/** @ingroup Base
 *  Thrown by the library, or by you, to report domain errors (domain in
 *  the mathematical sense).  */
class domain_error : public logic_error
{
  public:
    explicit domain_error(const String& __arg);
    virtual ~domain_error() throw();
};

/** @ingroup Base
 *  Thrown to report invalid arguments to functions.  */
class invalid_argument : public logic_error
{
  public:
    explicit invalid_argument(const String& __arg);
    virtual ~invalid_argument() throw();
};

/** @ingroup Base
 *  Thrown when an object is constructed that would exceed its maximum
 *  permitted size (e.g., a basic_String instance).  */
class length_error : public logic_error
{
  public:
    explicit length_error(const String& __arg);
    virtual ~length_error() throw();
};

/** @ingroup Base
 *  This represents an argument whose value is not within the expected
 *  range (e.g., boundary checks in basic_String).  */
class out_of_range : public logic_error
{
  public:
    explicit out_of_range(const String& __arg);
    virtual ~out_of_range() throw();
};

/** @ingroup Base
 *  Runtime errors represent problems outside the scope of a program;
 *  they cannot be easily predicted and can generally only be caught as
 *  the program executes.
 *  @brief One of two subclasses of Exception.
 */
class runtime_error : public Exception
{
    String _M_msg;

  public:
    /** Takes a character String describing the error.  */
    explicit runtime_error(const String& __arg);
    virtual ~runtime_error() throw();

    /** Returns a C-style character String describing the general cause of
     *  the current error (the same String passed to the ctor).  */
    virtual const Char* error() const throw();
};

/** @ingroup Base
 *  Thrown to indicate range errors in internal computations.  */
class range_error : public runtime_error
{
  public:
    explicit range_error(const String& __arg);
    virtual ~range_error() throw();
};

/** @ingroup Base
 *  Thrown to indicate arithmetic overflow.  */
class overflow_error : public runtime_error
{
  public:
    explicit overflow_error(const String& __arg);
    virtual ~overflow_error() throw();
};

/** @ingroup Base
 *  Thrown to indicate arithmetic underflow.  */
class underflow_error : public runtime_error
{
  public:
    explicit underflow_error(const String& __arg);
    virtual ~underflow_error() throw();
};

} // namespace STK

#endif /* STK_EXCEPTIONS_H */
