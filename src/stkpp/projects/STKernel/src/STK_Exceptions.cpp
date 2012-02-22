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

/** @file STK_Exceptions.cpp
 *  @brief In this file we implement the standard exceptions of the stk++
 *  library.
 **/

#include "../include/STK_Exceptions.h"

namespace STK
{

// Exception
Exception::Exception() throw() {}
Exception::~Exception() throw() {}
const Char* Exception::error() const throw()
{ return _T("Exception thrown\n");}

// logic_error
logic_error::logic_error( String const& __arg)
                        : _M_msg(__arg)
{}
logic_error::~logic_error() throw() {}
const Char* logic_error::error() const throw()
{ return _M_msg.c_str();}

// domain_error
domain_error::domain_error( String const& __arg)
                          : logic_error(__arg)
{}
domain_error::~domain_error() throw() {}

// invalid_argument
invalid_argument::invalid_argument( String const& __arg)
                                  : logic_error(__arg)
{}
invalid_argument::~invalid_argument() throw() {}

// length_error
length_error::length_error( String const& __arg)
                          : logic_error(__arg)
{}
length_error::~length_error() throw() {}

// out_of_range
out_of_range::out_of_range( String const& __arg)
                          : logic_error(__arg)
{}
out_of_range::~out_of_range() throw() {}

// runtime error
runtime_error::runtime_error( String const& __arg)
                            : _M_msg(__arg)
{}
runtime_error::~runtime_error() throw() {}
const Char* runtime_error::error() const throw()
{ return _M_msg.c_str();}

// range_error
range_error::range_error( String const& __arg)
                        : runtime_error(__arg)
{}
range_error::~range_error() throw() {}

// overflow_error
overflow_error::overflow_error( String const& __arg)
                              : runtime_error(__arg)
{}
overflow_error::~overflow_error() throw() {}

// underflow_error
underflow_error::underflow_error( String const& __arg)
                                : runtime_error(__arg)
{}
underflow_error::~underflow_error() throw() {}

}
