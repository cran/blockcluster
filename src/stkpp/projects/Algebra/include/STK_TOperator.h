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

    Contact : Serge.Iovleff@stkpp.org                                   */

/*
 * Project:  Algebra
 * Purpose:  Allow to inline binary and unary operators.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_TOperator.h
 *  @brief Define the templated classes for handling binary and
 *  unary operators. This file should not be include by other files.
 **/
 

#ifndef STK_TOPERATOR_H
#define STK_TOPERATOR_H

#include "../../STKernel/include/STK_Real.h"

namespace STK
{
/* forward declaration.*/
template<class Op, class ExpLeft, class ExpRight>
class BinOp;

template<class Op, class Exp>
class UnOp;


/** @ingroup Algebra
 *  @brief Binary Operator Plus
 **/
struct Plus
{
   inline static Real apply(Real const& x, Real const& y)
  { return x+y;}
};

/** @ingroup Algebra
 *  @brief Binary Operator Minus
 **/
struct Minus
{
   inline static Real apply(Real const& x, Real const& y)
  { return x-y;}
};
 
/** @ingroup Algebra
 *  @brief Binary Operator Mult
 **/
struct Mult
{
  inline static Real apply(Real const& x, Real const& y)
  { return x*y;}
};

/** @ingroup Algebra
 *  @brief Binary Operator Div
 **/
struct Div
{
  inline static Real apply(Real const& x, Real const& y)
  { return x/y;}
};

/** @ingroup Algebra
 *  @brief Unary Operator Plus
 **/
struct Uplus
{ 
  inline static Real apply(Real const& y)
  { return +y;}
};

/** @ingroup Algebra
 *  @brief Unary Operator Minus
 **/
struct Uminus
{
  inline static Real apply(Real const& y)
  { return -y;}
};
 
/* These Operator construct the recursion.                              */
/** Adding two Expressions.                                             */
template< class ExpLeft, class ExpRight>
BinOp<Plus, ExpLeft, ExpRight >
     operator+( const ExpLeft& lhs, const ExpRight& rhs)
{ return BinOp<Plus, ExpLeft, ExpRight >(lhs, rhs);}

/** soustracting A ExpRight to An Expression.                               */
template< class ExpLeft, class ExpRight>
BinOp<Minus, ExpLeft, ExpRight >
     operator-( const ExpLeft& lhs, const ExpRight& rhs)
{ return BinOp<Minus, ExpLeft, ExpRight >(lhs, rhs);}

/** Mult A ExpRight to An Expression.                                  */
template< class ExpLeft, class ExpRight>
BinOp<Mult, ExpLeft, ExpRight >
     operator*( const ExpLeft& lhs, const ExpRight& rhs)
{ return BinOp<Mult, ExpLeft, ExpRight >(lhs, rhs);}

/** Div A ExpRight to An Expression.                                   */
template< class ExpLeft, class ExpRight>
BinOp<Div, ExpLeft, ExpRight >
     operator/( const ExpLeft& lhs, const ExpRight& rhs)
{ return BinOp<Div, ExpLeft, ExpRight >(lhs, rhs);}

/** Unary Minus Operator.                                               */
template<class ExpRight>
UnOp<Uminus, ExpRight > operator-(const ExpRight& rhs)
{ return UnOp<Uminus, ExpRight >(rhs);}

/** Unary Plus Operator.                                                */
template<class ExpRight>
UnOp<Uplus, ExpRight > operator+(const ExpRight& rhs)
{ return UnOp<Uplus, ExpRight >(rhs);}

} // namespace STK

#endif /*STK_TOPERATOR_H*/
