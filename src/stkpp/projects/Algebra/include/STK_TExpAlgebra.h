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
 * Project:  Algebra
 * Purpose:  Allow to inline Algebraic expressions using template.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_TExpAlgebra.h
 *  @brief Define the main classes for inlining algebraic expressions
 *  with Reals.
 **/
#ifndef STK_TEXPALGEBRA_H
#define STK_TEXPALGEBRA_H

#include "../../Arrays/include/STK_MatrixSquare.h"
#include "../../Arrays/include/STK_MatrixLowerTriangular.h"
#include "../../Arrays/include/STK_MatrixUpperTriangular.h"

#include "STK_TOperator.h"


namespace STK
{
/** @ingroup Algebra
 *  @brief Binary Operator bbase class.
 * 
 * This class allow to handle the operations
 * Exp Op Exp, Real Op Exp or Exp Op Real in the
 * general case
**/
template<class Op, class ExpLeft, class ExpRight>
class BinOpBase
{
  protected:
    const ExpLeft&  lhs_;        ///< Reference on the Left term
    const ExpRight& rhs_;        ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase(const ExpLeft& X, const ExpRight& Y) : lhs_(X), rhs_(Y)
    { ;}
  
    /** Protected destructor.                                               */
    ~BinOpBase() { ;}

  public:
    /** The operator [i] will perform lhs_[i] Op rhs_[i]
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator[](Integer i) const
    { return Op::apply(lhs_[i], rhs_[i]);}

    /** The operator (i,j) will perform lhs_(i,j) Op rhs_(i,j)
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator()(Integer i, Integer j) const
    { return Op::apply(lhs_(i,j), rhs_(i,j));}
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  Real Op Exp.
**/
template<class Op, class ExpRight>
class BinOpBase<Op, Real, ExpRight>
{
  protected:
    Real const&     lhs_;      ///< Reference on the Left term
    const ExpRight& rhs_;      ///< Reference on the Right Side

    /** Protected constructor                                               */
    BinOpBase(Real const& y, const ExpRight& X) : lhs_(y), rhs_(X)
    { ;}
  
    /**Protected destructor.                                                */
    ~BinOpBase(){ ;}

  public:
    /** The operator[] will perform lhs_ Op rhs_[i]
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator[](Integer i) const
    { return Op::apply(lhs_,rhs_[i]);}

    /** The operator(i,j) will perform lhs_ Op rhs_(i,j)
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator()(Integer i, Integer j) const
    { return Op::apply(lhs_,rhs_(i,j));}
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case Exp Op Real
**/
template<class Op, class ExpLeft>
class BinOpBase<Op, ExpLeft, Real>
{
  protected:
    const ExpLeft& lhs_;      ///< Reference on the Left Side
    Real const&    rhs_;      ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase(const ExpLeft& X, Real const& y) : lhs_(X), rhs_(y)
    { ;}

    /**Protected destructor.                                                */
    ~BinOpBase(){ ;}

 public:
    /** The operator[i] will perform lhs_[i] Op rhs_
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator[](Integer i) const
    { return Op::apply(lhs_[i],rhs_);}

    /** The operator(i,j) will perform lhs_(i) Op rhs_
     *  where Op will be Plus (+), Minus (-), Mult (*) or Div (/)
     **/
    inline Real operator()(Integer i, Integer j) const
    { return Op::apply(lhs_(i,j),rhs_);}
};


/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  Matrix * ExpRight.
**/
template<class ExpRight>
class BinOpBase<Mult, Matrix, ExpRight>
{
  protected:
    Matrix const&   lhs_;      ///< Reference on the Left term
    const ExpRight& rhs_;      ///< Reference on the Right Side

    /** Protected constructor                                               */
    BinOpBase(const Matrix& y, const ExpRight& X) : lhs_(y), rhs_(X)
    { ;}
  
    /**Protected destructor.                                                */
    ~BinOpBase(){ ;}

  public:
    /** The operator[] will perform Matrix(i) * rhs_
     *  when Op is Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = lhs_.firstCol(); j <=lhs_.lastCol(); j++)
        aux += Mult::apply(lhs_(i,j),rhs_[j]);
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  ExpLeft * Matrix
**/
template<class ExpLeft>
class BinOpBase<Mult, ExpLeft, Matrix>
{
  protected:
    const ExpLeft& lhs_;      ///< Reference on the Left Side
    Matrix const&  rhs_;      ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase(const ExpLeft& X, Matrix const& y) : lhs_(X), rhs_(y)
    { ;}

    /** Protected destructor.                                               */
    ~BinOpBase(){ ;}

 public:
    /** The operator[] will perform lhs_ Op Matrix(j)
     *  where Op Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = rhs_.firstRow(); j<= rhs_.lastRow(); j++)
        aux += Mult::apply(lhs_[j],rhs_(j,i));
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  Matrix * ExpRight.
**/
template<class ExpRight>
class BinOpBase<Mult, MatrixSquare, ExpRight>
{
  protected:
    MatrixSquare const&   lhs_;      ///< Reference on the Left term
    const ExpRight& rhs_;      ///< Reference on the Right Side

    /** Protected constructor                                               */
    BinOpBase(MatrixSquare const& y, const ExpRight& X) : lhs_(y), rhs_(X)
    { ;}
  
    /** Protected destructor.                                               */
    ~BinOpBase(){ ;}

  public:
    /** The operator[] will perform MatrixSquare(i) * rhs_
     *  when Op is Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = lhs_.firstCol(); j <=lhs_.lastCol(); j++)
        aux += Mult::apply(lhs_(i,j),rhs_[j]);
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  ExpLeft * MatrixSquare
**/
template<class ExpLeft>
class BinOpBase<Mult, ExpLeft, MatrixSquare>
{
  protected:
    const ExpLeft& lhs_;      ///< Reference on the Left Side
    MatrixSquare const&  rhs_;      ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase(const ExpLeft& X, MatrixSquare const& y) : lhs_(X), rhs_(y)
    { ;}

    /** Protected destructor.                                               */
    ~BinOpBase(){ ;}

 public:
    /** The operator[] will perform lhs_ Op MatrixSquare(j)
     *  where Op Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = rhs_.firstRow(); j<= rhs_.lastRow(); j++)
        aux += Mult::apply(lhs_[j],rhs_(j,i));
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  MatrixUpperTriangular * ExpRight.
**/
template<class ExpRight>
class BinOpBase<Mult, MatrixUpperTriangular, ExpRight>
{
  protected:
    const MatrixUpperTriangular&   lhs_; ///< Reference on the Left term
    const ExpRight& rhs_;      ///< Reference on the Right Side

    /** Protected constructor                                               */
    BinOpBase( const MatrixUpperTriangular& y, const ExpRight& X)
             : lhs_(y), rhs_(X)
    { ;}
  
    /**Protected destructor.                                                */
    ~BinOpBase(){ ;}

  public:
    /** The operator[] will perform MatrixUpperTriangular(i) * rhs_
     *  when Op is Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = lhs_.compFirstHo(i); j <=lhs_.lastCol(); j++)
        aux += Mult::apply(lhs_(i,j),rhs_[j]);
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  ExpLeft * MatrixUpperTriangular
**/
template<class ExpLeft>
class BinOpBase<Mult, ExpLeft, MatrixUpperTriangular>
{
  protected:
    const ExpLeft& lhs_;      ///< Reference on the Left Side
    const MatrixUpperTriangular&  rhs_;  ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase( const ExpLeft& X, const MatrixUpperTriangular& y)
             : lhs_(X), rhs_(y)
    { ;}

    /** Protected destructor.                                               */
    ~BinOpBase(){ ;}

 public:
    /** The operator[] will perform lhs_ Op MatrixUpperTriangular(j)
     *  where Op Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = rhs_.firstRow(); j<= rhs_.compLastVe(i); j++)
        aux += Mult::apply(lhs_[j], rhs_(j,i));
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  MatrixLowerTriangular * ExpRight.
**/
template<class ExpRight>
class BinOpBase<Mult, MatrixLowerTriangular, ExpRight>
{
  protected:
    const MatrixLowerTriangular&   lhs_; ///< Reference on the Left term
    const ExpRight& rhs_;      ///< Reference on the Right Side

    /** Protected constructor                                               */
    BinOpBase( const MatrixLowerTriangular& y, const ExpRight& X)
             : lhs_(y), rhs_(X)
    { ;}
  
    /**Protected destructor.                                                */
    ~BinOpBase(){ ;}

  public:
    /** The operator[] will perform MatrixLowerTriangular(i) * rhs_
     *  when Op is Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = lhs_.firstCol(); j <=lhs_.compLastHo(i); j++)
        aux += Mult::apply(lhs_(i,j),rhs_[j]);
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Specialized Binary Operator Base class for the case
 *  ExpLeft * MatrixLowerTriangular
**/
template<class ExpLeft>
class BinOpBase<Mult, ExpLeft, MatrixLowerTriangular>
{
  protected:
    const ExpLeft& lhs_;      ///< Reference on the Left Side
    const MatrixLowerTriangular&  rhs_;  ///< Reference on the Right term

    /** Protected constructor                                               */
    BinOpBase( const ExpLeft& X, const MatrixLowerTriangular& y)
             : lhs_(X), rhs_(y)
    { ;}

    /** Protected destructor.                                               */
    ~BinOpBase(){ ;}

 public:
    /** The operator[] will perform lhs_ Op MatrixLowerTriangular(j)
     *  where Op Mult (*)
     **/
    inline Real operator[](Integer i) const
    {
      Real aux = 0.0;
      for (Integer j = rhs_.compFirstVe(i); j<= rhs_.lastRow(); j++)
        aux += Mult::apply(lhs_[j], rhs_(j,i));
      return aux;
    }
};

/** @ingroup Algebra
 *  @brief Binary Operator class derived from BinOpBase class which
 *  is specialized.
**/
template<class Op, class ExpLeft, class ExpRight>
class BinOp : public BinOpBase<Op, ExpLeft, ExpRight>
{
  typedef BinOp<Op, ExpLeft, ExpRight> _This;
 
  public:
    /* constructor                                                          */
    BinOp( const ExpLeft& X, const ExpRight& Y)
         : BinOpBase<Op, ExpLeft, ExpRight>(X,Y)
    { ;}

    /** destructor.                                                         */
    ~BinOp(){ ;}
};  

/** @ingroup Algebra
 *  @brief UnOp class for unary operators.
 **/
template<class Op, class Exp>
class UnOp
{
  typedef UnOp<Op, Exp> _This;

  protected:
    const Exp& rhs_;    ///< Right hand Side expression

  public:
    /** Default constructor                                                 */
    UnOp(const Exp& X) : rhs_(X)
    { ;}

    /** destructor.                                                         */
    ~UnOp(){}

    /** Op will be Uplus or Uminus.                                   */
    Real operator[](Integer i) const
    { return Op::apply(rhs_[i]);}

    /** Op will be Uplus or Uminus.                                   */
    Real operator()(Integer i) const
    { return Op::apply(rhs_(i));}

   /** Operator Plus.                                                 */
    template<class Right>
    BinOp<Plus, _This, Right> operator+(const Right& rhs)
    { return BinOp<Plus, _This, Right>(*this, rhs); }

    /** Operator Minus.                                               */
    template<class Right>
    BinOp<Minus, _This, Right> operator-(const Right& rhs)
    { return BinOp<Minus, _This, Right>(*this, rhs);}

    /** Operator Mult.                                                */
    template<class Right>
    BinOp<Mult, _This, Right> operator*(const Right& rhs)
    { return BinOp<Mult, _This, Right>(*this, rhs);}

    /** Operator Div.                                                 */
    template<class Right>
    BinOp<Div, _This, Right>
    operator/(const Right& rhs)
    { return BinOp<Div, _This, Right>(*this, rhs);}
};

} // Namespace STK

#endif //STK_TEXPALGEBRA_H
