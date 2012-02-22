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
 * Project: stkpp::Arrays
 * Purpose:  Define the templated CellBase, CellVe and CellHo classes.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Cell.h
  * @brief This file define the cell classes for the list classes.
  **/

#ifndef STK_CELL_H
#define STK_CELL_H

namespace STK
{
/** @ingroup Arrays
  * @brief Templated Base class for the Cell in a list .
  * 
  * The class TYPE should have surdefined the operator = if
  * necessary.                                                          
 **/
template<class TYPE>
class CellBase
{
  protected:
    TYPE data_;        ///< Data contained by the Cell

  public:
    /** constructor with a pointer on the data.                                */
    CellBase()
    { ;}

    /** constructor with a reference to the data.                              */
    CellBase( TYPE const& data)
            : data_(data)
    { ;}
    
    /** Copy constructor                                                      */
    CellBase( const CellBase<TYPE> &C)
            : data_(C.data_)
    { ;}

    /** virtual destructor.                                             */
    virtual ~CellBase() { ;}
    
    /** Operator = : overwrite the Cell with C.                         */
    CellBase<TYPE>& operator=(const CellBase<TYPE> &C)
    { data_ = C.data_;}

    /** Operator = : write a value on the cell.                         */
    CellBase<TYPE>& operator=(TYPE const& v)
    { data_ = v;}

    /** give a reference to the data for const.                         */
    TYPE const& data() const { return data_;}

    /** give a reference to the data.                                   */
    TYPE& data() { return data_;}

    /** Set the data.                                                  */
    void setData(TYPE const& data) { data_ = data;}
};

/** @ingroup Arrays
  * @brief Templated class for the Vertical Cell of a Vertical List.
 * 
 * The class CellVe is used by the class ListVe.
 **/
template<class TYPE>
class CellVe : virtual public CellBase<TYPE>
{
  using STK::CellBase<TYPE>::data_;

  protected:
    CellVe<TYPE>* up_;        ///< pointer on the upper cell
    CellVe<TYPE>* down_;        ///< pointer on the cell down

  public:
    /** Default constructor                                                   */
    CellVe( CellVe<TYPE>* up =NULL, CellVe<TYPE>* down =NULL);

    /** constructor with a reference to the data.                              */
    CellVe( CellVe<TYPE>* up
          , CellVe<TYPE>* down
          , TYPE const& data);
    
    /** Copy constructor                                                      */
    CellVe(const CellVe<TYPE> &C);

    /** Virtual destructor.                                                   */
    virtual ~CellVe();

    /** Operator = : overwrite the Cell with C.                         */
    CellVe<TYPE>& operator=(const CellVe<TYPE> &C);

    /** operator = : write a value on the cell.                         */
    CellVe<TYPE>& operator=(TYPE const& v);

    /** Give the adress of the cell up.                                 */
    CellVe<TYPE>* getUp()  const;

    /** Give the adress of the cell down.                               */
    CellVe<TYPE>* getDown()  const;

    /** Set the adress of the cell up.                                  */
    void setUp(CellVe<TYPE>* pcell);

    /** Set the adress of the cell down.                                */
    void setDown(CellVe<TYPE>* pcell);
};

/** @ingroup Arrays
  * @brief Templated class for the Horizontal Cell of a Horizontal List.
  * 
  * The class CellHo is used by the class List1D.
 **/
template<class TYPE>
class CellHo : virtual public CellBase<TYPE>
{
  // needed for templated classes
  using STK::CellBase<TYPE>::data_;
  
  protected:
    CellHo<TYPE> *left_;        ///< pointer on the left cell
    CellHo<TYPE> *right_;       ///< pointer on the right cell

  public:
    /** Default constructor                                                   */
    CellHo(CellHo<TYPE> *left =NULL, CellHo<TYPE> *right =NULL)
    { left_ = left; right_ = right;}
    
    /** constructor with a reference to the data.                              */
    CellHo( CellHo<TYPE>* left, CellHo<TYPE>* right, TYPE const& data)
          : CellBase<TYPE>(data)
    { left_ = left; right_ = right; }
    
    /** Copy constructor                                                      */
    CellHo( const CellHo<TYPE> &C)
          : CellBase<TYPE>(C)
          , left_(C.left_)
          , right_(C.right_)
    { ;}

    /** virtual destructor.                                                   */
    virtual ~CellHo() { ;}
    
    /** Operator = : overwrite the Cell with C.                         */
    CellHo<TYPE>& operator=(const CellHo<TYPE> &C)
    { data_ = C.data_; 
      left_ = C.left_;
      right_ = C.right_;
      return *this;
    }
    
    /** operator = : write a value on the cell.                         */
    CellHo<TYPE>& operator=(TYPE const& v)
    { data_ = v;
      return *this;
    }

    /** Give the adress of the cell left.                               */
    CellHo<TYPE>* getLeft()  const
    { return left_;}

    /** Give the adress of the cell right.                              */
    CellHo<TYPE>* getRight()  const
    { return right_;}

    /** Set the left cell adress.                                       */
    void setLeft(CellHo<TYPE>* left)
    { left_ = left;}
    
    /** Set the right cell adress.                                      */
    void setRight(CellHo<TYPE>* right)
    { right_ = right;}
};

/** @ingroup Arrays
  * @brief Templated class for the 2 Dimensional Cells.
  * The class Cell2D is not used.
 **/
template<class TYPE>
class Cell2D : virtual public CellVe<TYPE>
             , virtual public CellHo<TYPE>
{
  using STK::CellBase<TYPE>::data_;
  using STK::CellHo<TYPE>::up_;
  using STK::CellHo<TYPE>::down_;
  using STK::CellVe<TYPE>::left_;
  using STK::CellVe<TYPE>::right_;

  public:
    /** Default constructor with a pointer on the data.                        */
    Cell2D( const Cell2D<TYPE>* left  =NULL
          , const Cell2D<TYPE>* up    =NULL
          , const Cell2D<TYPE>* right =NULL
          , const Cell2D<TYPE>* down  =NULL);

    /** constructor with a reference to the data.                              */
    Cell2D( const Cell2D<TYPE>* left
          , const Cell2D<TYPE>* up
          , const Cell2D<TYPE>* right
          , const Cell2D<TYPE>* down
          , TYPE const& data);
    
    /** Copy constructor                                                      */
    Cell2D(const Cell2D<TYPE> &C);

    /** virtual destructor.                                             */
    virtual ~Cell2D();

    /** operator = : overwrite the Cell with C.                         */
    Cell2D<TYPE>& operator=(const Cell2D<TYPE> &C);

    /** Operator = : write a value on the cell.                         */
    Cell2D<TYPE>& operator=(TYPE const& v);
};

// Constructor with a pointer on the data
template<class TYPE>
CellVe<TYPE>::CellVe( CellVe<TYPE>* up
                        , CellVe<TYPE>* down)
                        : CellBase<TYPE>()
{ up_ = up; down_ = down;}

// Constructor with a reference to the data
template<class TYPE>
CellVe<TYPE>::CellVe( CellVe<TYPE>* up
                        , CellVe<TYPE>* down
                        , TYPE const& data)
                        : CellBase<TYPE>(data)
{ up_ = up; down_ = down; }

// Copy Constructor.
template<class TYPE>
CellVe<TYPE>::CellVe(const CellVe<TYPE> &C)
{ data_ = C.data_;  up_   = C.up_; down_ = C.down_; }

/* virtual destructor.                                                  */
template<class TYPE>
CellVe<TYPE>::~CellVe() { ;}

/* operator = : overwrite the Cell with C.                              */
template<class TYPE>
CellVe<TYPE>&
CellVe<TYPE>::operator=(const CellVe<TYPE> &C)
{ data_ = C.data_;  up_ = C.up_; down_ = C.down_;
  return *this;
}

/* operator = : overwrite the data with v.                              */
template<class TYPE>
CellVe<TYPE>& CellVe<TYPE>::operator=(TYPE const& v)
{ data_ = v;
  return *this;
}

/* Accessors.                                                          */
template<class TYPE>
CellVe<TYPE>* CellVe<TYPE>::getUp()  const
{ return up_;}

// Give the adress of the cell down
template<class TYPE>
CellVe<TYPE>* CellVe<TYPE>::getDown()  const
{ return down_;}

// Set the adress of the cell up
template<class TYPE>
void CellVe<TYPE>::setUp(CellVe<TYPE>* up)
{ up_ = up;}

// Set the adress of the cell down
template<class TYPE>
void CellVe<TYPE>::setDown(CellVe<TYPE>* down)
{ down_ = down;}

// Constructor with a pointer on the data
template<class TYPE>
Cell2D<TYPE>::Cell2D( const Cell2D<TYPE>* left
                    , const Cell2D<TYPE>* up
                    , const Cell2D<TYPE>* right
                    , const Cell2D<TYPE>* down)
                    : CellBase<TYPE>()
                    , CellHo<TYPE>(left, right)
                    , CellVe<TYPE>(up, down) 
{ ;}

// Constructor with a reference to the data
template<class TYPE>
Cell2D<TYPE>::Cell2D( const Cell2D<TYPE>* left
                    , const Cell2D<TYPE>* up
                    , const Cell2D<TYPE>* right
                    , const Cell2D<TYPE>* down
                    , TYPE const& data)
                    : CellBase<TYPE>(data)
                    , CellHo<TYPE>(left, right)
                    , CellVe<TYPE>(up, down) 
{ ;}

// Copy Constructor.
template<class TYPE>
Cell2D<TYPE>::Cell2D(const Cell2D<TYPE>& C)
{ data_  = C.data_;
  up_    = C.up_;
  down_  = C.down_;
  left_  = C.left_;
  right_ = C.right_;
}

/* virtual destructor.                                                  */
template<class TYPE>
Cell2D<TYPE>::~Cell2D() { ;}

/* operator = : overwrite the Cell with C.                              */
template<class TYPE>
Cell2D<TYPE>&
Cell2D<TYPE>::operator=(const Cell2D<TYPE> &C)
{ data_  = C.data_;
  up_    = C.up_;
  down_  = C.down_;
  left_  = C.left_;
  right_ = C.right_;
}

/* operator = : overwrite the data with v.                              */
template<class TYPE>
Cell2D<TYPE>& Cell2D<TYPE>::operator=(TYPE const& v)
{ data_ = v;}


} // namespace STK

#endif
// STK_CELL_H
