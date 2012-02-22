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
 * Project: stkpp::Arrays
 * Purpose:  Define the List1D class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_List1D.h
  * @brief This is an implementation of the Interface ITContainer1D.
 **/

#ifndef STK_LIST_H
#define STK_LIST_H

#include "../../Sdk/include/STK_ITContainer1D.h"
#include "../../Sdk/include/STK_IContainerRef.h"
#include "../../Arrays/include/STK_Display1D.h"
#include "STK_Cell.h"

namespace STK
{

/** @ingroup Arrays
  * @brief Templated One dimensional Horizontal List.
  *
  * A List1D is an implementation of the Interface ITContainer1D for list.
  **/
template<class TYPE>
class List1D : public ITContainer1D<TYPE, List1D<TYPE> >
             , public IContainerRef
{
  protected:
    CellHo<TYPE> *p_first_;       ///< First Element of the List
    CellHo<TYPE> *p_last_;        ///< Last Element of the List

  public:
    /** Default constructor : beg_ =1 and end_ =0.
     *  @param I range of the container
     **/
    List1D( Range const& I = Range())
          : ITContainer1D<TYPE, List1D>(I)
          , IContainerRef(false)
    { initialize(I); }

    /** Misc constructor, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    List1D( Range const& I, TYPE const& v)
          : ITContainer1D<TYPE, List1D >(I)
          , IContainerRef(false)
    { initialize(I);

      CellHo<TYPE>* p1  = p_first_;
      for ( Integer j=this->first(); j<=this->last(); j++)
      { (*p1) = v;             // overwrite the value of the current cell
        p1    = p1->getRight();   // Goto Right place
      }
    }

    /** Copy constructor
     *  @param T the container to copy
     **/
    List1D( const List1D<TYPE> &T)
          : ITContainer1D<TYPE, List1D>(T)
          , IContainerRef(false)
    {
      // initialize container
      initialize(T.range());
      // copy the container
      CellHo<TYPE>* p1  = p_first_;
      CellHo<TYPE>* pt1 = T.p_first_;

      for (Integer j=T.first(); j<=T.last(); j++)
      { (*p1) = pt1->data();   // write the value of the current cell
        p1    = p1->getRight();   // Goto Right
        pt1   = pt1->getRight();  // Goto Right
      }
    }

  protected:
    /** constructor by reference, ref_=1.
     *  This constructor does not copy physically the elements contained
     *  in the Container. The List1D is wrapped by a reference List1D reduced
     *  to the range J.
     *
     *  @param p_first the first cell of the container to wrap
     *  @param p_last the last cell of the container to wrap
     *  @param J range of the data to wrap
     **/
    List1D( CellHo<TYPE>* const & p_first
          , CellHo<TYPE>* const & p_last
          , Range const& J
          )
          : ITContainer1D<TYPE, List1D>(J)
          , IContainerRef(true)
          , p_first_(p_first)
          , p_last_(p_last)
    {
      // Current position
      currentPosition_  = this->first();
      p_current_ = p_first;
    }

  public:
    /** virtual dtor. */
    virtual ~List1D() { if (!this->isRef()) freeMem();}

    /** access to one element.
     *  @param pos index of the element
     *  @return a reference on the element @c pos
     **/
     inline TYPE& elt(Integer const& pos)
    {
      moveCurr(pos);
      return (p_current_->data());
    }

    /** access to one element const.
     *  @param pos index of the const element
     *  @return a constant reference on the element @c pos
     **/
    inline TYPE const& elt(Integer const& pos) const
    {
      moveCurr(pos);
      return (p_current_->data());
    }

    /** access to many elements.
     *  @param J the range of the elements
     *  @return a list with a reference to the elements in the given range
     **/
    inline List1D<TYPE> elt(Range const& J) const
    {
#ifdef STK_BOUNDS_CHECK
      if ((J.first()<this->first()))
      { throw out_of_range("List1D::elt(J) "
                           "J.first()<this->first()");
      }
      if ((J.last()>this->last()))
      { throw out_of_range("List1D::elt(J) "
                           "J.last()>this->last()");
      }
#endif
      // get J.first() cell adress
      moveCurr(J.first());
      CellHo<TYPE>* p_first = p_current_;
      // get J.last() cell adress
      moveCurr(J.last());
      CellHo<TYPE>* p_last = p_current_;
      // return the reference
      return List1D<TYPE>(p_first, p_last, J);
    }

    /** Clear the object. Memory is liberated and the
     *  range of the Container is set to 0:-1.
     **/
    void clear()
    {
      if (this->isRef()) return;   // Nothing to do for ref
      freeMem();        // Free mem
      this->setRange(); // Set to default the dimension
    }

    /** New first index for the object.
     *  @param beg new first index of the Container.
     **/
    void shift(Integer const& beg =1)
    {
      if (this->first() == beg) return;
#ifdef STK_DEBUG
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("List1D::shift(pos, n) "
                            "can't operate on references.");
      }
#endif
      //compute increment
      Integer inc = beg - this->first();
      this->incRange(inc);  // update this->range_()
      currentPosition_ += inc;         // update current position
    }

    /** Add n Elts to the container.
     *  @param n number of elements to add
     **/
    void pushBack(Integer const& n=1)
    {
      // if n==0 nothing to do
      if (n <= 0) return;
#ifdef STK_DEBUG
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("List1D::pushBack(pos, n) "
                            "can't operate on references.");
      }
#endif
      // If the container is empty : create it
      if (this->empty())
      {
        initialize(Range(this->first(), this->first()+ n -1));
      }
      else  // else adjust the beginning and the sizes
      {
        CellHo<TYPE> *p1, *p2;            // Auxiliary cells;
        try
        { p1 = new CellHo<TYPE>(p_last_);} // Create the end+1 cell
        catch (std::bad_alloc & error)   // if an alloc error is catched
        { // throw the Exception
          throw runtime_error("List1D::pushBack(n) "
                              "memory allocation failed.");
        }
        // if no error is intercepted
        p_last_->setRight(p1);             // Set the right ending cell
        for (Integer j=2; j<=n; j++)    // main loop for the other cells
        { try
          { p2 = new CellHo<TYPE>(p1);}  // try to allocate memory
          catch (std::bad_alloc & error) // if an alloc error occur
          {
            while ( p1 != p_last_)         // for all cells allocated
            { p2 = p1->getLeft();        // get the cell left
              delete p1;                 // delete the curent cell
              p1 = p2;                   // iterate
            }
            // set the original right side of cend
            p_last_->setRight(p_first_);
            // and throw an Exception
            throw runtime_error("List1D::pushBack(n) "
                                "memory allocation failed.");
          } // end catch
          // if no error is intercepted
          p1->setRight(p2);  // Set the right cell of the current cell
          p1 = p2;           // Set the current cell to the the next cell
        }
        p1->setRight(p_first_);    // the last cell point on the first cell
        p_first_->setLeft(p1);     // the first cell point on the last cell
        p_last_ = p1;             // the last cell adress
        this->incLast(n); // Update size of the container
      }
    }

    /** Insert element @c v in the range @c I of the List1D.
     *  @param I range of the index where to insert elements
     *  @param v the value tu insert
     **/
    void insert( Range const& I, TYPE const& v)
    {
      insertElt(I.first(), I.size());
      for (Integer i=I.first(); i<=I.last(); i++)
        elt(i) = v;
    }

    /** Insert n elts at the position pos of the container.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    void insertElt( Integer const& pos, Integer const& n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("List1D::insertElt(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check indices
      if (this->first() > pos)
      { throw out_of_range("List1D::insertElt(pos, n) "
                                "this->first() > pos");
      }
      if (this->last()+1 < pos)
      { throw out_of_range("List1D::insertElt(pos, n) "
                                "this->last()+1 < pos");
      }
#endif
      // Move the current position to j
      moveCurr(pos);
      CellHo<TYPE> *p0 = p_current_->getLeft(); // Get the j-1 cell
      CellHo<TYPE> *p1 = p0;                    // Auxiliary cell;
      // main loop for the other cells
      for (Integer j1=1; j1<=n; j1++)
      {
        CellHo<TYPE> *p2;        // Auxiliary cell;
        try
        { p2 = new CellHo<TYPE>(p1);}  // try to allocate memory
        catch (std::bad_alloc & error) // if an alloc error occur
        { while ( p1 != p0)            // for all cells allocated
          { p2 = p1;                   // get the cell left
            delete p1;                 // delete the curent cell
            p1 = p2->getLeft();        // iterate
          }
          p0->setRight(p_current_);
          // and throw an Exception
          throw runtime_error("List1D::insert(j, n) "
                              "memory allocation failed.");
        } // catch block
        // if no error is intercepted
        p1->setRight(p2);  // Set the right cell of the current cell
        p1 = p2;           // iterate
      }
      p1->setRight(p_current_);     // the last cell point on the first cell
      p_current_->setLeft(p1);      // the first cell point on the last cell
      if ( pos==this->first() )     // if the beginning was modified
      { p_first_ = p0->getRight();} // set new beginning
      this->incLast(n);             // Update the size of the container
      currentPosition_ +=n;         // Update the current position
    }

    /** Delete last elts of the container.
     *  @param n number of elts to delete
     **/
    void popBack(Integer const& n=1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("List1D::popBack() "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // if there is elts to erase
      if (this->size()<n)
      { throw out_of_range("List1D::popBack(n) "
                                "this->size() < n");
      }
#endif
      // erase elts with pos = end -n +1
      erase(this->last() - n +1, n);
    }

    /** Delete n elts at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void erase(Integer const& pos, Integer const& n=1)
    {
      // if n==0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { throw runtime_error("List1D::erase(pos, n) "
                                 "can't operate on references.");
      }
#ifdef STK_BOUNDS_CHECK
      // check bounds
      if (this->first() > pos)
      { throw out_of_range("List1D::erase(pos, n) "
                           "this->first() > pos");
      }
      if (this->last() < pos)
      { throw out_of_range("List1D::erase(pos, n) "
                           "this->last() < pos");
      }
      if (this->last() < pos+n-1)
      { throw out_of_range("List1D::erase(pos, n) "
                           "this->last() < pos+n-1");
      }
#endif
      // Move the current position to pos
      moveCurr(pos);
      CellHo<TYPE>* p2 = p_current_;  // get pos-th cell
      moveCurrLeft();                // set current to (pos-1)th position
      // delete n cells
      for (Integer l=1; l<=n; l++)
      { CellHo<TYPE>* p3 = p2->getRight();  // get right cell in p3
        delete p2;                          // delete current cell
        p2 = p3;                            // next
      }
      // If the last column have been erased update p_last_
      if (pos+n-1 == this->last()) { p_last_ = p_current_;}
      // Update the dimension of the container
      this->decLast(n);
      // If we have erased all cols
      if (this->size() == 0)
      { setDefault();}
      else
      { p2->setLeft(p_current_);        // p2 is the j+n cell
        p_current_->setRight(p2);       // p_current_ is on j-1 cell
        // If the first column has been erased
        if (pos == this->first())
        { p_first_  = p2;   // Set the new beg cell
          p_current_ = p2;   // p_current_
          currentPosition_++;       // and current position
        }
      }
    }

    /** Swapping the j1th column and the j2th column.
     *  @param j1 index of the first element to swap
     *  @param j2 index of the second element to swap
     **/
    void swap(Integer const& j1, Integer const& j2)
    {
#ifdef STK_BOUNDS_CHECK
      if (j1<this->first())
      { throw out_of_range("List1D::swap(j1, j2) "
                           "j1<this->first()");
      }
      if (j1>this->last())
      { throw out_of_range("List1D::swap(j1, j2) "
                           "j1>this->last()");
      }
      if (j2<this->first())
      { throw out_of_range("List1D::swap(j1, j2) "
                           "j2<this->first()");
      }
      if (j2>this->last())
      { throw out_of_range("List1D::swap(j1, j2) "
                           "j2>this->last()");
      }
#endif
      // get j1th value in aux
      moveCurr(j1);
      CellHo<TYPE> *p1 = p_current_;
      TYPE aux = p1->data();
      // set j2th value in j1th position
      moveCurr(j2);
      (*p1) = p_current_->data();
      // set j2th value to aux
      (*p_current_) = aux;
    }

    /** operator = : overwrite the Array1D with t.
     *  We resize the object if this and T does not have the same size
     *  but if they have the same size, we don't modify the range
     *  of the object.
     *  @param T the container to copy
     **/
     List1D<TYPE>& operator=(const List1D<TYPE> &T)
    {
      // We have to resize if this and T have not the same size
      // but if they have the same size, we don't scale the index
      if (this->size()!=T.size()) { this->resize(T.range());}

      /* copy without ovelapping.                                     */
      if (this->first() < T.first())
      { CellHo<TYPE> *p1 = p_first_, *pt1 = T.p_first_;
        for (Integer j=1; j<=this->size(); j++)
        { (*p1) = pt1->data();   // overwrite the value
          p1    = p1->getRight();   // Goto Right for this
          pt1   = pt1->getRight();  // Goto Right for T
        }
      }
      else
      { CellHo<TYPE> *p1 = p_last_, *pt1 = T.p_last_;
        for (Integer j=this->size(); j>=1; j--)
        { (*p1) = pt1->data();   // overwrite the value
          p1    = p1->getLeft();    // Goto Left for this
          pt1   = pt1->getLeft();   // Goto Left for T
        }
      }
      return *this;
    }

    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    List1D<TYPE>& operator=(TYPE const& v)
    {
      CellHo<TYPE>* p1 = p_first_;
      for (Integer j=1; j<=this->size(); j++)
      { p1->setData(v);      // overwrite the value of the current cell
        p1    = p1->getRight(); // Goto Right
      }
      return *this;
    }

  protected:
    /** Protected function for initialization.                        */
    void initialize(Range const& I)
    {
       // set new dimensions
      this->setRange(I);
      if (this->empty())
      {
        setDefault(); return;
      }
      // Allocate memory for the cells
      CellHo<TYPE> *p1, *p2;        // Auxiliary pointer for cells

      p1 = new CellHo<TYPE>();        // pointer on the first cell
      p_first_ = p1;                     // set the first cell
      // main loop for the other cells
      for (Integer j=this->first()+1; j<=this->last(); j++)
      { try
        { p2 = new CellHo<TYPE>(p1);}        // try to allocate memory
        catch (std::bad_alloc & error)       // if an alloc error occur
        { while (p1 != (CellHo<TYPE>*)NULL)  // for all cells allocated
          { p2 = p1->getLeft();              // get the cell left
            delete p1;                       // delete the cell
            p1 = p2;                         // and iterate
          }
          // set default
          setDefault();
          this->setRange();
          // and throw the Exception
          throw runtime_error("List1D::initialize(beg, end) "
                              "memory allocation failed.");
        }
        // if no error is catched
        p1->setRight(p2);      // Set the right cell
        p1 = p2;               // and iterate
      }
      p_last_ = p1;              // Set the last cell
      p_last_->setRight(p_first_);  // the last cell point on the first cell
      p_first_->setLeft(p_last_);   // the first cell point on the last cell

      currentPosition_  = this->first();    // current position is first position
      p_current_ = p_first_;            // Current cell is first cell
    }

    /** Protected function for deallocation.                         */
    void freeMem()
    {
      if (this->isRef()) return;   // Nothing to do for ref
      CellHo<TYPE> *p2, *p1 =p_first_;   // Auxiliary pointers for cells
      // for all cells
      for (Integer j=this->first(); j<=this->last(); j++)
      { p2 = p1->getRight();               // get the right cell
        delete p1;                         // delete the curent cell
        p1 = p2;                           // and iterate
      }
      setDefault();
    }

    /** Protected function for setting members to default.
     **/
    void setDefault()
    { p_first_  = (CellHo<TYPE>*)NULL;
      p_last_  = (CellHo<TYPE>*)NULL;
      p_current_ = (CellHo<TYPE>*)NULL;
      currentPosition_  = this->first();
    }

  private:
    /**
     *  Current position of pointer p_current_ in the List1D
     */
    mutable Integer currentPosition_;

    /**
     *  Current position pointed in the List1D
     */
    mutable CellHo<TYPE> *p_current_;

    /** Private function for moving the current position :
     *  it's very low level manipulation functions, and no check is done
     *  at this level.
     *  Note that it is constant functions, as the object itself is not
     *  modified, the modified members are mutable.
     *  @{
     **/
    /** Move Current position to left
     **/
    void moveCurrLeft()  const
    {
      p_current_ = p_current_->getLeft();
      currentPosition_--;
    }
    /** move Current position to right
     **/
    void moveCurrRight() const
    {
      p_current_ = p_current_->getRight();
      currentPosition_++;
    }
    /** Move the current position to pos
     *  @param pos the position to move
     **/
    void moveCurr(Integer const& pos) const
    {
      // if j is greater than the current
      if (pos>currentPosition_)
      {
        if ((pos-currentPosition_) <= (this->last()-pos))    // if j is near the current
          for( ;currentPosition_!=pos; ) moveCurrRight(); // move to right
        else                                 // else we are near the end
          for( currentPosition_ = this->last(), p_current_ = p_last_
             ; currentPosition_!=pos
             ; )
             moveCurrLeft();
      }
      else  // else j is less than the current
      {
        if ((currentPosition_-pos) <= (pos-this->first()))  // if j is near the current
          for( ;currentPosition_!=pos; ) moveCurrLeft();  // move to left
        else                                 // else we are near the beg
          for( currentPosition_ = this->first(), p_current_ = p_first_
             ; currentPosition_!=pos
             ; )
             moveCurrRight();
      }
    }
 };

/** ostream for List1D.
 *  @param s the output stream
 *  @param V the List1D to write
 **/
template<class TYPE>
ostream& operator<<(ostream& s, const List1D<TYPE>& V)
{ return out1D(s,V);}

} // namespace STK

#endif
// STK_LIST_H
