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
 * Purpose:  Implement the DataFrame (Table) class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
**/

/** @file STK_DataFrame.cpp
 *  @brief In this file we implement the DataFrame (Table) class.
 **/

#include "../include/STK_DataFrame.h"
#include "../include/STK_ExportToCsv.h"

namespace STK
{
/* Default constructor . */
DataFrame::DataFrame() : _BaseList()
                       , IContainer2D(Range(), Range())
{ ;}

/* Copy constructor                                                         */
DataFrame::DataFrame( DataFrame const& T, bool ref)
                    : _BaseList(T)
                    , IContainer2D(T)
{
  // the adress of the variables are copied in List1D
  // but we need to copy explicitly the data
  for (Integer j=first(); j<=last(); j++) // for all columns
    if (T[j])                            // if there is data
      at(j) = T[j]->clone(ref); // set the adress of a clone
}

/* destructor.                                                              */
DataFrame::~DataFrame()
{
  // free the rows as the destructor of _BaseList will not free the mem
  freeRows();
}

/* clear the object.                                                  */
void DataFrame::clear()
{
  freeMem();
  // set default range for list
  _BaseList::setRange();
  // set default range for container2D
  IContainer2D::setRange();
}

/* Operator = : overwrite the DataFrame with T. */
DataFrame& DataFrame::operator=(DataFrame const& T)
{
  // Resize if necessary.
  if (sizeHo() != T.sizeHo()) _BaseList::resize(T.rangeHo());

  // Copy without overlapping.
  if ((T.first()>=first()))
  {
    for (Integer jt=T.first(), j=first(); jt<=T.last(); j++, jt++)
    {
      // clear old mem if any
      if (at(j)) delete at(j);
      // if there is a variable, create a copy
      if (T[jt]) at(j) = T[jt]->clone();
      else       at(j) = (IVariable*)NULL;
    }
  }
  else
  {
    for (Integer jt=T.last(), j=last(); jt>=T.first(); j--, jt--)
    {
      // clear old mem if any
      if (at(j)) delete at(j);
      // if there is a variable, create a copy
      if (T[jt]) at(j) = T[jt]->clone();
      else       at(j) = (IVariable*)NULL;
    }
  }
  return *this;
}

/* New beginning index for the object. */
void DataFrame::shift(Integer const& cbeg)
{
  // list1D shift
  _BaseList::shift(cbeg);
  // IContainer2D shift for Col
  setFirstHo(cbeg);
}

/* New beginning index for the object. */
void DataFrame::shift(Integer const& rbeg, Integer const& cbeg)
{
  // if there is something to do
  if ((rbeg - firstRow() != 0)||(cbeg - firstCol() != 0))
  {
    // list1D shift
    shift(cbeg);
    // For each column update Variable
    for (Integer j=first(); j<=last(); j++)
      if (at(j)) { at(j)->shift(rbeg);}
    // update range of the rows
    setFirstVe(rbeg);
  }
}


/* Del last column of the container.                                     */
void DataFrame::popBackCols(Integer const& n)
{
  // if n<=0 nothing to do
  if (n<=0) return;
  // if there is columns to erase
  if (sizeHo()<n)
  { throw out_of_range("DataFrame::popBackCols(n) "
                       "sizeHo() < n");
  }
  // for all columns, delete variables
  for (Integer j=lastCol() - n +1; j<=lastCol(); j++)
    if (at(j)) delete at(j);
  // popBack() of List1D
  _BaseList::popBack(n);
  // update IContainer2D
  decLastHo(n);
  // if it was the last elt, free mem
  if (this->sizeHo() == 0) freeMem();
}

/* Delete n columns at the nuber pos of the container. */
void DataFrame::eraseCols( Integer const& pos, Integer const& n)
{
  // if n<=0 nothing to do
  if (n<=0) return;
  // check conditions
  if (pos<first())
  { throw out_of_range("DataFrame::eraseCols(pos, n) "
                       "pos < first()");
  }
  if (pos>last())
  { throw out_of_range("DataFrame::eraseCols(pos, n) "
                       "pos > last()");
  }
  if (last() < pos+n-1)
  { throw out_of_range("DataFrame::eraseCols(pos, n) "
                       "last() < pos+n-1");
  }
  // for all columns, delete variables
  for (Integer j=pos+n-1; j>=pos; j--)
    if (at(j)) { delete at(j);}
  // delete elements of the List1D
  erase(pos, n);
  // update rangeHo_
  decLastHo(n);
  // if it was the last col, free mem
  if (this->sizeHo() == 0) freeMem();
}

/* Insert variable V at the position i to the container.              */
void DataFrame::insertVariable(Integer const& pos, IVariable* const & V)
{
  // List1D
  List1D<IVariable*>::insertElt(pos);
  at(pos) = V;
  // the variable have to be in the same range
  at(pos)->shift(firstRow());
  // update horizontal range (the number of column)
  incLastHo(1);

  // update rows with NA values
  Integer inc = sizeVe() - V->size();
  if (inc == 0) return; // same size
  if (inc > 0) // V has less rows
  { // put NA values to the inserted column
    at(pos)->pushBackNAValues(inc);
  }
  else
  { // put NA values to the other columns
    for (Integer i=this->firstCol(); i <pos; i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    for (Integer i=pos+1; i <=this->lastCol(); i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    // update LastVe
    incLastVe(-inc);
  }
}

/* Merge this with a dataframe (horizontally). */
void DataFrame::pushBackVariable( IVariable* const &V)
{
  // List1D
  push_back(V);
  // update horizontal range (the number of column)
  incLastHo(1);
  // adjust the first index of the inserted variable
  at(lastCol())->shift(firstRow());
  // update rows with NA values
  Integer inc = sizeVe() - V->size();
  if (inc == 0) return; // same size
  if (inc > 0) //V has less rows
  { // put NA values to the inserted columns
    at(lastCol())->pushBackNAValues(inc);
  }
  else
  { // put NA values to the oter columns
    for (Integer i=this->first(); i <lastCol(); i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    // update LastVe
    incLastVe(-inc);
  }
}

/* Insert the DatatFrame D at the column pos to the container.           */
void DataFrame::insertDataFrame( Integer const& pos, const DataFrame& D)
{
  // List1D
  insertElt(pos, D.sizeHo());
  // insert all columns of D
  for (Integer i = D.firstCol(), icol = pos; i <=D.lastCol(); i++, icol++)
  {
    if (D.at(i))
    {
      at(icol) = D.at(i)->clone();
      at(icol)->shift(firstRow());
    }
  }
  // update LastHo
  incLastHo(D.sizeHo());
  // update rows with NA values
  Integer inc = sizeVe() - D.sizeVe();
  if (inc == 0) return; // same size
  if (inc > 0) // D has less rows
  { // put NA values to the inserted columns
    for (Integer i= pos+D.sizeHo()-1; i >=pos; i--)
      if (at(i)) { at(i)->pushBackNAValues(inc);}
  }
  else
  { // put NA values to the oter columns
    for (Integer i=this->firstCol(); i <pos; i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    for (Integer i=pos+D.sizeHo(); i <=this->last(); i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    // update LastVe
    incLastVe(-inc);
  }
}

/* Merge this with a dataframe (horizontally).
*/
void DataFrame::pushBackDataFrame( DataFrame const &D)
{
  // compute pos
  Integer pos(last()+1);
  // List1D
  pushBack(D.sizeHo());
  // insert all columns of D
  for (Integer i = D.first(), icol = pos; i <=D.last(); i++, icol++)
  {
    if (D.at(i))
    {
      at(icol) = D.at(i)->clone();
      at(icol)->shift(firstRow());
    }
  }
  // update LastHo
  incLastHo(D.sizeHo());
  // update rows with NA values
  Integer inc = sizeVe() - D.sizeVe();
  if (inc == 0) return; // same size
  if (inc > 0) // D has less rows
  { // put NA values to the inserted columns
    for (Integer i= last(); i >=pos; i--)
      if (at(i)) { at(i)->pushBackNAValues(inc);}
  }
  else
  { // put NA values to the oter columns
    for (Integer i=this->first(); i <pos; i++)
      if (at(i)) { at(i)->pushBackNAValues(-inc);}
    // update LastVe
    incLastVe(-inc);
  }
}

/* Add columns to the container.                                         */
void DataFrame::pushBackCols(Integer const& n)
{
  // if n<=0 nothing to do
  if (n <= 0) return;
  // add n columns to list1D
  insert(Range(last()+1, last()+n), (IVariable*)NULL);
  // update IContainer2D
  incLastHo(n);
}

/* Insert columns at the specified position to the container.            */
void DataFrame::insertCols( Integer const& pos, Integer const& n)
{
  if (n <= 0) return;        // if n<=0 nothing to do
#ifdef STK_BOUNDS_CHECK
  // check conditions
  if (pos<first())
  { throw out_of_range("DataFrame::insertCols(pos, n) "
                            "pos<first()");
  }
  if (pos>last())
  { throw out_of_range("Dataframe::insertCols(pos, n) "
                            "pos>last()");
  }
#endif
  // insert n elements in list1D
  insert(Range(pos, pos+n-1), (IVariable*)NULL);
  // update IContainer2D
  incLastHo(n);
}

/* Add n rows to the container.                                       */
void DataFrame::pushBackRows(Integer const& n)
{
  // if n<=0 nothing to do
  if (n<=0) return;
  // for each column append row
  for (Integer j=first(); j<=last(); j++)
  {
    if (at(j))
    { at(j)->pushBack(n);}
  }
  // update range of the container
  incRangeVe(n);
}

/* Insert n rows at the ith position of the container.                */
void DataFrame::insertRows( Integer const& pos, Integer const& n)
{
  // if n<=0 nothing to do
  if (n<=0) return;
#ifdef STK_BOUNDS_CHECK
  if (firstRow() > pos)
  { throw out_of_range("DataFrame::insertRows(pos, n) "
                            "firstRow() > pos");
  }
  if (lastRow()+1 < pos)
  { throw out_of_range("DataFrame::insertRows(pos, n) "
                            "lastRow()+1 < pos");
  }
#endif
  // insert rows to each variables
  for (Integer j=first(); j<=last(); j++)
  {
    // if there is a variable
    if (at(j))
    { at(j)->insertElt(pos, n);}
  }
  // update rangeVe_
  incLastVe(n);
}

/* Dell last row of the container.                                    */
void DataFrame::popBackRows(Integer const& n)
{
  if (sizeVe() < n)
  { throw out_of_range("DataFrame::popBackRows(n) "
                            "sizeVe() < n");
  }
  // del last row to each variable
  for (Integer j=first(); j<=last(); j++)
    if (at(j)) { at(j)->popBack(n);}
  // update rangeVe_
  decLastVe(n);
}

/* Dell n rows at the ith position to the container.                  */
void DataFrame::eraseRows( Integer const& pos, Integer const& n)
{
  // if n<=0 nothing to do
  if (n<=0) return;
#ifdef STK_BOUNDS_CHECK
  if (firstRow() > pos)
  { throw out_of_range("DataFrame::eraseRows(pos, n) "
                            "firstRow() > pos");
  }
  if (lastRow() < pos)
  { throw out_of_range("DataFrame::eraseRows(pos, n) "
                            "lastRow() < pos");
  }
  if (lastRow() < pos+n-1)
  { throw out_of_range("DataFrame::eraseRows(pos, n) "
                            "lastRow() < pos+n-1");
  }
#endif
  // for each variable erase elts
  for (Integer j=first(); j<=last(); j++)
    if (at(j)) { at(j)->erase(pos, n);}
  // update rangeVe_
  decLastVe(n);
}

/* Protected function for memory deallocation.                       */
void DataFrame::freeMem()
{
  // liberate variables
  freeRows();
  // call base freeMem
  _BaseList::freeMem();
  // set range to default
  setRangeVe();
  setRangeHo();
}

/* Protected function for rows memory deallocation.                  */
void DataFrame::freeRows()
{
  // for all columns
  for (Integer j=first(); j<=last(); j++)
    if (at(j))          // if there is mem allocated
    {
      delete at(j);     // erase
      at(j) = 0;        // set default
    }
  // set default range
  setRangeVe();
}

// write a selection
void DataFrame::writeDataFrame( ostream& os, Integer const& left
                              , Integer const& right
                              ) const
{
  // Export  to csv the DataFrame
  ExportToCsv csv(*this);
  // get the csv
  ReadWriteCsv* pData = csv.p_readWriteCsv();
  // set delimiters to blank
  pData->setDelimiters(STRING_BLANK);
  // write the csv
  pData->writeSelection(os, firstRow(), lastRow(), left, right);
}

/* Print a DataFrame.                                                 */
ostream& operator<< (ostream& s, const DataFrame& V)
{
  s << std::right;
  V.writeDataFrame(s, V.firstCol(), V.lastCol());

  return s;
}

} // Namespace STK
