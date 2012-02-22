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
 * Purpose:  Define the DataFrame class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_DataFrame.h
 *  @brief In this file we declare the DataFrame class.
 **/

#ifndef STK_DATAFRAME_H
#define STK_DATAFRAME_H

#include "../../Sdk/include/STK_IContainer2D.h"

#include "STK_List1D.h"
#include "STK_IVariable.h"

namespace STK
{
/** @ingroup DManager
  * @brief DataFrame is a List of Variable with the same number of rows.
  *  This is thus also a 2D container.
  *
  * A DataFrame inherit from the class @c List1D and and @c IContainer2D.
  * It is essentially a List, except that each column (the variables)
  * possess the same dimension.
  *
  * Each Cell of the List1D contain a pointer on a Variable.
 **/
class DataFrame : protected List1D<IVariable* >
                , public IContainer2D
{
  protected:
    /** Type for the list container. */
    typedef List1D<IVariable*> _BaseList;

  public:
    /** Default Constructor, empty table. */
    DataFrame();

    /** Copy constructor If ref is true, only references of the variables
     *  are copied into the DataFrame.
     *  @param T the DataFrame to copy
     *  @param ref true if we want to wrap the variables of T
     **/
    DataFrame( DataFrame const& T, bool ref = false);

    /** Destructor. */
    virtual ~DataFrame();

    /** Clear the object. */
    void clear();

    /** access to an element. Set the method elt as a public method. */
    inline IVariable* & elt(Integer const& i)
    { return _BaseList::elt(i);}

    /** access to a constant element. Set the method elt as a public method. */
    inline IVariable* const & elt(Integer const& i) const
    { return _BaseList::elt(i);}

    /** Operator = : overwrite the DataFrame with T.                  */
    DataFrame& operator=(DataFrame const& T);

    /** New beginning index for the object.                          */
    void shift(Integer const& cbeg =1);

     /** New beginning index for the object.                          */
    void shift(Integer const& rbeg, Integer const& cbeg);

    /** Del n column of the container.                                   */
    void popBackCols(Integer const& n);

    /** Del n cols at the position of the container.                  */
    void eraseCols( Integer const& pos, Integer const& n=1);

    /** Swapping the j1th column and the j2th column.
     *  @param j1 index of the first column to swap
     *  @param j2 index of the second column to swap
     * */
    inline void swapCols(Integer j1, Integer j2)
    { _BaseList::swap(j1, j2);}

    /** Dell last rows of the container.
     *  @param n number of rows to delete
     **/
    void popBackRows(Integer const& n);

    /** Delete n rows at the position @c pos to the container.
     *  @param pos position of the rows to delete
     *  @param  n  number of rows to delete
     **/
    void eraseRows( Integer const& pos, Integer const& n=1);

    /** Insert a Vartiable at the specified position to the container.
     *  @param pos the position in the container
     *  @param V the Variable to insert
     **/
    void insertVariable( Integer const& pos, IVariable* const & V);

    /** Append a DataFrame back.
     *  @param V The variable to append to the DataFrame
     **/
    void pushBackVariable( IVariable* const & V);

    /** Append a DataFrame front.                                     */
    inline void pushFrontVariable( IVariable* const & V)
    { insertVariable(firstCol(), V);}

    /** Insert a DataFrame at the specified position
     *  to the container.
     **/
    void insertDataFrame( Integer const& pos, const DataFrame& D);

    /** Append a DataFrame back.                                      */
    void pushBackDataFrame( DataFrame const &D);

    /** Append a DataFrame front.                                     */
    inline void pushFrontDataFrame( DataFrame const &D)
    { insertDataFrame(firstCol(), D);}

    /** write a DataFrame to the output stream os.                    */
    void writeDataFrame( ostream  &os
                       , Integer const& left
                       , Integer const& right
                       ) const;
  protected:
    /** function for memory deallocation.                            */
    void freeMem();

    /** function for row memory deallocation.                        */
    void freeRows();

    /** Add cols to the container.                                    */
    void pushBackCols(Integer const& n=1);

    /** Insert cols at the specified position to the container.       */
    void insertCols( Integer const& pos, Integer const& n=1);

    /** Add n rows to the container.                                  */
    void pushBackRows(Integer const& n=1);

    /** Insert n rows at the ith position of the container. */
    void insertRows( Integer const& pos, Integer const& n =1);
};

/** Print a DataFrame.
 * @param s ooutput stream
 * @param V the Dataframe to write
 **/
ostream& operator<< (ostream& s, const DataFrame& V);

} // Namespace STK

#endif
//STK_DATAFRAME_H
