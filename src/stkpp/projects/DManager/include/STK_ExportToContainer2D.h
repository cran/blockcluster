/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * created on: 9 juin 2011
 * Purpose:  Create an utility class in order to transfer the Data from
 * a DataFrame in a TContainer2D.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_ExportToContainer2D.h
 *  @brief In this file we define the ExportToContainer2D class.
 **/

#ifndef STK_EXPORTTOCONTAINER2D_H
#define STK_EXPORTTOCONTAINER2D_H

#include "STK_DataFrame.h"
#include "../../Sdk/include/STK_ITContainer2D.h"

namespace STK
{

/** @brief The ExportToContainer2D class allow to export the data of some @c TYPE
 *  stored in a @c DataFrame to be exported in a 2D container.
 *
 *  The 2D container is created on the stack and will be deleted with the
 *  @c ExportToContainer2D structure. It is possible to release the 2D container
 *  by calling explicitly the @c release()
 */
template<class TYPE, class TContainer2D>
class ExportToContainer2D
{
  public:
    /** Constructor. Create an empty 2D Container and export all the variables
     * of type @c TYPE in it.
     *  @param df the DataFrame to export
     */
    ExportToContainer2D( DataFrame const& df)
                       : p_data_(new TContainer2D())
    {
      // for each field Try a type conversion
      for(Integer iVar = df.firstCol(); iVar<=df.lastCol(); iVar++)
      {
        IVariable* const p_var = df.elt(iVar);
        // if there is a variable
        if (p_var)
        {
          // check the type of the variable
          if (p_var->getType() == IdTypeImpl<TYPE>::returnType())
          {
            Variable<TYPE>* p_variable = static_cast<Variable<TYPE>* >(p_var);
            p_data_->pushBackCol(*p_variable);
          }
        }
      }
    }

    /** virtual destructor */
    virtual ~ExportToContainer2D()
    { if (p_data_) delete p_data_;}

    /** Accesor. get the 2D container. This method is not constant
     * in order to the user to modified directly the 2D container.
     * @return a ptr on the 2D container constructed
     **/
    inline TContainer2D* p_container2D() { return p_data_->asPtrLeaf();}

    /** release the Array2D. It will be freed by the user. */
    inline void release() { p_data_ =0;}

    /** remove NA values from the Array2D
     * @param byRow @c true if the user want to delete the rows with NA values,
     * @c false if the user want to remove the column with NA values.
     **/
    inline void eraseNAValues(bool byRow)
    {
      // check if there exists data
      if (!p_data_) return;
      // get the first index
      const Integer firstCol = p_data_->firstCol(), firstRow = p_data_->firstRow();
      // get the last index
      Integer lastCol = p_data_->lastCol(), lastRow = p_data_->lastRow();
      if (byRow)
      {
        // loop on the rows
        for (Integer i= firstRow; i <= lastRow; ++i)
        {
          bool asNA = false;
          // loop on the element of the row
          for (Integer j= firstCol; j <= lastCol; ++j)
          {
            if (Arithmetic<TYPE>::isNA(p_data_->elt(i,j)))
            {
              asNA = true;
              break;
            }
          }
          // remove current row
          if (asNA)
          {
            p_data_->eraseRows(i);
            i--;
            lastRow--;
          }
        }
      }
      else
      {
        // loop on the column
        for (Integer j= firstCol; j <= lastCol; ++j)
        {
          bool asNA = false;
          // loop on the element of the column
          for (Integer i= firstRow; i <= lastRow; ++i)
          {
            if (Arithmetic<TYPE>::isNA(p_data_->elt(i,j)))
            {
              asNA = true;
              break;
            }
          }
          // remove current column
          if (asNA)
          {
            p_data_->eraseCols(j);
            j--;
            lastCol--;
          }
        }
      }
    }

  private:
    /** ptr on the 2D Container the class will export. */
    ITContainer2D<TYPE, TContainer2D>* p_data_;
};

}

#endif /* STK_EXPORTTOCONTAINER2D_H */
