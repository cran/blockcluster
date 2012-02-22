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
 * Project:  stkpp::
 * created on: 17 févr. 2012
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Arrays_Util.h
 *  @brief In this file we define utilities for the Array classes.
 **/


#ifndef STK_ARRAYS_UTIL_H

#ifndef BLOCKSIZE
#define BLOCKSIZE 8
#endif

#ifndef PANELSIZE
#define PANELSIZE 256
#endif

namespace STK
{

namespace Arrays
{

/** Define the different type of Array that can be handle */
enum TypeArray
{
  unknown_ =0, //!< unknown_ unknown type
  dense_,      //!< denses_  dense array
  sparse_      //!< sparse_  sparse  a one row array
};

/** Define the different structure of Array that can be handle */
enum StructureArray
{
  unknown_ =0,//!< unknown_ unknown type
  vector_,    //!< vector_  a one column array
  point_,     //!< point_   a one row array
  rect_,      //!< rect_    a rectangular array
  tri_sup_,   //!< tri_sup_ a triangular superior array
  tri_inf_,   //!< tri_inf_ a triangular inferior array
  square_,    //!< square_  a square array
  band_,      //!< band_    a band array
  diag_       //!< diag_    a diagonal array
};

/**Define the Storage Orientation of the containe */
enum StorageOrientation
{
  unknown_ =0,//!< unknown_ unknown storage
  by_row_,    //!< by_row_  storage by row
  by_col_     //!< by_col_  storage by column
};

} // namespace Arrays

} // namespace STK

#define STK_ARRAYS_UTIL_H


#endif /* STK_ARRAYS_UTIL_H */
