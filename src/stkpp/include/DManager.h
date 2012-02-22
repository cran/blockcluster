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
 * Purpose:  Include all files of the DManager project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file DManager.h
 *  @brief This file include all the other header files of the
 *  project DManager.
 **/

/**
 *  @defgroup DManager Data Management
 *  @brief The DManager project propose classes and functions for managing the
 *  data.
 *
 * The DManager project contains all the class related to data management.
 * It provides
 * <ul>
 *    <li> an abstract base class for statistical variables,
 *    <li> a templated implementation for arbitrary data,
 *    <li> a Dataframe (Table) class,
 *    <li> classes for read and write csv and (TODO) dbf files,
 *    <li> various  classes for importing/exporting data from different
 *    containers,
 *    <li> methods for sorting one dimensional containers,
 *    <li> abstract classes for reading pages of options.
 *</ul>
 *
 **/

#ifndef DMANAGER_H
#define DMANAGER_H

/* Utilities used in the DManager Project */
#include "../projects/DManager/include/STK_DManager_Util.h"
#include "../projects/DManager/include/STK_Import_Util.h"

/* Interface Variable class and Variable class */
#include "../projects/DManager/include/STK_IVariable.h"
#include "../projects/DManager/include/STK_Variable.h"

/* DataFrame. */
#include "../projects/DManager/include/STK_Cell.h"
#include "../projects/DManager/include/STK_List1D.h"
#include "../projects/DManager/include/STK_DataFrame.h"

/* Export data from a Data Frame to an Array2D */
#include "../projects/DManager/include/STK_ExportToContainer2D.h"

/* HeapSort utilities. */
#include "../projects/DManager/include/STK_HeapSort.h"

/* main classes for option files. */
#include "../projects/DManager/include/STK_Option.h"
#include "../projects/DManager/include/STK_IPage.h"
#include "../projects/DManager/include/STK_ReadWritePages.h"

/* main classes for managing Csv data. */
#include "../projects/DManager/include/STK_ReadWriteCsv.h"
#include "../projects/DManager/include/STK_ExportToCsv.h"
#include "../projects/DManager/include/STK_ImportFromCsv.h"

#endif  /* DMANAGER_H */
