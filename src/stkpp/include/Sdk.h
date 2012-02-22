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
 * Project:  stkpp::Sdk
 * Purpose:  main include file for Sdk project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file Sdk.h
 *  @brief This file include all the other header files of the
 *  project Sdk.
 *
 **/

/**
 * @defgroup Sdk Software Development Kit.
 * @brief The Sdk project propose a set of high level Interfaces that are
 * implemented in the stk++ projects.
 *
 * In Sdk we define the pure Interface classes than can be used throughout the
 * stkpp whole project. The aim is to unified the syntax and the treatment
 * realized by the statistical methods.
 *
 * The sub-project TContainer define the API that have to be used by concrete
 * implementation of the 1D and 2D containers.
 **/

/** @ingroup Sdk
 *  @defgroup TContainer  Subproject Sdk::TContainer
 *  In the TContainer subproject, we define and implement
 *  <ul>
 *     <li> abstract and templated interface base classes for storing arbitrary
 *     data in one dimensional and two dimensional containers.
 * </ul>
 **/

#ifndef SDK_H
#define SDK_H

#include "../projects/Sdk/include/STK_IRecursiveTemplate.h"

/* Interface for all containers that can be constructed as a reference.*/
#include "../projects/Sdk/include/STK_IContainerRef.h"

/* Interface Base for all containers. */
#include "../projects/Sdk/include/STK_IContainerBase.h"

/* Interface for all 1D and 2D containers. */
#include "../projects/Sdk/include/STK_IContainer1D.h"
#include "../projects/Sdk/include/STK_IContainer2D.h"

/* Interface for all 1D and 2D containers storing a single type of data. */
#include "../projects/Sdk/include/STK_ITContainer1D.h"
#include "../projects/Sdk/include/STK_ITContainer2D.h"

/* Interface for all runners */
#include "../projects/Sdk/include/STK_IRunner.h"

#include "../projects/Sdk/include/STK_IRunnerPtr2D.h"


#endif  /* SDK_H */
