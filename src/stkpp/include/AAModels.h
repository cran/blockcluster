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

    Contact : Serge.Iovleff@stkpp.org                                  */

/*
 * Project:  stkpp::AAModels
 * Purpose:  Main include file for the AAModels project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file AAModels.h
 *  @brief This file include all the header files of the project AAModels.
 *
 *  @defgroup AAModels Auto-Associative Models.
 *  @brief The project AAM provides classes and tools for unsupervised learning
 *  and data analysis using Auto-Associative models.
 *
 * A function \f$ g \f$ is an auto-associative function
 * of dimension d if it is a map from \f$ \mathbb{R}^p \f$ to
 * \f$ \mathbb{R}^p \f$ that can be written \f$ g=r\circ p \f$ where
 * p (the ``Reduction'') is a map from \f$ \mathbb{R}^p \f$ to
 * \f$ \mathbb{R}^d \f$ (generally d<p) and r (the ``Regression'') is
 * a map from \f$ \mathbb{R}^d \f$ to \f$ \mathbb{R}^p \f$ .
 **/

#ifndef AAMODELS_H
#define AAMODELS_H

#include "../projects/AAModels/include/STK_IAAModel.h"
#include "../projects/AAModels/include/STK_GaussianAAModel.h"
#include "../projects/AAModels/include/STK_LinearAAModel.h"

#endif // AAMODELS_H
