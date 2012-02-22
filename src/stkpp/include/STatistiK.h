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
 * Project:  stkpp::STatistiK
 * Purpose:  Primary include file for STatistiK project.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STatistiK.h
 *  @brief This file include all the header files of the project STatistiK.
 * 
 * @defgroup STatistiK STatistiK (Tools for usual statistics).
 * @brief The StatistiK project contains the main tools for computing the usual
 * statistics.
 *
 * It is divided in various sub-projects:
 * - @ref Analysis
 * - @ref Laws
 * - @ref StatDesc
 **/

/** @ingroup STatistiK
 * @defgroup Analysis The Analysis sub-project
 * In this sub-project, we compute usual and special functions: gamma,
 * gammmaRatio, betaRatio,...
 *
 * @ingroup STatistiK
 * @defgroup Laws     The Laws sub-project.
 * In this sub-project, we compute and simulate usual probabilities laws:
 * normal, binomial, cacuchy,...
 *
 * @ingroup STatistiK
 * @defgroup StatDesc The descriptive statistics sub-project.
 * In this sub-project, we compute the usual descriptive statistics of
 * variables.
 *
 * 
 * @ingroup STatistiK
 * @namespace STK::Const
 * @brief This is the namespace containing the usual mathematical constants
 *
 * @ingroup STatistiK
 * @namespace STK::Law
 * @brief This is the namespace containing the class handling the usual
 * probabilities Laws.
 * The namespace Law is the domain space of the probabilities laws
 * like normal law, gamma law, binomial law...
 *
 * @ingroup STatistiK
 * @namespace STK::Funct
 * @brief The namespace Funct enclose all usual and special functions.
 * The namespace Funct is the domain space of the special function
 * like gamma function, beta function, incomplete gamma function,
 * incomplete beta function... It include also some useful raw
 * functions like log1p...
 * 
 * @ingroup STatistiK
 * @namespace STK::Stat
 * @brief this is the namespace for the statistical treatment.
 * The namespace Stat is the domain space for the usual statistical
 * treatment of the variable like mean, variance, covariance, ...
 *
 **/

#ifndef STATISTIK_H
#define STATISTIK_H

// templated generic algorithms
#include "../projects/STatistiK/include/STK_Algo.h"

// namespace Const
// Mathematical constant
#include "../projects/STatistiK/include/STK_Const_Math.h"

// namespace Funct
// usual fonctions
#include "../projects/STatistiK/include/STK_Funct_util.h"
// raw functions
#include "../projects/STatistiK/include/STK_Funct_raw.h"
// gamma function
#include "../projects/STatistiK/include/STK_Funct_gamma.h"
// gamma ratio function
#include "../projects/STatistiK/include/STK_Funct_gammaRatio.h"
// beta Ratio function
#include "../projects/STatistiK/include/STK_Funct_betaRatio.h"

// namespace Law
// probabilities laws
#include "../projects/STatistiK/include/MersenneTwister.h"
#include "../projects/STatistiK/include/STK_RandBase.h"

#include "../projects/STatistiK/include/STK_Law_ITUnivariate.h"
#include "../projects/STatistiK/include/STK_Law_Normal.h"
#include "../projects/STatistiK/include/STK_Law_ILawBase.h"

#include "../projects/STatistiK/include/STK_Law_ITMultivariate.h"
#include "../projects/STatistiK/include/STK_Law_MultivariateNormal.h"

// namespace Stat
// Univariate Statistics
#include "../projects/STatistiK/include/STK_Stat_Univariate.h"
#include "../projects/STatistiK/include/STK_Stat_UnivariateReal.h"

// bivariate Statistics
#include "../projects/STatistiK/include/STK_Stat_Bivariate.h"
#include "../projects/STatistiK/include/STK_Stat_BivariateRealReal.h"

// Multivariate Statistics
#include "../projects/STatistiK/include/STK_Stat_Multivariate.h"
#include "../projects/STatistiK/include/STK_Stat_MultivariateReal.h"

// perform the usual transformations on variables
#include "../projects/STatistiK/include/STK_Stat_Transform.h"

#endif /*STATISTIK_H*/

