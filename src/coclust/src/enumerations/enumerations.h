
/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2011  Parmeet Singh Bhatia

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
*/

/*
 * Project:  cocluster
 * created on: Dec 26, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file enumerations.h
 *  @brief Defines all the enumerations used in CoClust Project.
 **/


#ifndef ENUMERATIONS_H_
#define ENUMERATIONS_H_

/*
 * Enumeration for Data-type.
 */
enum DataType
{
    Binary = 1,
    Contingency = 2,
    Continuous = 3
};

/*
 * Enumeration for Algorithms.
 */
enum Algorithm
{
    BEM2 = 1,
    BCEM = 2,
    XEMStrategy = 3,
    XCEMStrategy = 4
};

/**
 * Enumeration for Stopping Criteria.
 */
enum StopCriteria
{
    Parameter = 1,
    Likelihood = 2
};

/**
 * Enumeration for Model Initialization.
 */
enum Initialization
{
    e_CEMInit = 1,
    e_FuzzyCEMInit = 2,
    e_RandomInit = 3
};

/**
 * Enumeration for all data Models.
 */
enum Model
{
    pi_rho_epsilon = 1,
    pik_rhol_epsilon = 2,
    pi_rho_epsilonkl = 3,
    pik_rhol_epsilonkl = 4,
    pi_rho_unknown = 5,
    pik_rhol_unknown = 6,
    pi_rho_known = 7,
    pik_rhol_known = 8,
    pi_rho_sigma2 = 9,
    pik_rhol_sigma2 = 10,
    pi_rho_sigma2kl = 11,
    pik_rhol_sigma2kl = 12
};

#endif /* ENUMERATIONS_H_ */
