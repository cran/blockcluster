/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

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
 * Project:  CoCluster-Debug@Build
 * created on: Jan 6, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file CEMAlgo.h
 *  @brief Declares CEM algorithm class CEMAlgo derived from IAlgo.
 **/

#ifndef CEMALGO_H_
#define CEMALGO_H_

#include "IAlgo.h"

/** @brief Concrete Algorithm class derived from IAlgo and declares CEM algorithm. It contains only one function
 * that is required to run CEM algorithm.
 */
class CEMAlgo : public IAlgo
{
  public:
    /**Constructor*/
    CEMAlgo(AlgoParameters const& Aparam)
    :IAlgo(Aparam)
    {};
    /**This function runs the CEM algorithm for model pointed by IAlgo::p_Model_.*/
    virtual bool run();
    /**Destructor*/
    virtual ~CEMAlgo(){};
};

#endif /* CEMALGO_H_ */
