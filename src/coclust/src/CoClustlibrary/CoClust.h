
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
 * Project:  cocluster
 * created on: Feb 21, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file CoClust.h
 *  @brief 
 **/

#ifndef RPACKAGE
#ifndef COCLUST_H_
#define COCLUST_H_

#include "CoCluster.h"
#include "../Algorithms/IAlgo.h"
#include "../DataHandling/BinaryData.h"
#include "../DataHandling/ContinuousData.h"
#include "../DataHandling/ContingencyData.h"
#include "../InputParameters/InputParameters.h"
#include "../InputParameters/BinaryInputParameters.h"
#include "../InputParameters/ContinuousInputParameters.h"
#include "../InputParameters/ContingencyInputParameters.h"
#include "../Models/ICoClustModel.h"
#include "../StoppingCriteria/IStopCriteria.h"
#include "../Initialization/IInit.h"
#include "../Algorithms/EMAlgo.h"
#include "../Algorithms/CEMAlgo.h"
#include "../Models/BinaryLBModel.h"
#include "../Models/ContinuousLBModel.h"
#include "../Models/ContingencyLBModel.h"
#include "../Models/ContingencyLBModel_mu_i_nu_j.h"
#include "../StoppingCriteria/ParameterCriteria.h"
#include "../StoppingCriteria/LikelihoodCriteria.h"
#include "../Initialization/CEMInit.h"
#include "../Initialization/FuzzyCEMInit.h"
#include "../Initialization/RandomInit.h"
#include "../../../stkpp/include/STKernel.h"

#endif /* COCLUST_H_ */
#endif /* RPACKAGE*/
