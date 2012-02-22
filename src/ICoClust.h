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
 * Project:  RCocluster
 * created on: Feb 22, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ICoClust.h
 *  @brief 
 **/

#ifndef ICOCLUST_H_
#define ICOCLUST_H_
#include <iostream>
#include "coclust/src/CoClustlibrary/CoCluster.h"
#include "coclust/src/Initialization/CEMInit.h"
#include "coclust/src/Initialization/FuzzyCEMInit.h"
#include "coclust/src/Initialization/RandomInit.h"
#include "coclust/src/Algorithms/EMAlgo.h"
#include "coclust/src/Algorithms/CEMAlgo.h"
#include "coclust/src/Algorithms/XEMStrategyAlgo.h"
#include "coclust/src/Algorithms/XCEMStrategyAlgo.h"
#include "coclust/src/StoppingCriteria/ParameterCriteria.h"
#include "coclust/src/StoppingCriteria/LikelihoodCriteria.h"
#include <Rcpp.h>

class IAlgo;
class ICoClustModel;
class IStopCriteria;
class IInit;
/** @brief
 *
 */
class ICoClust
{
  public:
    ICoClust(){};
    virtual void ModelInput(Rcpp::S4 & obj,Rcpp::S4 & strategy) = 0;
    virtual void SetModel() = 0;
    virtual void Output(Rcpp::S4& obj,bool) = 0;
    void InitializeEnum();
    void SetStrategy(Rcpp::S4 & obj,Rcpp::S4 & strategy);
    bool Run();
    virtual ~ICoClust();
  protected:
    Strategy strategy_;
    AlgoParameters Aparam_;
    ModelParameters Mparam_;
    std::map<std::string,DataType> S_DataType;
    std::map<std::string,Algorithm> S_Algorithm;
    std::map<std::string,StopCriteria> S_StopCriteria;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model> S_Model;
    IAlgo * p_Algo_;
    ICoClustModel * p_Model_;
    IInit * p_Init_;
    IStopCriteria * p_StopCriteria_;
    CoCluster * p_CoClust_;
};





#endif /* ICOCLUST_H_ */
