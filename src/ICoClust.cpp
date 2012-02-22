/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

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

    Contact : bhatia.parmeet@gmail.com , parmeet.bhatia@inria.fr
*/

/*
 * Project:  Rcoclust::
 * created on: Apr 11, 2012
 * Purpose:  .
 * Author:   modal-laptop, parmeet.bhatia@inria.fr
 *
 **/

/** @file ICoClust.cpp
 *  @brief In this file .
 **/

#include "ICoClust.h"

ICoClust::~ICoClust()
{
  if(p_Algo_)delete p_Algo_;
  if(p_StopCriteria_)delete p_StopCriteria_;
  if(p_Init_)delete p_Init_;
  if(p_Model_)delete p_Model_;
  if(p_CoClust_)delete p_CoClust_;
}

bool ICoClust::Run()
{
  p_CoClust_ = new CoCluster();
  p_CoClust_->SetAlgo(p_Algo_);
  p_CoClust_->SetModel(p_Model_);
  p_CoClust_->SetInit(p_Init_);
  p_CoClust_->SetStopCriteria(p_StopCriteria_);
  return p_CoClust_->run();
}

void ICoClust::SetStrategy(Rcpp::S4 & obj,Rcpp::S4 & strategy)
{
  //Get Strategy
  strategy_.Algo_ = S_Algorithm[Rcpp::as<std::string>(strategy.slot("algo"))];
  strategy_.stopcriteria_ = S_StopCriteria[Rcpp::as<std::string>(strategy.slot("stopcriteria"))];
  strategy_.Init_ = S_Init[Rcpp::as<std::string>(strategy.slot("initmethod"))];
  strategy_.DataType_ = S_DataType[Rcpp::as<std::string>(obj.slot("datatype"))];
  strategy_.Model_ = S_Model[Rcpp::as<std::string>(obj.slot("model"))];
  //Get Algo parameters
  Aparam_.nbiter_xem_ = Rcpp::as<int>(strategy.slot("nbiterationsxem"));
  Aparam_.nbiter_XEM_ = Rcpp::as<int>(strategy.slot("nbiterationsXEM"));
  Aparam_.nbtry_ = Rcpp::as<int>(strategy.slot("nbtry"));
  Aparam_.nbxem_ = Rcpp::as<int>(strategy.slot("nbxem"));

  //set strategy
  switch (strategy_.stopcriteria_)
  {
    case Parameter:
      p_StopCriteria_ = new ParameterCriteria();
      break;
    case Likelihood:
      p_StopCriteria_ = new LikelihoodCriteria();
    default:
      p_StopCriteria_ = new ParameterCriteria();
      break;
  }
  // Set Algorithm
  switch (strategy_.Algo_)
  {
    case BCEM:
      p_Algo_ = new CEMAlgo(Aparam_);
      break;
    case BEM2:
      p_Algo_ = new EMAlgo(Aparam_);
      break;
    case XEMStrategy:
      p_Algo_ = new XEMStrategyAlgo(Aparam_);
      break;
    case XCEMStrategy:
      p_Algo_ = new XCEMStrategyAlgo(Aparam_);
      break;
    default:
      p_Algo_ = new XEMStrategyAlgo(Aparam_);
      break;
  }

  // set Initialization method
  switch (strategy_.Init_) {
    case e_CEMInit:
      p_Init_ = new CEMInit();
      break;
    case e_FuzzyCEMInit:
      p_Init_ = new FuzzyCEMInit();
      break;
    case e_RandomInit:
      p_Init_ = new RandomInit();
      break;
    default:
      p_Init_ = new CEMInit();
      break;
  }
}

void ICoClust::InitializeEnum()
{
  //Datatype
  S_DataType["binary"] = Binary;
  S_DataType["contingency"] = Contingency;
  S_DataType["continuous"] = Continuous;

  //Algorithm
  S_Algorithm["XEMStrategy"] = XEMStrategy;
  S_Algorithm["BEM2"] = BEM2;
  S_Algorithm["XCEMStrategy"] = XCEMStrategy;
  S_Algorithm["BCEM"] = BCEM;

  //StopCriteria
  S_StopCriteria["Parameter"] = Parameter;
  S_StopCriteria["Likelihood"] = Likelihood;

  //Initialization
  S_Init["CEMInit"] = e_CEMInit;
  S_Init["FuzzyCEMInit"] = e_FuzzyCEMInit;
  S_Init["RandomInit"] = e_RandomInit;

  //Models
  S_Model["pi_rho_epsilon"] = pi_rho_epsilon;
  S_Model["pik_rhol_epsilon"] = pik_rhol_epsilon;
  S_Model["pi_rho_epsilonkl"] = pi_rho_epsilonkl;
  S_Model["pik_rhol_epsilonkl"] = pik_rhol_epsilonkl;
  S_Model["pi_rho_unknown"] = pi_rho_unknown;
  S_Model["pik_rhol_unknown"] = pik_rhol_unknown;
  S_Model["pi_rho_known"] = pi_rho_known;
  S_Model["pik_rhol_known"] = pik_rhol_known;
  S_Model["pi_rho_sigma2"] = pi_rho_sigma2;
  S_Model["pik_rhol_sigma2"] = pik_rhol_sigma2;
  S_Model["pi_rho_sigma2kl"] = pi_rho_sigma2kl;
  S_Model["pik_rhol_sigma2kl"] = pik_rhol_sigma2kl;
}
