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

#include "IDataExchange.h"

IDataExchange::~IDataExchange()
{}

void IDataExchange::SetInput(Rcpp::S4 & obj)
{
  //Get Strategy
  Rcpp::S4 strategy(obj.slot("strategy"));
  strategy_.Algo_ = S_Algorithm[Rcpp::as<std::string>(strategy.slot("algo"))];
  strategy_.stopcriteria_ = S_StopCriteria[Rcpp::as<std::string>(strategy.slot("stopcriteria"))];
  strategy_.Init_ = S_Init[Rcpp::as<std::string>(strategy.slot("initmethod"))];
  strategy_.DataType_ = S_DataType[Rcpp::as<std::string>(obj.slot("datatype"))];
  strategy_.Model_ = S_Model[Rcpp::as<std::string>(obj.slot("model"))];
  strategy_.SemiSupervised = Rcpp::as<bool>(obj.slot("semisupervised"));

  //Get strategy parameters
  Stratparam_.nbiter_xem_ = Rcpp::as<int>(strategy.slot("nbiterationsxem"));
  Stratparam_.nbiter_XEM_ = Rcpp::as<int>(strategy.slot("nbiterationsXEM"));
  Stratparam_.nbtry_ = Rcpp::as<int>(strategy.slot("nbtry"));
  Stratparam_.nbxem_ = Rcpp::as<int>(strategy.slot("nbxem"));

  //get row/column labels if semisupervised
  if(strategy_.SemiSupervised)
  {
    v_rowlabels_ = convertvector<VectorInteger,Rcpp::NumericVector>((SEXP(obj.slot("rowlabels"))));
    v_collabels_ = convertvector<VectorInteger,Rcpp::NumericVector>((SEXP(obj.slot("collabels"))));
  }

  // Set stopping-criteria
  switch (strategy_.stopcriteria_)
  {
    case Parameter:
      Stratparam_.Stop_Criteria = &ICoClustModel::ParameterStopCriteria;
      break;
    case Likelihood:
      Stratparam_.Stop_Criteria = &ICoClustModel::likelihoodStopCriteria;
      break;
    default:
      Stratparam_.Stop_Criteria = &ICoClustModel::ParameterStopCriteria;
      break;
  }
  //Set various  model parameters

  //get various iterations
  Mparam_.nbinititerations_ = Rcpp::as<int>(strategy.slot("nbinititerations"));
  Mparam_.nbiterations_int_ = Rcpp::as<int>(strategy.slot("nbiterations_int"));

  //get threshold values
  Mparam_.eps_xem_ = Rcpp::as<double>(strategy.slot("epsilonxem"));
  Mparam_.eps_XEM_ = Rcpp::as<double>(strategy.slot("epsilonXEM"));
  Mparam_.initepsilon_ = Rcpp::as<double>(strategy.slot("initepsilon"));
  Mparam_.epsilon_int_ = Rcpp::as<double>(strategy.slot("epsilon_int"));

  //get row and column clusters
  Rcpp::IntegerVector nbcocluster(SEXP(obj.slot("nbcocluster")));
  Mparam_.nbrowclust_ = nbcocluster(0);
  Mparam_.nbcolclust_ = nbcocluster(1);
}

void IDataExchange::InitializeParamEnum()
{
  //Datatype
  S_DataType["binary"] = Binary;
  S_DataType["contingency"] = Contingency;
  S_DataType["continuous"] = Continuous;

  //Algorithm
  S_Algorithm["BEM"] = BEM;
  S_Algorithm["BCEM"] = BCEM;
  S_Algorithm["BSEM"] = BSEM;

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
