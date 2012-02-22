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
 * created on: Jan 10, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file Rcoclustermain.cpp
 *  @brief 
 **/

#include <Rcpp.h>
#include <time.h>
#include <exception>
#include <map>
//CoClust::library
#include "coclust/src/CoClustFacade/CoCluster.h"
//CoClust::strategy
#include "coclust/src/Strategy/XStrategyAlgo.h"
#include "coclust/src/Strategy/XStrategyforSEMAlgo.h"
//CoClust::Initialization
#include "coclust/src/Initialization/CEMInit.h"
#include "coclust/src/Initialization/FuzzyCEMInit.h"
#include "coclust/src/Initialization/RandomInit.h"
//CoClust::Models
#include "coclust/src/Models/BinaryLBModel.h"
#include "coclust/src/Models/BinaryLBModelequalepsilon.h"
#include "coclust/src/Models/ContingencyLBModel.h"
#include "coclust/src/Models/ContingencyLBModel_mu_i_nu_j.h"
#include "coclust/src/Models/ContinuousLBModel.h"
#include "coclust/src/Models/ContinuousLBModelequalsigma.h"
//CoClust::Algorithms
#include "coclust/src/Algorithms/EMAlgo.h"
#include "coclust/src/Algorithms/CEMAlgo.h"
#include "coclust/src/Algorithms/SEMAlgo.h"
//CoClust::Enumeration
#include "coclust/src/enumerations/enumerations.h"
//Package::DataExchangeFiles
#include "IDataExchange.h"
#include "ContinuousDataExchange.h"
#include "BinaryDataExchange.h"
#include "ContingencyDataExchange.h"



#ifdef SUPPORT_OPENMP
#include <omp.h>

RcppExport SEXP CoClustmain(SEXP robj)
{
  BEGIN_RCPP
  //Initialize Rcpp object
  Rcpp::S4 CoClustobj(robj);

  //set number of threads equal to number of processors
  omp_set_num_threads(omp_get_num_procs());

  std::map<std::string,DataType> S_DataType;
  S_DataType["binary"] = Binary;
  S_DataType["contingency"] = Contingency;
  S_DataType["continuous"] = Continuous;
  DataType datatype = S_DataType[Rcpp::as<std::string>(CoClustobj.slot("datatype"))];

  //Various pointers declarations
  IStrategy * p_Strategy_ = NULL;
  IAlgo * p_Algo_ = NULL;
  ICoClustModel * p_Model_ = NULL;
  ICoClustModel * p_FinalModel_ = NULL;
  IInit * p_Init_ = NULL;
  CoCluster* p_CoCluster_ = NULL;
  IDataExchange * p_DataExchange_;

  //set p_DataExchange_ to required case
  switch (datatype) {
    case Binary:
      p_DataExchange_ = new BinaryDataExchange();
      break;
    case Contingency:
      p_DataExchange_ = new ContingencyDataExchange();
      break;
    case Continuous:
      p_DataExchange_ = new ContinuousDataExchange();
      break;
    default:
      break;
  }

  //Get input parameters
  p_DataExchange_->SetInput(CoClustobj);
  //Get data
  p_DataExchange_->DataInput(CoClustobj);
  // Set Aparam_ for parallel use
  StrategyParameters Sparam_ = p_DataExchange_->GetStrategyParameters();
  int nbtry = Sparam_.nbtry_;
  Sparam_.nbtry_ = 1;

  //set Model
  ModelParameters Mparam_ = p_DataExchange_->GetModelParameters();
  if(!p_DataExchange_->GetStrategy().SemiSupervised){
    switch (datatype) {
      case Binary:{
        BinaryDataExchange * p_tempBinary_ = dynamic_cast<BinaryDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_epsilonkl:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new BinaryLBModel(p_tempBinary_->GetData(),Mparam_);
            break;
          case pik_rhol_epsilon:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),Mparam_);
            break;
          case pi_rho_epsilonkl:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new BinaryLBModel(p_tempBinary_->GetData(),Mparam_);
            break;
          case pi_rho_epsilon:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),Mparam_);
            break;
        }
      }
        break;
      case Contingency:{
        ContingencyDataExchange * p_tempContingency_ = dynamic_cast<ContingencyDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_unknown:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContingencyLBModel(p_tempContingency_->GetData(),Mparam_);
            break;
          case pik_rhol_known:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
          case pi_rho_unknown:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new ContingencyLBModel(p_tempContingency_->GetData(),Mparam_);
            break;
          case pi_rho_known:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
        }
      }
        break;
      case  Continuous:{
        ContinuousDataExchange * p_tempContinuous_ = dynamic_cast<ContinuousDataExchange *>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_sigma2kl:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContinuousLBModel(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pik_rhol_sigma2:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pi_rho_sigma2kl:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new ContinuousLBModel(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pi_rho_sigma2:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),Mparam_);
            break;
        }
      }
        break;
      default:
        break;
    }
  }else{
    switch (datatype) {
      case Binary:{
        BinaryDataExchange * p_tempBinary_ = dynamic_cast<BinaryDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_epsilonkl:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new BinaryLBModel(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                         p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_epsilon:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                                     p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_epsilonkl:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new BinaryLBModel(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                         p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_epsilon:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                                     p_DataExchange_->GetColLabels(),Mparam_);
            break;
        }
      }
        break;
      case Contingency:{
        ContingencyDataExchange * p_tempContingency_ = dynamic_cast<ContingencyDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_unknown:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContingencyLBModel(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                              p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_known:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                                        p_DataExchange_->GetColLabels(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
          case pi_rho_unknown:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new ContingencyLBModel(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                              p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_known:
            Mparam_.fixedproportions_ = true;
            p_FinalModel_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                                        p_DataExchange_->GetColLabels(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
        }
      }
        break;
      case  Continuous:{
        ContinuousDataExchange * p_tempContinuous_ = dynamic_cast<ContinuousDataExchange *>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_sigma2kl:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContinuousLBModel(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                             p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_sigma2:
            Mparam_.fixedproportions_ = false;
            p_FinalModel_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                                       p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_sigma2kl:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new ContinuousLBModel(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                             p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_sigma2:
            Mparam_.fixedproportions_ =true;
            p_FinalModel_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                                       p_DataExchange_->GetColLabels(),Mparam_);
            break;
        }
      }
        break;
      default:
        break;
    }
  }

  //use to update p_FinalModel_ in case the thread model is better then it.
  float Lmax = -RealMax;
  //to measure global success
  bool globalsuccess = false;
  //start measuring time after exchange of data and various initializations
  //double starttime=omp_get_wtime();
#pragma omp parallel shared(Lmax,globalsuccess,Sparam_,nbtry) private(p_Model_,p_CoCluster_,p_Algo_,p_Init_,p_Strategy_)
      {
        // set Initialization method
        switch (p_DataExchange_->GetStrategy().Init_) {
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

        // Set p_Algo_ and corresponding strategy p_Strategy_
        switch (p_DataExchange_->GetStrategy().Algo_)
        {
          case BEM:
            p_Algo_ = new EMAlgo();
            p_Strategy_ = new XStrategyAlgo(Sparam_);
            break;
          case BCEM:
            p_Algo_ = new CEMAlgo();
            p_Strategy_ = new XStrategyAlgo(Sparam_);
            break;
          case BSEM:
            p_Algo_ = new SEMAlgo();
            p_Strategy_ = new XStrategyforSEMAlgo(Sparam_);
            break;
          default:
            p_Algo_ = new EMAlgo();
            p_Strategy_ = new XStrategyAlgo(Sparam_);
            break;
        }

        // get copy of model for each thread
        p_Model_ = p_FinalModel_->Clone();

        //set cocluster
        p_CoCluster_ = new CoCluster();
        p_CoCluster_->SetStrategy(p_Strategy_);
        p_CoCluster_->SetAlgo(p_Algo_);
        p_CoCluster_->SetModel(p_Model_);
        p_CoCluster_->SetInit(p_Init_);

#pragma omp for schedule(dynamic,1)

      for (int i = 0; i < nbtry; ++i) {
        bool success = p_CoCluster_->run();
#pragma omp critical
        {
          if (Lmax<p_Model_->GetLikelihood()&& success) {
            Lmax = p_Model_->GetLikelihood();
            globalsuccess = true;
            delete p_FinalModel_;
            p_FinalModel_ = p_Model_->Clone();
          }
        }
      }

      p_Model_->ConsoleOut();
      delete p_Model_;
      delete p_Algo_;
      delete p_Init_;
      delete p_CoCluster_;

      }

  p_DataExchange_->Output(CoClustobj,p_FinalModel_,globalsuccess);
  //CoClustobj.slot("time") = omp_get_wtime()-starttime;

  //release memory
  delete p_DataExchange_;
  delete p_FinalModel_;

  return CoClustobj;
  END_RCPP
}
#else
RcppExport SEXP CoClustmain(SEXP robj)
{
  BEGIN_RCPP
  //Initialize Rcpp object
  Rcpp::S4 CoClustobj(robj);

  std::map<std::string,DataType> S_DataType;
  S_DataType["binary"] = Binary;
  S_DataType["contingency"] = Contingency;
  S_DataType["continuous"] = Continuous;

  DataType datatype = S_DataType[Rcpp::as<std::string>(CoClustobj.slot("datatype"))];

  //Various pointers declarations
  IStrategy * p_Strategy_ = NULL;
  IAlgo * p_Algo_ = NULL;
  ICoClustModel * p_Model_ = NULL;
  IInit * p_Init_ = NULL;
  CoCluster* p_CoCluster_ = NULL;
  IDataExchange * p_DataExchange_ = NULL;

  //set p_DataExchange_ to required case
  switch (datatype) {
    case Binary:
      p_DataExchange_ = new BinaryDataExchange();
      break;
    case Contingency:
      p_DataExchange_ = new ContingencyDataExchange();
      break;
    case Continuous:
      p_DataExchange_ = new ContinuousDataExchange();
      break;
    default:
      break;
  }

  //Get input parameters
  p_DataExchange_->SetInput(CoClustobj);
  //Get data
  p_DataExchange_->DataInput(CoClustobj);
  // set Initialization method
  switch (p_DataExchange_->GetStrategy().Init_) {
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
  // Set Strategy to be used to run algorithm method
  StrategyParameters Sparam_ = p_DataExchange_->GetStrategyParameters();

  // Set p_Algo_ and corresponding strategy p_Strategy_
  switch (p_DataExchange_->GetStrategy().Algo_)
  {
    case BEM:
      p_Algo_ = new EMAlgo();
      p_Strategy_ = new XStrategyAlgo(Sparam_);
      break;
    case BCEM:
      p_Algo_ = new CEMAlgo();
      p_Strategy_ = new XStrategyAlgo(Sparam_);
      break;
    case BSEM:
      p_Algo_ = new SEMAlgo();
      p_Strategy_ = new XStrategyforSEMAlgo(Sparam_);
      break;
    default:
      p_Algo_ = new EMAlgo();
      p_Strategy_ = new XStrategyAlgo(Sparam_);
      break;
  }

  //set Model
  ModelParameters Mparam_ = p_DataExchange_->GetModelParameters();
  if(!p_DataExchange_->GetStrategy().SemiSupervised){
    switch (datatype) {
      case Binary:{
        BinaryDataExchange * p_tempBinary_ = dynamic_cast<BinaryDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_epsilonkl:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new BinaryLBModel(p_tempBinary_->GetData(),Mparam_);
            break;
          case pik_rhol_epsilon:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),Mparam_);
            break;
          case pi_rho_epsilonkl:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new BinaryLBModel(p_tempBinary_->GetData(),Mparam_);
            break;
          case pi_rho_epsilon:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),Mparam_);
            break;
        }
      }
        break;
      case Contingency:{
        ContingencyDataExchange * p_tempContingency_ = dynamic_cast<ContingencyDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_unknown:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContingencyLBModel(p_tempContingency_->GetData(),Mparam_);
            break;
          case pik_rhol_known:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
          case pi_rho_unknown:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new ContingencyLBModel(p_tempContingency_->GetData(),Mparam_);
            break;
          case pi_rho_known:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
        }
      }
        break;
      case  Continuous:{
        ContinuousDataExchange * p_tempContinuous_ = dynamic_cast<ContinuousDataExchange *>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_sigma2kl:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContinuousLBModel(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pik_rhol_sigma2:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pi_rho_sigma2kl:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new ContinuousLBModel(p_tempContinuous_->GetData(),Mparam_);
            break;
          case pi_rho_sigma2:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),Mparam_);
            break;
        }
      }
        break;
      default:
        break;
    }
  }else{
    switch (datatype) {
      case Binary:{
        BinaryDataExchange * p_tempBinary_ = dynamic_cast<BinaryDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_epsilonkl:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new BinaryLBModel(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                         p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_epsilon:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                                     p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_epsilonkl:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new BinaryLBModel(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                         p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_epsilon:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new BinaryLBModelequalepsilon(p_tempBinary_->GetData(),p_DataExchange_->GetRowLabels(),
                                                     p_DataExchange_->GetColLabels(),Mparam_);
            break;
        }
      }
        break;
      case Contingency:{
        ContingencyDataExchange * p_tempContingency_ = dynamic_cast<ContingencyDataExchange*>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_unknown:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContingencyLBModel(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                              p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_known:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                                        p_DataExchange_->GetColLabels(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
          case pi_rho_unknown:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new ContingencyLBModel(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                              p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_known:
            Mparam_.fixedproportions_ = true;
            p_Model_ = new ContingencyLBModel_mu_i_nu_j(p_tempContingency_->GetData(),p_DataExchange_->GetRowLabels(),
                                                        p_DataExchange_->GetColLabels(),p_tempContingency_->GetMui()
                                                        ,p_tempContingency_->GetNuj(),Mparam_);
            break;
        }
      }
        break;
      case  Continuous:{
        ContinuousDataExchange * p_tempContinuous_ = dynamic_cast<ContinuousDataExchange *>(p_DataExchange_);
        switch (p_DataExchange_->GetStrategy().Model_)
        {
          case pik_rhol_sigma2kl:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContinuousLBModel(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                             p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pik_rhol_sigma2:
            Mparam_.fixedproportions_ = false;
            p_Model_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                                       p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_sigma2kl:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new ContinuousLBModel(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                             p_DataExchange_->GetColLabels(),Mparam_);
            break;
          case pi_rho_sigma2:
            Mparam_.fixedproportions_ =true;
            p_Model_ = new ContinuousLBModelequalsigma(p_tempContinuous_->GetData(),p_DataExchange_->GetRowLabels(),
                                                       p_DataExchange_->GetColLabels(),Mparam_);
            break;
        }
      }
        break;
      default:
        break;
    }
  }

  //create cocluster object and run coclustering
  double begin =clock();
  p_CoCluster_ = new CoCluster();
  p_CoCluster_->SetStrategy(p_Strategy_);
  p_CoCluster_->SetAlgo(p_Algo_);
  p_CoCluster_->SetModel(p_Model_);
  p_CoCluster_->SetInit(p_Init_);
  bool success = p_CoCluster_->run();
  p_DataExchange_->Output(CoClustobj,p_Model_,success);
  //CoClustobj.slot("time") = double(clock()-begin)/CLOCKS_PER_SEC;

  //release memory
  delete p_DataExchange_;
  delete p_Init_;
  delete p_Algo_;
  delete p_Model_;
  delete p_CoCluster_;

  return CoClustobj;
  END_RCPP
}
#endif
