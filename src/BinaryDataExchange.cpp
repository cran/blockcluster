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

/** @file BinaryCoClust.cpp
 *  @brief 
 **/

#include "BinaryDataExchange.h"
#include "coclust/src/Models/BinaryLBModel.h"
#include "coclust/src/Models/BinaryLBModelequalepsilon.h"

void BinaryDataExchange::Output(Rcpp::S4& obj,ICoClustModel* model,bool successful)
{
  if(!successful)
  {
    obj.slot("successful") = false;
    obj.slot("message") = model->GetErrormsg();
  }
  else
  {
    obj.slot("successful") = true;
    obj.slot("message") = "Co-Clustering successfully terminated!";
    BinaryLBModel* ptrLBM;
    BinaryLBModelequalepsilon* ptrLBMeq;
    MatrixReal dispersion;
    switch (strategy_.Model_)
    {
      case pik_rhol_epsilonkl:
         ptrLBM = dynamic_cast<BinaryLBModel*>(model);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->Getdispersion());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->GetArrangedDataClusters());
        break;
      case pik_rhol_epsilon:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(model);
        dispersion = ptrLBMeq->Getdispersion()*MatrixReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(dispersion);
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->GetArrangedDataClusters());
        break;
      case pi_rho_epsilonkl:
        ptrLBM = dynamic_cast<BinaryLBModel*>(model);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->Getdispersion());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->GetArrangedDataClusters());
        break;
      case pi_rho_epsilon:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(model);
        dispersion = ptrLBMeq->Getdispersion()*MatrixReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(dispersion);
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->GetArrangedDataClusters());
        break;
    }
    obj.slot("rowclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetRowClassificationVector());
    obj.slot("colclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetColumnClassificationVector());
    obj.slot("rowproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetRowProportions());
    obj.slot("columnproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetColProportions());
    obj.slot("rowposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetRowPosteriorprob());
    obj.slot("colposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetColPosteriorprob());
    obj.slot("likelihood") = model->GetLikelihood();
    obj.slot("ICLvalue") = model->ICLCriteriaValue();
  }
}

void BinaryDataExchange::DataInput(Rcpp::S4 & obj)
{
  Rcpp::NumericMatrix data(SEXP(obj.slot("data")));
  convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(data,m_Dataij_);
  Mparam_.nbrowdata_ = m_Dataij_.rows();
  Mparam_.nbcoldata_ = m_Dataij_.cols();

  //Get Strategy
  Rcpp::S4 strategy(obj.slot("strategy"));
  //get hyper-parameters
  if(Rcpp::as<bool>(strategy.slot("bayesianform")))
  {
    Rcpp::IntegerVector hyperparam(SEXP(strategy.slot("hyperparam")));
    a_ = hyperparam(0);
    b_ = hyperparam(1);
  }
}

void BinaryDataExchange::instantiateModel(ICoClustModel*& model){
  if(!strategy_.Bayesianform_){
    if(!strategy_.SemiSupervised){
      switch (strategy_.Model_)
      {
        case pik_rhol_epsilonkl:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModel(m_Dataij_,Mparam_);
          break;
        case pik_rhol_epsilon:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_);
          break;
        case pi_rho_epsilonkl:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModel(m_Dataij_,Mparam_);
          break;
        case pi_rho_epsilon:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_);
          break;
        default:
          break;
      }
    }
    else{
      switch (strategy_.Model_)
      {
        case pik_rhol_epsilonkl:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModel(m_Dataij_,v_rowlabels_,
                                       v_collabels_,Mparam_);
          break;
        case pik_rhol_epsilon:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_,
                                                   v_collabels_,Mparam_);
          break;
        case pi_rho_epsilonkl:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModel(m_Dataij_,v_rowlabels_,
                                       v_collabels_,Mparam_);
          break;
        case pi_rho_epsilon:
          Mparam_.fixedproportions_ =true;
          model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_,
                                                   v_collabels_,Mparam_);
          break;
        default:
          break;
      }
    }
  }else{
    if(!strategy_.SemiSupervised){
      switch (strategy_.Model_)
      {
        case pik_rhol_epsilonkl:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModel(m_Dataij_,Mparam_,a_,b_);
          break;
        case pik_rhol_epsilon:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_,a_,b_);
          break;
        case pi_rho_epsilonkl:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModel(m_Dataij_,Mparam_,a_,b_);
          break;
        case pi_rho_epsilon:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_,a_,b_);
          break;
        default:
          break;
      }
    }
    else{
      switch (strategy_.Model_)
      {
        case pik_rhol_epsilonkl:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModel(m_Dataij_,v_rowlabels_,
                                       v_collabels_,Mparam_,a_,b_);
          break;
        case pik_rhol_epsilon:
          Mparam_.fixedproportions_ = false;
          model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_,
                                                   v_collabels_,Mparam_,a_,b_);
          break;
        case pi_rho_epsilonkl:
          Mparam_.fixedproportions_ = true;
          model = new BinaryLBModel(m_Dataij_,v_rowlabels_,
                                       v_collabels_,Mparam_,a_,b_);
          break;
        case pi_rho_epsilon:
          Mparam_.fixedproportions_ =true;
          model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_,
                                                   v_collabels_,Mparam_,a_,b_);
          break;
        default:
          break;
      }
    }
  }
}
