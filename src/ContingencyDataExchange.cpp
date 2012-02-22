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

/** @file ContingencyCoClust.cpp
 *  @brief 
 **/

#include "ContingencyDataExchange.h"
#include "coclust/src/Models/ContingencyLBModel.h"
#include "coclust/src/Models/ContingencyLBModel_mu_i_nu_j.h"
void ContingencyDataExchange::Output(Rcpp::S4& obj,ICoClustModel* model,bool successful)
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
    ContingencyLBModel* ptrLBM;
    ContingencyLBModel_mu_i_nu_j* ptrLBMknown;
    switch (strategy_.Model_)
    {
      case pik_rhol_unknown:
         ptrLBM = dynamic_cast<ContingencyLBModel*>(model);
         obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->GetGamma());
         obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBM->GetArrangedDataClusters());
        break;
      case pik_rhol_known:
        ptrLBMknown = dynamic_cast<ContingencyLBModel_mu_i_nu_j*>(model);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBMknown->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBMknown->GetArrangedDataClusters());
        break;
      case pi_rho_unknown:
        ptrLBM = dynamic_cast<ContingencyLBModel*>(model);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBM->GetArrangedDataClusters());
        break;
      case pi_rho_known:
        ptrLBMknown = dynamic_cast<ContingencyLBModel_mu_i_nu_j*>(model);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBMknown->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBMknown->GetArrangedDataClusters());
        break;
    }
    obj.slot("rowclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetRowClassificationVector());
    obj.slot("colclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(model->GetColumnClassificationVector());
    obj.slot("rowproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetRowProportions());
    obj.slot("columnproportions") = convertvector<Rcpp::NumericVector,VectorReal>(model->GetColProportions());
    obj.slot("rowposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetRowPosteriorprob());
    obj.slot("colposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(model->GetColPosteriorprob());
    obj.slot("likelihood") = model->GetLikelihood();
  }
}

void ContingencyDataExchange::DataInput(Rcpp::S4 & obj)
{
  Rcpp::NumericMatrix data(SEXP(obj.slot("data")));
  convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(data,m_Dataij_);
  Mparam_.nbrowdata_ = m_Dataij_.rows();
  Mparam_.nbcoldata_ = m_Dataij_.cols();
  if(strategy_.Model_ == pik_rhol_known||strategy_.Model_ == pi_rho_known )
  {
    Rcpp::NumericVector datamui(SEXP(obj.slot("datamui")));
    v_Mui_ = convertvector<VectorReal,Rcpp::NumericVector>(datamui);
    Rcpp::NumericVector datanuj(SEXP(obj.slot("datanuj")));
    v_Nuj_ = convertvector<VectorReal,Rcpp::NumericVector>(datanuj);
  }
}

void ContingencyDataExchange::instantiateModel(ICoClustModel*& model){
  if(!strategy_.SemiSupervised){
    switch (strategy_.Model_)
    {
      case pik_rhol_unknown:
        Mparam_.fixedproportions_ = false;
        model = new ContingencyLBModel(m_Dataij_,Mparam_);
        break;
      case pik_rhol_known:
        Mparam_.fixedproportions_ = false;
        model = new ContingencyLBModel_mu_i_nu_j(m_Dataij_,v_Mui_
                                                    ,v_Nuj_,Mparam_);
        break;
      case pi_rho_unknown:
        Mparam_.fixedproportions_ = true;
        model = new ContingencyLBModel(m_Dataij_,Mparam_);
        break;
      case pi_rho_known:
        Mparam_.fixedproportions_ = true;
        model = new ContingencyLBModel_mu_i_nu_j(m_Dataij_,v_Mui_
                                                    ,v_Nuj_,Mparam_);
        break;
      default:
        break;
    }
  }
  else{
    switch (strategy_.Model_)
    {
      case pik_rhol_unknown:
        Mparam_.fixedproportions_ = false;
        model = new ContingencyLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_);
        break;
      case pik_rhol_known:
        Mparam_.fixedproportions_ = false;
        model = new ContingencyLBModel_mu_i_nu_j(m_Dataij_,v_rowlabels_,v_collabels_,v_Mui_
                                                    ,v_Nuj_,Mparam_);
        break;
      case pi_rho_unknown:
        Mparam_.fixedproportions_ = true;
        model = new ContingencyLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_);
        break;
      case pi_rho_known:
        Mparam_.fixedproportions_ = true;
        model = new ContingencyLBModel_mu_i_nu_j(m_Dataij_,v_rowlabels_,v_collabels_,v_Mui_
                                                    ,v_Nuj_,Mparam_);
        break;
      default:
        break;
    }
  }
}
