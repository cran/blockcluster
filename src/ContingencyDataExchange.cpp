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

void ContingencyDataExchange::Output(Rcpp::S4& obj,ICoClustModel* p_Model_,bool successful)
{
  if(!successful)
  {
    obj.slot("successful") = false;
    obj.slot("message") = p_Model_->GetErrormsg();
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
         ptrLBM = dynamic_cast<ContingencyLBModel*>(p_Model_);
         obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->GetGamma());
         obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBM->GetArrangedDataClusters());
        break;
      case pik_rhol_known:
        ptrLBMknown = dynamic_cast<ContingencyLBModel_mu_i_nu_j*>(p_Model_);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBMknown->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBMknown->GetArrangedDataClusters());
        break;
      case pi_rho_unknown:
        ptrLBM = dynamic_cast<ContingencyLBModel*>(p_Model_);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBM->GetArrangedDataClusters());
        break;
      case pi_rho_known:
        ptrLBMknown = dynamic_cast<ContingencyLBModel_mu_i_nu_j*>(p_Model_);
        obj.slot("classgamma") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBMknown->GetGamma());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixInteger>(ptrLBMknown->GetArrangedDataClusters());
        break;
    }
    obj.slot("rowclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(p_Model_->GetRowClassificationVector());
    obj.slot("colclass") = convertvector<Rcpp::IntegerVector,VectorInteger>(p_Model_->GetColumnClassificationVector());
    obj.slot("rowproportions") = convertvector<Rcpp::NumericVector,VectorReal>(p_Model_->GetRowProportions());
    obj.slot("columnproportions") = convertvector<Rcpp::NumericVector,VectorReal>(p_Model_->GetColProportions());
    obj.slot("rowposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(p_Model_->GetRowPosteriorprob());
    obj.slot("colposteriorprob") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(p_Model_->GetColPosteriorprob());
    obj.slot("likelihood") = p_Model_->GetLikelihood();
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
