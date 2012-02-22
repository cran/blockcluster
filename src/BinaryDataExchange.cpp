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

void BinaryDataExchange::Output(Rcpp::S4& obj,ICoClustModel* p_Model_,bool successful)
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
    BinaryLBModel* ptrLBM;
    BinaryLBModelequalepsilon* ptrLBMeq;
    MatrixReal dispersion;
    switch (strategy_.Model_)
    {
      case pik_rhol_epsilonkl:
         ptrLBM = dynamic_cast<BinaryLBModel*>(p_Model_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->Getdispersion());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->GetArrangedDataClusters());
        break;
      case pik_rhol_epsilon:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(p_Model_);
        dispersion = ptrLBMeq->Getdispersion()*MatrixReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(dispersion);
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->GetArrangedDataClusters());
        break;
      case pi_rho_epsilonkl:
        ptrLBM = dynamic_cast<BinaryLBModel*>(p_Model_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(ptrLBM->Getdispersion());
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBM->GetArrangedDataClusters());
        break;
      case pi_rho_epsilon:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(p_Model_);
        dispersion = ptrLBMeq->Getdispersion()*MatrixReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        obj.slot("classmean") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->Getmean());
        obj.slot("classdispersion") = convertMatrix<Rcpp::NumericMatrix,MatrixReal>(dispersion);
        obj.slot("coclusterdata") = convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(ptrLBMeq->GetArrangedDataClusters());
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

void BinaryDataExchange::DataInput(Rcpp::S4 & obj)
{
  Rcpp::NumericMatrix data(SEXP(obj.slot("data")));
  convertMatrix<Rcpp::NumericMatrix,MatrixBinary>(data,m_Dataij_);
  Mparam_.nbrowdata_ = m_Dataij_.rows();
  Mparam_.nbcoldata_ = m_Dataij_.cols();
}
