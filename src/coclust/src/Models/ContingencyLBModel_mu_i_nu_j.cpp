/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2013  Parmeet Singh Bhatia

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

/** @file ContingencyLBModel_mu_i_nu_j.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j derived from ICoClustModel.
 **/

#include "ContingencyLBModel_mu_i_nu_j.h"
#include <iostream>

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j(MatrixInteger const& m_Dataij,VectorReal const& v_Mui,
                                                           VectorReal const& v_Nuj,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij),v_Mui_(v_Mui),v_Nuj_(v_Nuj)
{
  DataSum_ = m_Dataij.array().sum();
};

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j(MatrixInteger const& m_Dataij,VectorInteger const & rowlabels,VectorInteger const & collabels,
                                                           VectorReal const& v_Mui,VectorReal const& v_Nuj,ModelParameters const& Mparam)
                             : ICoClustModel(Mparam,rowlabels,collabels)
                             , m_Dataij_(m_Dataij),v_Mui_(v_Mui),v_Nuj_(v_Nuj)
{
  DataSum_ = m_Dataij.array().sum();
};

bool ContingencyLBModel_mu_i_nu_j::EMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ERows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepRows();
    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::SEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  if(!SERows()) return false;
  //M-step
  MStepRows();
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::CEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    if(!CERows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepRows();
    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::EMCols()
{
  //Initialization
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ECols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepCols();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::SEMCols()
{
  //Initialization
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;

  if(!SECols()) return false;
  //M-step
  MStepCols();

  return true;
}

bool ContingencyLBModel_mu_i_nu_j::CEMCols()
{
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!CECols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepCols();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

void ContingencyLBModel_mu_i_nu_j::LogSumRows(MatrixReal & m_ik)
{
  m_ik = VectorReal::Ones(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.array().log()).matrix().transpose())
      -v_Mui_*(m_Gammakl_*v_nul_).transpose();
}

void ContingencyLBModel_mu_i_nu_j::LogSumCols(MatrixReal & m_jl)
{
  m_jl = VectorReal::Ones(nbVar_)*v_logRhol_.transpose() + m_Vjk_*((m_Gammakl_.array().log()).matrix())
      -v_Nuj_*(v_muk_.transpose()*m_Gammakl_);
}

float ContingencyLBModel_mu_i_nu_j::EstimateLikelihood()
{
  Likelihood_ = (m_Ykl_.array()*(m_Gammakl_.array()).log()).sum() - DataSum_ + v_Tk_.transpose()*v_logPiek_ + v_Rl_.transpose()*v_logRhol_
            -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
            -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum();

  return Likelihood_;
}
void ContingencyLBModel_mu_i_nu_j::likelihoodStopCriteria()
{
  Likelihood_ = EstimateLikelihood();

  if(std::abs(1-Likelihood_/Likelihood_old)<Mparam_.epsilon_)
      StopAlgo = true;
  else
  {
    Likelihood_old = Likelihood_;
    StopAlgo = false;
  }
}

void ContingencyLBModel_mu_i_nu_j::ParameterStopCriteria()
{
  float relativechange = (((m_Gammakl1_.array()-m_Gammakl1old_.array()).abs()/m_Gammakl1_.array()).sum());
  if(relativechange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

bool ContingencyLBModel_mu_i_nu_j::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Running Initialization.."<<"\n";
#endif
  //Initialize random row and column partition
  VectorReal probarows = (1.0/Mparam_.nbrowclust_)*VectorReal::Ones(Mparam_.nbrowclust_);
  VectorReal probacols = (1.0/Mparam_.nbcolclust_)*VectorReal::Ones(Mparam_.nbcolclust_);
  v_Zi_ = PartRnd(nbSample_,probarows);
  v_Wj_ = PartRnd(nbVar_,probacols);
  m_Zik_ = MatrixInteger::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Wjl_ = MatrixInteger::Zero(nbVar_,Mparam_.nbcolclust_);
  m_Tik_.setZero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_.setZero(nbVar_,Mparam_.nbcolclust_);
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for( Label_pair : knownLabelsRows_){
    m_Tik_(Label_pair.first,Label_pair.second) = 1.0;
  }

  for ( int i : UnknownLabelsRows_) {
    m_Tik_(i,v_Zi_(i)-1) = 1.0;
  }

  for (Label_pair : knownLabelsCols_) {
    m_Rjl_(Label_pair.first,Label_pair.second)=1.0;
  }

  for ( int j : UnknownLabelsCols_) {
    m_Rjl_(j,v_Wj_(j)-1) = 1.0;
  }
#else
  for(int i =0;i<knownLabelsRows_.size();i++){
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1.0;
  }

  for ( int i =0;i< UnknownLabelsRows_.size();i++) {
    m_Tik_(i,v_Zi_(UnknownLabelsRows_[i])-1) = 1.0;
  }

  for ( int j=0;j<knownLabelsCols_.size();j++) {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1.0;
  }

  for ( int j =0;j< UnknownLabelsCols_.size();j++) {
    m_Rjl_(j,v_Wj_(UnknownLabelsCols_[j])-1) = 1.0;
  }
#endif

  v_Rl_ = m_Rjl_.colwise().sum().transpose();
  //Check for empty cluster
  if((v_Rl_.array()<.00001).any()){
    Error_msg_  = "Column clustering failed while running initialization.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  empty_cluster = true;
  return false;
  }
  else
  {empty_cluster = false;}
  m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_).array()/(m_Tik_.transpose()*v_Mui_*v_Nuj_.transpose()*m_Rjl_).array();

  //Initializing model parameters
  m_Gammakl1_ = m_Gammakl_;
  m_Gammakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  m_Gammaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  m_Uil_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
  m_Vjk_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Tk_ = VectorReal::Zero(Mparam_.nbrowclust_);
  v_Zi_ = VectorInteger::Zero(nbSample_);
  v_Wj_ = VectorInteger::Zero(nbVar_);
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(VectorReal::Ones(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(VectorReal::Ones(Mparam_.nbcolclust_));

#ifdef COVERBOSE
  std::cout<<"Initialization over."<<"\n";
#endif

  return true;
}

VectorInteger ContingencyLBModel_mu_i_nu_j::PartRnd(int n,VectorReal proba)
{
  int clusters = proba.rows();
  VectorInteger v_Z = VectorInteger::Zero(n);
  VectorInteger v_Randperm = RandSample(n,n);
  VectorInteger remainingclusters(n-clusters);

  for ( int ind = 0; ind < clusters; ++ind) {
    v_Z(v_Randperm(ind)) = ind+1;
  }
  remainingclusters = (clusters+1)*MatrixInteger::Ones(n-clusters,1) - ((MatrixReal::Ones(clusters,1)*(Unifrnd(0,1,1,n-clusters))).array() < (Cumsum(proba)*MatrixReal::Ones(1,n-clusters)).array()).matrix().cast<int>().colwise().sum().transpose();

  for ( int ind = clusters; ind < n; ++ind) {
    v_Z(v_Randperm(ind)) = remainingclusters(ind-clusters)==(clusters+1)?clusters:remainingclusters(ind-clusters);
  }
  return v_Z;

}

VectorReal ContingencyLBModel_mu_i_nu_j::Cumsum(VectorReal proba)
{
  int size = proba.rows();
  VectorReal v_temp = VectorReal::Zero(size);
  v_temp(0) = proba(0);
  for ( int itr = 1; itr < size; ++itr) {
    v_temp(itr) = v_temp(itr-1) + proba(itr);
  }
  return v_temp;
}

MatrixReal ContingencyLBModel_mu_i_nu_j::Unifrnd(float a,float b, int row, int col)
{
  MatrixReal m_temp(row,col);
  for ( int r = 0; r < row; ++r) {
    for ( int c = 0; c < col; ++c) {
      m_temp(r,c) = (b-a)*(std::rand()/float(RAND_MAX)) + a;
    }
  }
  return m_temp;
}

void ContingencyLBModel_mu_i_nu_j::FinalizeOutput()
{
  CommonFinalizeOutput();
}

void ContingencyLBModel_mu_i_nu_j::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

const MatrixInteger& ContingencyLBModel_mu_i_nu_j::GetArrangedDataClusters()
{
  Eigen::ArrayXi v_Zi = (GetRowClassificationVector()).array();
  Eigen::ArrayXi v_Wj = (GetColumnClassificationVector()).array();
  m_ClusterDataij_ = MatrixInteger::Zero(nbSample_,nbVar_);

  //Rearrange data into clusters

  VectorInteger rowincrement = MatrixInteger::Zero(Mparam_.nbrowclust_,1);

  VectorInteger nbindrows = MatrixInteger::Zero(Mparam_.nbrowclust_,1);
  for ( int k = 1; k < Mparam_.nbrowclust_; ++k) {
    nbindrows(k) = (v_Zi==(k-1)).count()+nbindrows(k-1);
  }

  VectorInteger colincrement = MatrixInteger::Zero(Mparam_.nbcolclust_,1);

  VectorInteger nbindcols = MatrixInteger::Zero(Mparam_.nbcolclust_,1);
  for ( int l = 1; l < Mparam_.nbcolclust_; ++l) {
    nbindcols(l)=(v_Wj==(l-1)).count()+nbindcols(l-1);
  }

  for ( int j = 0; j < nbVar_; ++j) {
    m_ClusterDataij_.col(colincrement(v_Wj(j)) + nbindcols(v_Wj(j))) = m_Dataij_.col(j);
    colincrement(v_Wj(j))+=1;
  }
  MatrixInteger temp = m_ClusterDataij_;

  for ( int i = 0; i < nbSample_; ++i) {
    m_ClusterDataij_.row( rowincrement(v_Zi(i)) + nbindrows(v_Zi(i))) = temp.row(i);
    rowincrement(v_Zi(i))+=1;
  }

  return m_ClusterDataij_;

}

void ContingencyLBModel_mu_i_nu_j::Modify_theta_start()
{
  m_Gammaklstart_ = m_Gammakl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

void ContingencyLBModel_mu_i_nu_j::Copy_theta_start()
{
  m_Gammakl_ = m_Gammaklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
  m_Gammakl1_ = m_Gammakl_;
}

void ContingencyLBModel_mu_i_nu_j::Copy_theta_max()
{
  m_Gammakl_ = m_Gammaklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

void ContingencyLBModel_mu_i_nu_j::Modify_theta_max()
{
  m_Gammaklmax_ = m_Gammakl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}

void ContingencyLBModel_mu_i_nu_j::MStepFull()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_;
  m_Gammakl_ = m_Ykl_.array()/(m_Tik_.transpose()*v_Mui_*v_Nuj_.transpose()*m_Rjl_).array();
}
