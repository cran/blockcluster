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
 * created on: Jan 30, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

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
  nbSample_ = m_Dataij.rows();
  nbVar_ = m_Dataij.cols();
  Mparam_.nbrowdata_ = nbSample_;
  Mparam_.nbcoldata_ = nbVar_;
  DataSum_ = m_Dataij.array().sum();
};

bool ContingencyLBModel_mu_i_nu_j::Estep()
{
  if (EMRows()) {
    if (EMCols()) {
      return true;
    }
  }
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::CEstep()
{
  if (CEMRows()) {
    if (CEMCols()) {
      return true;
    }
  }
  return false;
}

void ContingencyLBModel_mu_i_nu_j::Mstep()
{
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
}


bool ContingencyLBModel_mu_i_nu_j::EMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  //Temporary variables
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);

  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      PoissonLogsumRows(m_sumik);
      Zsumk.setZero();
      for (int i = 0; i < nbSample_; ++i) {
        v_sumikmax(i) = m_sumik.row(i).maxCoeff(&maxIndex);
        Zsumk(maxIndex)+=1;
      }

      if((Zsumk.array()<0.00001).any())
      {
        empty_cluster = true;
        throw "Row clustering failed while running model.";
      }
      else
      {
        empty_cluster = false;
      }

      m_prodik=(m_sumik-v_sumikmax*Onesk.transpose()).array().exp();
      v_sumi = m_prodik.rowwise().sum();
      m_Tik_ = m_prodik.array()/(v_sumi*Onesk.transpose()).array();
      v_Tk_ = m_Tik_.colwise().sum();

      //M-step
      if(!Mparam_.fixedproportions_) {
        v_logPiek_=(v_Tk_.array()/nbSample_).log();
      }

      m_Ykl_ = m_Tik_.transpose()*m_Uil_;
      m_Gammaklold_ = m_Gammakl_;
      m_Gammakl_ = m_Ykl_.array()/(m_Tik_.transpose()*v_Mui_*v_nul_.transpose()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
        break;
      }
    }
  }
  catch (char const *e) {
#ifdef RPACKAGE
    R_errormsg = e;
#endif
#ifdef COVERBOSE
    std::cerr<<e<<"\n";
#endif
    return false;
  }

  catch (...) {
#ifdef RPACKAGE
    R_errormsg = "Unknown error occurred..Terminating Algorithm";
#endif

#ifdef COVERBOSE
    std::cerr<<"Unknown error occurred..Terminating Algorithm"<<"\n";
#endif
    return false;
  }
  return true;
}


bool ContingencyLBModel_mu_i_nu_j::CEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  //Temporary variables
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);

  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      PoissonLogsumRows(m_sumik);
      Zsumk.setZero();
      for (int i = 0; i < nbSample_; ++i) {
        v_sumikmax(i) = m_sumik.row(i).maxCoeff(&maxIndex);
        Zsumk(maxIndex)+=1;
      }

      if((Zsumk.array()<0.00001).any())
      {
        empty_cluster = true;
        throw "Row clustering failed while running model.";
      }
      else
      {
        empty_cluster = false;
      }

      m_prodik=(m_sumik-v_sumikmax*Onesk.transpose()).array().exp();
      v_sumi = m_prodik.rowwise().sum();
      m_Tik_ = m_prodik.array()/(v_sumi*Onesk.transpose()).array();
      v_Tk_ = m_Tik_.colwise().sum();

      //M-step
      //M-step
      if(!Mparam_.fixedproportions_) {
        v_logPiek_=(v_Tk_.array()/nbSample_).log();
      }
      m_Ykl_ = m_Tik_.transpose()*m_Uil_;
      m_Gammaklold_ = m_Gammakl_;
      m_Gammakl_ = m_Ykl_.array()/(m_Tik_.transpose()*v_Mui_*v_nul_.transpose()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
        break;
      }
    }
  }
  catch (char const *e) {
#ifdef RPACKAGE
    R_errormsg = e;
#endif
#ifdef COVERBOSE
    std::cerr<<e<<"\n";
#endif
    return false;
  }

  catch (...) {
#ifdef RPACKAGE
    R_errormsg = "Unknown error occurred..Terminating Algorithm";
#endif

#ifdef COVERBOSE
    std::cerr<<"Unknown error occurred..Terminating Algorithm"<<"\n";
#endif
    return false;
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::EMCols()
{
  //Initialization
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;

  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_),Wsuml(Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);

  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      PoissonLogsumCols(m_sumjl);
      Wsuml.setZero();
      for (int j = 0; j < nbVar_; ++j) {
        v_sumjlmax(j) = m_sumjl.row(j).maxCoeff(&maxIndex);
        Wsuml(maxIndex)+=1;
      }

      //Check for empty cluster
      if((Wsuml.array()<.00001).any()){
        empty_cluster = true;
        throw "Column clustering failed while running model.";
      }
      else
      {empty_cluster = false;}


      m_prodjl=(m_sumjl-v_sumjlmax*Onesl.transpose()).array().exp();
      v_sumj = m_prodjl.rowwise().sum();
      m_Rjl_ = m_prodjl.array()/(v_sumj*Onesl.transpose()).array();
      v_Rl_ = m_Rjl_.colwise().sum();

      //M-step
      if(!Mparam_.fixedproportions_) {
        v_logRhol_=(v_Rl_.array()/nbVar_).log();
      }

      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Gammaklold_ = m_Gammakl_;
      m_Gammakl_ = m_Ykl_.array()/(v_muk_*v_Nuj_.transpose()*m_Rjl_).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
        break;
      }
    }
  }
  catch (char const *e) {
#ifdef RPACKAGE
    R_errormsg = e;
#endif
#ifdef COVERBOSE
    std::cerr<<e<<"\n";
#endif
    return false;
  }

  catch (...) {
#ifdef RPACKAGE
    R_errormsg = "Unknown error occurred..Terminating Algorithm";
#endif

#ifdef COVERBOSE
    std::cerr<<"Unknown error occurred..Terminating Algorithm"<<"\n";
#endif
    return false;
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::CEMCols()
{
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  //CEMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      PoissonLogsumCols(m_sumjl);
      m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
      for (int j = 0; j < nbVar_; ++j) {
        m_sumjl.row(j).maxCoeff(&maxIndex);
        m_Rjl_(j,maxIndex)=1;
      }
      v_Rl_ = m_Rjl_.colwise().sum();
      if((v_Rl_.array()<.00001).any()){
        empty_cluster = true;
        throw "Column clustering failed while running model.";
      }
      else
      {empty_cluster = false;}

      //M-step
      if(!Mparam_.fixedproportions_) {
        v_logRhol_=(v_Rl_.array()/nbVar_).log();
      }
      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Gammaklold_ = m_Gammakl_;
      m_Gammakl_ = m_Ykl_.array()/(v_muk_*v_Nuj_.transpose()*m_Rjl_).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array()).abs()/m_Gammakl_.array()).sum())<Mparam_.epsilon_int_) {
        break;
      }
    }
  }
  catch (char const *e) {
#ifdef RPACKAGE
    R_errormsg = e;
#endif
#ifdef COVERBOSE
    std::cerr<<e<<"\n";
#endif
    return false;
  }

  catch (...) {
#ifdef RPACKAGE
    R_errormsg = "Unknown error occurred..Terminating Algorithm";
#endif

#ifdef COVERBOSE
    std::cerr<<"Unknown error occurred..Terminating Algorithm"<<"\n";
#endif
    return false;
  }
  return true;
}

void ContingencyLBModel_mu_i_nu_j::PoissonLogsumRows(MatrixReal & m_ik)
{
  m_ik = VectorReal::Ones(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.array().log()).matrix().transpose())
      -v_Mui_*(m_Gammakl_*v_nul_).transpose();
}

void ContingencyLBModel_mu_i_nu_j::PoissonLogsumCols(MatrixReal & m_jl)
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

bool ContingencyLBModel_mu_i_nu_j::FuzzyCEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Fuzzy CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Fuzzy CEM initialization is not valid for this model.";
#endif
  return false;
}
bool ContingencyLBModel_mu_i_nu_j::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "CEM initialization is not valid for this model.";
#endif
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Running Initialization.."<<"\n";
#endif
  //Initialize random row and column partition
  try {
    VectorReal probarows = (1.0/Mparam_.nbrowclust_)*VectorReal::Ones(Mparam_.nbrowclust_);
    VectorReal probacols = (1.0/Mparam_.nbcolclust_)*VectorReal::Ones(Mparam_.nbcolclust_);
    v_Zi_ = PartRnd(nbSample_,probarows);
    v_Wj_ = PartRnd(nbVar_,probacols);
    m_Tik_.setZero(nbSample_,Mparam_.nbrowclust_);
    m_Rjl_.setZero(nbVar_,Mparam_.nbcolclust_);

    for (int i = 0; i < nbSample_; ++i) {
      m_Tik_(i,v_Zi_(i)-1) = 1.0;
    }

    for (int j = 0; j < nbVar_; ++j) {
      m_Rjl_(j,v_Wj_(j)-1) = 1.0;
    }

    m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_).array()/(m_Tik_.transpose()*v_Mui_*v_Nuj_.transpose()*m_Rjl_).array();

    //Initializing model parameters
    v_Rl_ = m_Rjl_.colwise().sum().transpose();
    m_Gammakl1_ = m_Gammakl_;
    m_Gammakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m_Gammaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m_Uil_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
    m_Vjk_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
    v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
    v_Zi_ = MatrixInteger::Zero(nbSample_,1);
    v_Wj_ = MatrixInteger::Zero(nbVar_,1);
    v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(VectorReal::Ones(Mparam_.nbrowclust_));
    v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(VectorReal::Ones(Mparam_.nbcolclust_));
  }
  catch (char const *e) {
#ifdef RPACKAGE
    R_errormsg = e;
#endif
#ifdef COVERBOSE
    std::cerr<<e<<"\n";
#endif
    return false;
  }

  catch (...) {
#ifdef RPACKAGE
    R_errormsg = "Unknown error occurred during initialization.";
#endif

#ifdef COVERBOSE
    std::cerr<<"Unknown error occurred during initialization."<<"\n";
#endif
    return false;
  }
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

  for (int ind = 0; ind < clusters; ++ind) {
    v_Z(v_Randperm(ind)) = ind+1;
  }
  remainingclusters = (clusters+1)*MatrixInteger::Ones(n-clusters,1) - ((MatrixReal::Ones(clusters,1)*(Unifrnd(0,1,1,n-clusters))).array() < (Cumsum(proba)*MatrixReal::Ones(1,n-clusters)).array()).matrix().cast<int>().colwise().sum().transpose();

  for (int ind = clusters; ind < n; ++ind) {
    v_Z(v_Randperm(ind)) = remainingclusters(ind-clusters)==(clusters+1)?clusters:remainingclusters(ind-clusters);
  }
  return v_Z;

}

VectorReal ContingencyLBModel_mu_i_nu_j::Cumsum(VectorReal proba)
{
  int size = proba.rows();
  VectorReal v_temp = VectorReal::Zero(size);
  v_temp(0) = proba(0);
  for (int itr = 1; itr < size; ++itr) {
    v_temp(itr) = v_temp(itr-1) + proba(itr);
  }
  return v_temp;
}

MatrixReal ContingencyLBModel_mu_i_nu_j::Unifrnd(float a,float b, int row, int col)
{
  MatrixReal m_temp(row,col);
  for (int r = 0; r < row; ++r) {
    for (int c = 0; c < col; ++c) {
      m_temp(r,c) = (b-a)*(std::rand()/float(RAND_MAX)) + a;
    }
  }
  return m_temp;
}

void ContingencyLBModel_mu_i_nu_j::FinalizeOutput()
{
  // Calculate row and column proportions
  if(!Mparam_.fixedproportions_){
    v_Piek_ = v_logPiek_.array().exp();
    v_Rhol_ = v_logRhol_.array().exp();
  }else
  {
    v_Piek_ = (1.0/Mparam_.nbrowclust_)*MatrixReal::Ones(Mparam_.nbrowclust_,1);
    v_Rhol_ = (1.0/Mparam_.nbcolclust_)*MatrixReal::Ones(Mparam_.nbcolclust_,1);;
  }

}

void ContingencyLBModel_mu_i_nu_j::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

MatrixInteger ContingencyLBModel_mu_i_nu_j::GetArrangedDataClusters()
{
  Eigen::ArrayXi v_Zi = (GetRowClassificationVector()).array();
  Eigen::ArrayXi v_Wj = (GetColumnClassificationVector()).array();
  MatrixInteger m_clusterij = MatrixInteger::Zero(nbSample_,nbVar_);

  //Rearrange data into clusters

  VectorInteger rowincrement = MatrixInteger::Zero(Mparam_.nbrowclust_,1);

  VectorInteger nbindrows = MatrixInteger::Zero(Mparam_.nbrowclust_,1);
  for (int k = 1; k < Mparam_.nbrowclust_; ++k) {
    nbindrows(k) = (v_Zi==(k-1)).count()+nbindrows(k-1);
  }

  VectorInteger colincrement = MatrixInteger::Zero(Mparam_.nbcolclust_,1);

  VectorInteger nbindcols = MatrixInteger::Zero(Mparam_.nbcolclust_,1);
  for (int l = 1; l < Mparam_.nbcolclust_; ++l) {
    nbindcols(l)=(v_Wj==(l-1)).count()+nbindcols(l-1);
  }

  for (int j = 0; j < nbVar_; ++j) {
    m_clusterij.col(colincrement(v_Wj(j)) + nbindcols(v_Wj(j))) = m_Dataij_.col(j);
    colincrement(v_Wj(j))+=1;
  }
  MatrixInteger temp = m_clusterij;

  for (int i = 0; i < nbSample_; ++i) {
    m_clusterij.row( rowincrement(v_Zi(i)) + nbindrows(v_Zi(i))) = temp.row(i);
    rowincrement(v_Zi(i))+=1;
  }

  return m_clusterij;

}
