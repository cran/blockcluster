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
 * created on: Feb 3, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ContinuousLBModel.cpp
 *  @brief Implements concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#include "ContinuousLBModel.h"
#include <exception>
ContinuousLBModel::ContinuousLBModel( MatrixReal const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  nbSample_ = m_Dataij.rows();
  nbVar_ = m_Dataij.cols();
  Mparam_.nbrowdata_ = nbSample_;
  Mparam_.nbcoldata_ = nbVar_;
  dimprod_ = nbSample_*nbVar_;
};

bool ContinuousLBModel::Estep()
{
  if (EMRows()) {
    if (EMCols()) {
      return true;
    }
  }
  return false;
}

bool ContinuousLBModel::CEstep()
{
  if (CEMRows()) {
    if (CEMCols()) {
      return true;
    }
  }
  return false;
}


void ContinuousLBModel::Mstep()
{
  // Doing nothing here
}

bool ContinuousLBModel::EMRows()
{
  //Initialization
  m_Mukl2_ = m_Mukl_.array().square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;
  //Temporary variables
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_),m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);
  VectorReal tmp1(Mparam_.nbrowclust_);
  //EM algorithm
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      tmp1 = 0.5*(((m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix())*v_Rl_);
      tmp1 -= v_logPiek_;
      m_sumik = -VectorReal::Ones(nbSample_)*tmp1.transpose()
                -0.5*(m_Uil2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix().transpose()))
                +(m_Uil1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix().transpose()));
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


      m_Muklold2_ = m_Mukl_;
      m_trkl = (v_Tk_*v_Rl_.transpose());
      m_Mukl_ = (m_Tik_.transpose()*m_Uil1_).array()/m_trkl.array();
      m_Mukl2_ = m_Mukl_.array().square();
      m_Sigma2kl_ = (m_Tik_.transpose()*m_Uil2_).array()/m_trkl.array() - m_Mukl2_.array();

      //Termination check
      if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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

bool ContinuousLBModel::EMCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_),m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_),Wsuml(Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);
  VectorReal tmp1(Mparam_.nbcolclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //GaussianLogsumCols(m_sumjl);
      tmp1 = 0.5*((v_Tk_.transpose())*(m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix());
      tmp1 -= v_logRhol_;
      m_sumjl = -VectorReal::Ones(nbVar_)*tmp1.transpose()
              -0.5*(m_Vjk2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix()))
              +m_Vjk1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix());
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

      m_Muklold2_ = m_Mukl_;
      m_trkl = v_Tk_*v_Rl_.transpose();
      m_Mukl_ = (m_Vjk1_.transpose()*m_Rjl_).array()/m_trkl.array();
      m_Mukl2_ = m_Mukl_.array().square();
      m_Sigma2kl_ = (m_Vjk2_.transpose()*m_Rjl_).array()/m_trkl.array() - m_Mukl2_.array();

      //Termination check
      if((((m_Mukl_-m_Muklold2_).array()/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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


bool ContinuousLBModel::CEMRows()
{
  //Initialization
  m_Mukl2_ = m_Mukl_.array().square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;
  //Temporary variables
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_),m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);
  VectorReal tmp1(Mparam_.nbrowclust_);
  //EM algorithm
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      tmp1 = 0.5*(((m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix())*v_Rl_);
      tmp1 -= v_logPiek_;
      m_sumik = -VectorReal::Ones(nbSample_)*tmp1.transpose()
                -0.5*(m_Uil2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix().transpose()))
                +(m_Uil1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix().transpose()));
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


      m_Muklold2_ = m_Mukl_;
      m_trkl = (v_Tk_*v_Rl_.transpose());
      m_Mukl_ = (m_Tik_.transpose()*m_Uil1_).array()/m_trkl.array();
      m_Mukl2_ = m_Mukl_.array().square();
      m_Sigma2kl_ = (m_Tik_.transpose()*m_Uil2_).array()/m_trkl.array() - m_Mukl2_.array();

      //Termination check
      if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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


bool ContinuousLBModel::CEMCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

  //Temporary variables
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_),m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  VectorReal tmp1(Mparam_.nbcolclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //GaussianLogsumCols(m_sumjl);
      tmp1 = 0.5*((v_Tk_.transpose())*(m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix());
      tmp1 -= v_logRhol_;
      m_sumjl = -VectorReal::Ones(nbVar_)*tmp1.transpose()
              -0.5*(m_Vjk2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix()))
              +m_Vjk1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix());
      m_Rjl_.setZero();
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

      m_Muklold2_ = m_Mukl_;
      m_trkl = (v_Tk_*v_Rl_.transpose());
      m_Mukl_ = (m_Vjk1_.transpose()*m_Rjl_).array()/m_trkl.array();
      m_Mukl2_ = m_Mukl_.array().square();
      m_Sigma2kl_ = (m_Vjk2_.transpose()*m_Rjl_).array()/m_trkl.array(); - m_Mukl2_.array();

      //Termination check
      if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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

void ContinuousLBModel::ParameterStopCriteria()
{
    if((((m_Mukl_.array()-m_Muklold1_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_) {
    StopAlgo = true;
  }
  else
  {
    StopAlgo = false;
    m_Muklold1_ =m_Mukl_;
  }
}

float ContinuousLBModel::EstimateLikelihood()
{
  Likelihood_ = (-0.5*(dimprod_+v_Tk_.transpose()*(m_Sigma2kl_.array().log()).matrix()*v_Rl_) + v_Tk_.transpose()*v_logPiek_
          + v_Rl_.transpose()*v_logRhol_ -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
          -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum())/dimprod_;

  return Likelihood_;
}

void ContinuousLBModel::likelihoodStopCriteria()
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
// Initialization using CEM algo

bool ContinuousLBModel::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Running Initialization.."<<"\n";
#endif
  if (InitCEMRows()) {
    v_Rl_ = m_Rjl_.colwise().sum();
    v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(VectorReal::Ones(Mparam_.nbrowclust_));
    v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(VectorReal::Ones(Mparam_.nbcolclust_));
    m_Dataij2_ = m_Dataij_.array().square();
    m_Uil1_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
    m_Uil2_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
    m_Vjk1_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
    m_Vjk2_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
    v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
    m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
    v_Zi_ = MatrixInteger::Zero(nbSample_,1);
    m_Zik_ = MatrixInteger::Zero(nbSample_,Mparam_.nbrowclust_);
    m_Wjl_ = MatrixInteger::Zero(nbVar_,Mparam_.nbcolclust_);
    v_Wj_ = MatrixInteger::Zero(nbVar_,1);
    m_Muklold1_ = m_Mukl_;
    m_Mukl2_ = m_Mukl_.array().square();
    m_Muklold2_ = MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
#ifdef COVERBOSE
    std::cout<<"Initialization over."<<"\n";
#endif
    return true;
  }
    return false;
}

bool ContinuousLBModel::FuzzyCEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Fuzzy CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Fuzzy CEM initialization is not valid for this model.";
#endif
  return false;
}
bool ContinuousLBModel::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Random initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Random initialization is not valid for this model.";
#endif
  return false;
}


// Private initialization  functions

bool ContinuousLBModel::InitCEMRows()
{

  //Temporary variables
  MatrixReal m_Vjk(nbVar_,Mparam_.nbrowclust_);
  MatrixReal m_Vkj(Mparam_.nbrowclust_,nbVar_);
  MatrixReal m_Djl(nbVar_,Mparam_.nbcolclust_);
  MatrixReal m_Mulk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  MatrixReal m_Mul2(1,Mparam_.nbcolclust_);
  MatrixReal m_Rlk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  VectorReal::Index minIndex;
  float W1 = RealMax , W1_old;

  //Private members initialization
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
  m_Sigma2kl_ = MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  //Initializations
  SelectRandomRowsfromdata(m_Vkj);
  m_Vjk = m_Vkj.transpose();
  VectorReal v_Vj2 = (m_Vjk.array().square()).rowwise().sum();
  GenerateRandomMean(m_Vjk,m_Mulk);
  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      m_Mul2 = ((m_Mulk.array().square()).rowwise().sum()).matrix().transpose();
      m_Djl = v_Vj2*MatrixReal::Ones(1,Mparam_.nbcolclust_) + MatrixReal::Ones(nbVar_,1)*m_Mul2 - 2*(m_Vjk*m_Mulk.transpose());
      m_Rjl_.setZero(nbVar_,Mparam_.nbcolclust_);
      for (int j = 0; j < nbVar_; ++j) {
        m_Djl.row(j).minCoeff(&minIndex);
        m_Rjl_(j,minIndex)=1;
      }

      v_Rl_ = m_Rjl_.colwise().sum();
      //Check for empty cluster
      if((v_Rl_.array()<.00001).any()){
        empty_cluster = true;
        throw "Column clustering failed inside Initialization";
      }
      else
      {empty_cluster = false;}

      m_Mulk = (m_Rjl_.transpose()*m_Vjk).array()/(v_Rl_*MatrixReal::Ones(1,Mparam_.nbrowclust_)).array();
      //Check for termination
      W1_old = W1;
      W1 = (m_Rjl_.array()*m_Djl.array()).sum();
      //Initialization of model parameters
      m_Mukl_ = m_Mulk.transpose();
      m_Sigma2kl_ = (W1/(Mparam_.nbrowclust_*Mparam_.nbcolclust_))*MatrixReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      if (std::abs((W1-W1_old)/W1)<Mparam_.initepsilon_) {
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

void ContinuousLBModel::FinalizeOutput()
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

void ContinuousLBModel::ConsoleOut()
{
#ifndef RPACKAGE
    std::cout<<"Output Model parameter:"<<"\nBlock Mean:\n"<<m_Mukl_<<"\nBlock Sigma:\n"<<m_Sigma2kl_<<"\npiek: "<<
        v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

void ContinuousLBModel::SelectRandomRowsfromdata(MatrixReal& _m_kj)
{
  VectorInteger v_temp= RandSample(nbSample_,Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    _m_kj.row(k) = m_Dataij_.row(v_temp(k));
    //_m_kj.row(k) = m_Dataij_.row(k);
  }
}

void ContinuousLBModel::GenerateRandomMean(const MatrixReal & m_jk, MatrixReal & mean_lk)
{
  VectorInteger v_temp= RandSample(nbVar_,Mparam_.nbcolclust_);
  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    mean_lk.row(l) = m_jk.row(v_temp(l));
    //mean_lk.row(l) = m_jk.row(l);
  }
}

MatrixReal ContinuousLBModel::GetArrangedDataClusters()
{
  Eigen::ArrayXi v_Zi = (GetRowClassificationVector()).array();
  Eigen::ArrayXi v_Wj = (GetColumnClassificationVector()).array();
  MatrixReal m_clusterij = MatrixReal::Zero(nbSample_,nbVar_);

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
  MatrixReal temp = m_clusterij;

  for (int i = 0; i < nbSample_; ++i) {
    m_clusterij.row( rowincrement(v_Zi(i)) + nbindrows(v_Zi(i))) = temp.row(i);
    rowincrement(v_Zi(i))+=1;
  }

  return m_clusterij;
}
