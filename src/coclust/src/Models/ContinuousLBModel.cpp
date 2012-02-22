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


/** @file ContinuousLBModel.cpp
 *  @brief Implements concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#include "ContinuousLBModel.h"
#include <exception>
ContinuousLBModel::ContinuousLBModel( MatrixReal const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  dimprod_ = nbSample_*nbVar_;
};

ContinuousLBModel::ContinuousLBModel(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,
                                     VectorInteger const & collabels,ModelParameters const& Mparam)
               : ICoClustModel(Mparam,rowlabels,collabels)
               , m_Dataij_(m_Dataij)
{
  dimprod_ = nbSample_*nbVar_;
};

void ContinuousLBModel::LogSumRows(MatrixReal& m_sumik)
{
  VectorReal tmp1(Mparam_.nbrowclust_);
  tmp1 = 0.5*(((m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix())*v_Rl_);
  tmp1 -= v_logPiek_;
  m_sumik = -VectorReal::Ones(nbSample_)*tmp1.transpose()
            -0.5*(m_Uil2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix().transpose()))
            +(m_Uil1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix().transpose()));

}

void ContinuousLBModel::LogSumCols(MatrixReal & m_sumjl)
{
  VectorReal tmp1(Mparam_.nbcolclust_);
  tmp1 = 0.5*((v_Tk_.transpose())*(m_Sigma2kl_.array().log()+m_Mukl2_.array()/m_Sigma2kl_.array()).matrix());
  tmp1 -= v_logRhol_;
  m_sumjl = -VectorReal::Ones(nbVar_)*tmp1.transpose()
          -0.5*(m_Vjk2_*((Array2DReal::Ones(Mparam_.nbrowclust_,Mparam_.nbcolclust_)/m_Sigma2kl_.array()).matrix()))
          +m_Vjk1_*((m_Mukl_.array()/m_Sigma2kl_.array()).matrix());

}

bool ContinuousLBModel::EMRows()
{
  //Initialization
  m_Mukl2_ = m_Mukl_.array().square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ERows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    MStepRows();
    //Termination check
    if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }

  }
  return true;
}

bool ContinuousLBModel::SEMRows()
{
  //Initialization
  m_Mukl2_ = m_Mukl_.array().square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  if(!SERows()) return false;
  //M-step
  MStepRows();
  return true;
}

bool ContinuousLBModel::EMCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ECols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    MStepCols();
    //Termination check
    if((((m_Mukl_-m_Muklold2_).array()/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }

  }
  return true;
}

bool ContinuousLBModel::SEMCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  if(!SECols()) return false;
  //M-step
  MStepCols();
  return true;
}


bool ContinuousLBModel::CEMRows()
{
  //Initialization
  m_Mukl2_ = m_Mukl_.array().square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!CERows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    MStepRows();
    //Termination check
    if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }

  }
  return true;
}


bool ContinuousLBModel::CEMCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!CECols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    MStepCols();
    //Termination check
    if((((m_Mukl_.array()-m_Muklold2_.array())/m_Mukl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }

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
  if (InitCEMCols()) {
    v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(VectorReal::Ones(Mparam_.nbrowclust_));
    v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(VectorReal::Ones(Mparam_.nbcolclust_));
    m_Dataij2_ = m_Dataij_.array().square();
    m_Uil1_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
    m_Uil2_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
    m_Vjk1_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
    m_Vjk2_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
    v_Tk_ = VectorReal::Zero(Mparam_.nbrowclust_);
    m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
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

// Private initialization  functions

bool ContinuousLBModel::InitCEMCols()
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
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for ( Label_pair : knownLabelsCols_) {
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
#else
  for ( int j=0;j<knownLabelsCols_.size();j++) {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
#endif
  m_Sigma2kl_ = MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  //Initializations
  SelectRandomRowsfromdata(m_Vkj);
  m_Vjk = m_Vkj.transpose();
  VectorReal v_Vj2 = (m_Vjk.array().square()).rowwise().sum();
  GenerateRandomMean(m_Vjk,m_Mulk);

  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
    m_Mul2 = ((m_Mulk.array().square()).rowwise().sum()).matrix().transpose();
    m_Djl = v_Vj2*MatrixReal::Ones(1,Mparam_.nbcolclust_) + MatrixReal::Ones(nbVar_,1)*m_Mul2 - 2*(m_Vjk*m_Mulk.transpose());
#ifdef RANGEBASEDFORLOOP
    for ( int j : UnknownLabelsCols_) {
      m_Djl.row(j).minCoeff(&minIndex);
      m_Rjl_.row(j).setZero();
      m_Rjl_(j,minIndex)=1;
    }
#else
    for ( int j =0;j< UnknownLabelsCols_.size();j++) {
      m_Djl.row(UnknownLabelsCols_[j]).minCoeff(&minIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZero();
      m_Rjl_(UnknownLabelsCols_[j],minIndex)=1;
    }
#endif

    v_Rl_ = m_Rjl_.colwise().sum();
    //Check for empty cluster
    if((v_Rl_.array()<.00001).any()){
      Error_msg_  = "Column clustering failed while running initialization.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  empty_cluster = true;
  return false;
    }else{empty_cluster = false;}

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

  return true;
}

void ContinuousLBModel::FinalizeOutput()
{
  CommonFinalizeOutput();
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
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
    //_m_kj.row(k) = m_Dataij_.row(v_temp(k));
    _m_kj.row(k) = m_Dataij_.row(k);
  }
}

void ContinuousLBModel::GenerateRandomMean(const MatrixReal & m_jk, MatrixReal & mean_lk)
{
  VectorInteger v_temp= RandSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l) {
    //mean_lk.row(l) = m_jk.row(v_temp(l));
    mean_lk.row(l) = m_jk.row(l);
  }
}

const MatrixReal& ContinuousLBModel::GetArrangedDataClusters()
{
  ArrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContinuousLBModel::Modify_theta_start()
{
  m_Muklstart_ = m_Mukl_;
  m_Sigma2klstart_ = m_Sigma2kl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

void ContinuousLBModel::Copy_theta_start()
{
  m_Mukl_ = m_Muklstart_;
  m_Sigma2kl_ = m_Sigma2klstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  // initializations
  m_Mukl2_ = m_Mukl_.array().square();
  m_Muklold1_ = m_Mukl_;
  v_Rl_ = m_Rjl_.colwise().sum();
}

void ContinuousLBModel::Copy_theta_max()
{
  m_Mukl_ = m_Muklmax_;
  m_Sigma2kl_ = m_Sigma2klmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

void ContinuousLBModel::Modify_theta_max()
{
  m_Muklmax_ = m_Mukl_;
  m_Sigma2klmax_ = m_Sigma2kl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}

void ContinuousLBModel::MStepFull()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  MatrixReal m_trkl = (v_Tk_*v_Rl_.transpose());
  m_Mukl_ = (m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_).array()/m_trkl.array();
  m_Mukl2_ = m_Mukl_.array().square();
  m_Sigma2kl_ = (m_Tik_.transpose()*m_Dataij2_.cast<float>()*m_Rjl_).array()/m_trkl.array() - m_Mukl2_.array();
}

