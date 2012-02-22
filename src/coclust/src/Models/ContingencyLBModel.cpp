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

/** @file ContingencyLBModel.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#include "ContingencyLBModel.h"

ContingencyLBModel::ContingencyLBModel(MatrixInteger const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  DataSum_ = m_Dataij.array().sum();
};

ContingencyLBModel::ContingencyLBModel(MatrixInteger const& m_Dataij,VectorInteger const & rowlabels,
                                       VectorInteger const & collabels,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam,rowlabels,collabels)
                           , m_Dataij_(m_Dataij)
{
  DataSum_ = m_Dataij.array().sum();
};


bool ContingencyLBModel::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (InitCEMRows()) {
    if (InitCEMCols()) {
      m_Gammakl1_ = m_Gammakl_;
      m_Gammakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Gammaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Uil_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
      m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(VectorReal::Ones(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(VectorReal::Ones(Mparam_.nbcolclust_));
#ifdef COVERBOSE
  std::cout<<"Initialization over."<<std::endl;
#endif
      return true;
    }
  }
  return false;
}

bool ContingencyLBModel::EMRows()
{
  //Initialization
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_Yl_ = m_Uil_.colwise().sum().transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ERows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepRows();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel::SEMRows()
{
  //Initialization
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_Yl_ = m_Uil_.colwise().sum().transpose();

  if(!SERows()) return false;
  //M-step
  MStepRows();
  return true;
}

bool ContingencyLBModel::CEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_Yl_ = m_Uil_.colwise().sum().transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    if(!CERows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepRows();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel::EMCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_Yk_ = m_Vjk_.colwise().sum().transpose();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ECols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepCols();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel::SEMCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_Yk_ = m_Vjk_.colwise().sum().transpose();

  if(!SECols()) return false;
  //M-step
  MStepCols();
  return true;
}

bool ContingencyLBModel::CEMCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_Yk_ = m_Vjk_.colwise().sum().transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!CECols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    MStepCols();

    //Termination check
    if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

float ContingencyLBModel::EstimateLikelihood()
{
  Likelihood_ = (m_Ykl_.array()*(m_Gammakl_.array()).log()).sum() - DataSum_ + v_Tk_.transpose()*v_logPiek_ + v_Rl_.transpose()*v_logRhol_
            -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
            -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum();
  return Likelihood_;
}

void ContingencyLBModel::likelihoodStopCriteria()
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

void ContingencyLBModel::ParameterStopCriteria()
{
  float relativechange = (((m_Gammakl1_.array()-m_Gammakl1old_.array())/m_Gammakl1_.array()).abs().sum());
  if(relativechange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

void ContingencyLBModel::SelectRandomColsfromdata(MatrixReal& _m_il,int cols)
{
  if(cols==nbVar_)
    _m_il = m_Dataij_.cast<float>();
  else{
    //random shuffle Algorithm
    int random,index,temp;
    VectorInteger _v_temp(Mparam_.nbcoldata_);
    for ( int j = 0; j < Mparam_.nbcoldata_; ++j) {
      _v_temp(j)=j;
    }
    for ( int l = 0; l < cols; ++l){
      random=std::rand()%(Mparam_.nbcoldata_-l);
      index=_v_temp(random);
      _m_il.col(l)=m_Dataij_.cast<float>().col(index);
      //swap elements
      temp=_v_temp(Mparam_.nbcoldata_-l-1);
      _v_temp(Mparam_.nbcoldata_-l-1)=_v_temp(random);
      _v_temp(random)=temp;
    }
  }
}

bool ContingencyLBModel::InitCEMRows()
{
  //Temporary variables
  int cols=std::min(100,int(nbVar_));
  MatrixReal m_Akl(Mparam_.nbrowclust_,cols),m_Aklold(Mparam_.nbrowclust_,cols);
  MatrixReal m_sumik(Mparam_.nbrowdata_,Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;

  //Model parameters
  m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Rl_ = VectorReal::Ones(cols);
  v_Tk_ = VectorReal::Zero(Mparam_.nbrowclust_);
  m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);

  // Initializations

  m_Uil_ = MatrixReal::Zero(nbSample_,cols);
  SelectRandomColsfromdata(m_Uil_,cols);
  v_Ui_ = m_Uil_.rowwise().sum();
  GenerateRandomPoissonParameterRows(m_Akl,cols);
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for(Label_pair : knownLabelsRows_){
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }
#else
  for(int i=0;i<knownLabelsRows_.size();i++){
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }
#endif
  //Determine row partition using CEM algorithm with equal proportions
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
    m_sumik = m_Uil_*(m_Akl.transpose());
#ifdef RANGEBASEDFORLOOP
    for ( int i : UnknownLabelsRows_) {
      m_sumik.row(i).maxCoeff(&maxIndex);
      m_Tik_.row(i).setZero();
      m_Tik_(i,maxIndex)=1;
    }
#else
    for ( int i =0;i< UnknownLabelsRows_.size();i++) {
      m_sumik.row(UnknownLabelsRows_[i]).maxCoeff(&maxIndex);
      m_Tik_.row(UnknownLabelsRows_[i]).setZero();
      m_Tik_(UnknownLabelsRows_[i],maxIndex)=1;
    }
#endif
      v_Tk_ = m_Tik_.colwise().sum();
      if((v_Tk_.array()<.00001).any()){
        Error_msg_  = "Row clustering failed while running initialization.";
#ifdef COVERBOSE
    std::cout<<Error_msg_<<"\n";
#endif
    empty_cluster = true;
    return false;
      }else{empty_cluster = false;}

    // M-step
    m_Aklold = m_Akl;
    m_Akl = ((m_Tik_.transpose()*m_Uil_).array()/((m_Tik_.transpose()*v_Ui_)*MatrixReal::Ones(1,cols)).array()+RealMin).log();
    if((((m_Akl.array()-m_Aklold.array()).abs()/m_Akl.array()).sum())<Mparam_.initepsilon_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel::InitCEMCols()
{
  //Temporary variables
  MatrixReal m_Alk(Mparam_.nbcolclust_,Mparam_.nbrowclust_) , m_Alkold(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;

  //Initializations
  m_Vjk_ = (m_Dataij_.cast<float>().transpose())*m_Tik_;
  v_Vj_ = m_Vjk_.rowwise().sum();
  GenerateRandomPoissonParameterCols(m_Alk);
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for (Label_pair : knownLabelsCols_) {
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
#else
  for ( int j=0;j<knownLabelsCols_.size();j++) {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
#endif

  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
    // CE-step
    m_sumjl = m_Vjk_*(m_Alk.transpose());
#ifdef RANGEBASEDFORLOOP
    for ( int j : UnknownLabelsCols_) {
      m_sumjl.row(j).maxCoeff(&maxIndex);
      m_Rjl_.row(j).setZero();
      m_Rjl_(j,maxIndex)=1;
    }
#else
    for ( int j =0;j< UnknownLabelsCols_.size();j++) {
      m_sumjl.row(UnknownLabelsCols_[j]).maxCoeff(&maxIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZero();
      m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
    }
#endif

    v_Rl_ = m_Rjl_.colwise().sum();
    if((v_Rl_.array()<.00001).any()){
      Error_msg_  = "Column clustering failed while running initialization.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  empty_cluster = true;
  return false;
    }else{empty_cluster = false;}
    // M-step
    m_Alkold = m_Alk;
    m_Alk = ((m_Rjl_.transpose()*m_Vjk_).array()/((m_Rjl_.transpose()*v_Vj_)*MatrixReal::Ones(1,Mparam_.nbrowclust_)).array()+RealMin).log();
    if((((m_Alk.array()-m_Alkold.array()).abs()/m_Alk.array()).sum())<Mparam_.initepsilon_) {
      break;
    }
  }
  MatrixReal m_Ykl = m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_;
  m_Gammakl_ = MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  m_Gammakl_ = m_Ykl.array()/((m_Ykl.rowwise().sum())*(m_Ykl.colwise().sum())).array();
  return true;
}

void ContingencyLBModel::LogSumRows(MatrixReal & m_ik)
{
  m_ik = VectorReal::Ones(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.array().log()).matrix().transpose());
}

void ContingencyLBModel::LogSumCols(MatrixReal & m_jl)
{
  m_jl = VectorReal::Ones(nbVar_)*v_logRhol_.transpose() +
      m_Vjk_*((m_Gammakl_.array().log()).matrix());
}

void ContingencyLBModel::GenerateRandomPoissonParameterRows(MatrixReal& m_kl,int cols)
{
  int index;
  VectorInteger _v_temp = RandSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k){
    index=_v_temp(k);
    //index=k;
    for ( int l = 0; l < cols; ++l) {
      m_kl(k,l)=std::log(m_Uil_(index,l)/v_Ui_(index)+RealMin);
    }
  }
}

void ContingencyLBModel::GenerateRandomPoissonParameterCols(MatrixReal& m_lk)
{
  int index;
  VectorInteger _v_temp = RandSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l){
    index=_v_temp(l);
    //index=l;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
      m_lk(l,k) = std::log(m_Vjk_(index,k)/v_Vj_(index)+RealMin);
    }
  }
}

void ContingencyLBModel::FinalizeOutput()
{
  CommonFinalizeOutput();
}

void ContingencyLBModel::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}
const MatrixInteger& ContingencyLBModel::GetArrangedDataClusters()
{
  ArrangedDataCluster<MatrixInteger>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContingencyLBModel::Modify_theta_start()
{
  m_Gammaklstart_ = m_Gammakl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
}

void ContingencyLBModel::Copy_theta_start()
{
  m_Gammakl_ = m_Gammaklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
  m_Gammakl1_ = m_Gammakl_;
}

void ContingencyLBModel::Copy_theta_max()
{
  m_Gammakl_ = m_Gammaklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

void ContingencyLBModel::Modify_theta_max()
{
  m_Gammaklmax_ = m_Gammakl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}

void ContingencyLBModel::MStepFull()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_;
  m_Gammakl_ = m_Ykl_.array()/((m_Ykl_).rowwise().sum()*(m_Dataij_.cast<float>()*m_Rjl_).colwise().sum()).array();
}
