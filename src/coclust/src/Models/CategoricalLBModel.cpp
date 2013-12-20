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

#include <math.h>
#include "CategoricalLBModel.h"

CategoricalLBModel::CategoricalLBModel(MatrixInteger const& m_Dataij,ModelParameters const& Mparam
                                       ,int a,int b)
:ICoClustModel(Mparam) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  int maxr =  m_Dataij_.maxCoeff();
  int minr =  m_Dataij_.minCoeff();
  r_ = maxr-minr+1;
#ifdef COVERBOSE
  std::cout<<"\nNumber of categories: "<<r_<<"\n";
  std::cout<<"Min category value: "<<minr<<"\n";
  std::cout<<"Max category value: "<<maxr<<"\n";
#endif
  //initialize various data storages
  m3_Yhij_.resize(r_);
  m3_Yijh_.resize(Mparam.nbrowdata_);
  m3_Yjih_.resize(Mparam.nbcoldata_);

  for (int i = 0; i < Mparam.nbrowdata_; ++i){
    m3_Yijh_[i].resize(Mparam_.nbcoldata_,r_);
  }

  for (int j = 0; j < Mparam.nbcoldata_; ++j){
    m3_Yjih_[j].resize(Mparam_.nbrowdata_,r_);
  }

  for (int h = 0; h < r_; ++h) {
    m3_Yhij_[h] = (m_Dataij_.array() == minr+h);
    for (int i = 0; i < Mparam.nbrowdata_; ++i) {
      for (int j = 0; j < Mparam.nbcoldata_; ++j) {
        m3_Yijh_[i](j,h) = m3_Yhij_[h](i,j);
        m3_Yjih_[j](i,h) = m3_Yhij_[h](i,j);
      }
    }
  }
}

CategoricalLBModel::CategoricalLBModel(MatrixInteger const& m_Dataij,VectorInteger const & rowlabels,
                                       VectorInteger const & collabels,ModelParameters const& Mparam
                                       ,int a,int b)
: ICoClustModel(Mparam,rowlabels,collabels) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  int maxr =  m_Dataij_.maxCoeff();
  int minr =  m_Dataij_.minCoeff();
  r_ = maxr-minr+1;
#ifdef COVERBOSE
  std::cout<<"\nNumber of categories: "<<r_<<"\n";
  std::cout<<"Min category value: "<<minr<<"\n";
  std::cout<<"Max category value: "<<maxr<<"\n";
#endif
  //initialize various data storages
  m3_Yhij_.resize(r_);
  m3_Yijh_.resize(Mparam.nbrowdata_);
  m3_Yjih_.resize(Mparam.nbcoldata_);

  for (int i = 0; i < Mparam.nbrowdata_; ++i){
    m3_Yijh_[i].resize(Mparam_.nbcoldata_,r_);
  }

  for (int j = 0; j < Mparam.nbcoldata_; ++j){
    m3_Yjih_[j].resize(Mparam_.nbrowdata_,r_);
  }

  for (int h = 0; h < r_; ++h) {
    m3_Yhij_[h] = (m_Dataij_.array() == minr+h);
    for (int i = 0; i < Mparam.nbrowdata_; ++i) {
      for (int j = 0; j < Mparam.nbcoldata_; ++j) {
        m3_Yijh_[i](j,h) = m3_Yhij_[h](i,j);
        m3_Yjih_[j](i,h) = m3_Yhij_[h](i,j);
      }
    }
  }
}

CategoricalLBModel::~CategoricalLBModel()
{
  // TODO Auto-generated destructor stub
}

void CategoricalLBModel::LogSumRows(MatrixReal & m_sum){
  std::vector<MatrixReal> v_mlogtemphl(Mparam_.nbrowclust_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_mlogtemphl[k].resize(r_,Mparam_.nbcolclust_);
    for (int h = 0; h < r_; ++h) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        v_mlogtemphl[k](h,l) = m3_logAlhphahkl_[h](k,l);
      }
    }
  }

  for (int i = 0; i < Mparam_.nbrowdata_; ++i) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      m_sum(i,k) = v_logPiek_(k) + (m3_Yijh_[i].cast<float>()*v_mlogtemphl[k]).array().sum();
    }
  }
}

void CategoricalLBModel::LogSumCols(MatrixReal & m_sum){
  std::vector<MatrixReal> v_mlogtemphk(Mparam_.nbcolclust_);

  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    v_mlogtemphk[l].resize(r_,Mparam_.nbrowclust_);
    for (int h = 0; h < r_; ++h) {
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        v_mlogtemphk[l](h,k) = m3_logAlhphahkl_[h](k,l);
      }
    }
  }

  for (int j = 0; j < Mparam_.nbcoldata_; ++j) {
    for (int l = 0;l < Mparam_.nbcolclust_; ++l) {
      m_sum(j,l) = v_logRhol_(l) + (m3_Yjih_[j].cast<float>()*v_mlogtemphk[l]).array().sum();
    }
  }
}

void CategoricalLBModel::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=((v_Tk_.array()+a_-1)/(nbSample_+Mparam_.nbrowclust_*(a_-1))).log();
  }

  Array2DReal m_vtkrl = (v_Tk_*v_Rl_.transpose()).array()+r_*(b_-1);
  for (int h = 0; h < r_; ++h) {
    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<float>()*m_Rjl_).array()+b_-1)/(m_vtkrl+RealMin);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }
}

void CategoricalLBModel::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=((v_Rl_.array()+a_-1)/(nbVar_+Mparam_.nbcolclust_*(a_-1))).log();
  }

  Array2DReal m_vtkrl = (v_Tk_*v_Rl_.transpose()).array()+r_*(b_-1);
  for (int h = 0; h < r_; ++h) {
    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<float>()*m_Rjl_).array()+b_-1)/(m_vtkrl+RealMin);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }
}

void CategoricalLBModel::Modify_theta_start()
{
  m3_Alphahklstart_ = m3_Alphahkl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
  m_Tikstart_ = m_Tik_;
}

void CategoricalLBModel::Copy_theta_start()
{
  m3_Alphahkl_ = m3_Alphahklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;
  m_Rjl_ = m_Rjlstart_;
  m_Tik_ = m_Tikstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
}

void CategoricalLBModel::Copy_theta_max()
{
  m3_Alphahkl_ = m3_Alphahklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;
  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

void CategoricalLBModel::Modify_theta_max()
{
  m3_Alphahklmax_ = m3_Alphahkl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;
  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}


bool CategoricalLBModel::EMRows(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ERows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    MStepRows();

    float netchange = 0.0;
    for (int h = 0; h < r_; ++h) {
      netchange+= ((m3_Alphahkl_[h].array()-m3_Alphahklold_[h].array()).abs()/(m3_Alphahkl_[h].array()+RealMin)).sum();
    }

    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) {
      break;
    }
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::CEMRows(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!CERows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    MStepRows();

    float netchange = 0.0;
    for (int h = 0; h < r_; ++h) {
      netchange+= ((m3_Alphahkl_[h].array()-m3_Alphahklold_[h].array()).abs()/(m3_Alphahkl_[h].array()+RealMin)).sum();
    }

    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) {
      break;
    }
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::EMCols(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ECols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    MStepRows();

    float netchange = 0.0;
    for (int h = 0; h < r_; ++h) {
      netchange+= ((m3_Alphahkl_[h].array()-m3_Alphahklold_[h].array()).abs()/(m3_Alphahkl_[h].array()+RealMin)).sum();
    }

    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) {
      break;
    }
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::CEMCols(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!CECols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    MStepRows();

    float netchange = 0.0;
    for (int h = 0; h < r_; ++h) {
      netchange+= ((m3_Alphahkl_[h].array()-m3_Alphahklold_[h].array()).abs()/(m3_Alphahkl_[h].array()+RealMin)).sum();
    }

    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) {
      break;
    }
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::SEMRows(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  if(!SERows()) return false;
  MStepRows();

  return true;

}

bool CategoricalLBModel::SEMCols(){
  //Initializations
  for (int h = 0; h < r_; ++h) {
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

  if(!SECols()) return false;
  MStepRows();

  return true;
}

float CategoricalLBModel::EstimateLikelihood()
{

  Array2DReal m_Ukl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  Array2DReal m_vtkrl = v_Tk_*v_Rl_.transpose();
  float tempsum = -(m_vtkrl*(m_vtkrl+RealMin).log()).sum();
  for (int h = 0; h < r_; ++h) {
    m_Ukl = m_Tik_.transpose()*m3_Yhij_[h].cast<float>()*m_Rjl_;
    tempsum+= (m_Ukl.array()*(m_Ukl+RealMin).log()).sum()+(b_-1)*((m3_Alphahkl_[h].array()+RealMin).log()).sum();
  }
  Likelihood_ = tempsum
               +v_Tk_.dot(v_logPiek_) - Mparam_.nbrowdata_*log(float(Mparam_.nbrowdata_))
               +v_Rl_.dot(v_logRhol_) - Mparam_.nbcoldata_*log(float(Mparam_.nbcoldata_))
               -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
               -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum()
               +(a_-1)*(v_logPiek_.sum()+v_logRhol_.sum());

  return Likelihood_;
}
//Computer Fuzzy clustering criteria and set the terminate variable accordingly
void CategoricalLBModel::likelihoodStopCriteria()
{
  EstimateLikelihood();

  if(std::abs(1-Likelihood_/Likelihood_old)<Mparam_.epsilon_)
      StopAlgo = true;
  else
  {
    Likelihood_old = Likelihood_;
    StopAlgo = false;
  }
}

void CategoricalLBModel::ParameterStopCriteria(){
  float netchange = 0.0;
  for (int h = 0; h < r_; ++h) {
    netchange+= ((m3_Alphahkl1_[h].array()-m3_Alphahkl1old_[h].array()).abs()/(m3_Alphahkl1_[h].array()+RealMin)).sum();
  }

  netchange/=r_;

  if(netchange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

bool CategoricalLBModel::RandomInit(){
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

  v_Tk_ = m_Tik_.colwise().sum().transpose();
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

  //Initializing model parameters
  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);
  m3_Alphahkl1_.resize(r_);
  m3_Alphahklold_.resize(r_);
  m3_Alphahkl1old_.resize(r_);
  Array2DReal m_vtkrl = (v_Tk_*v_Rl_.transpose()).array()+r_*(b_-1);
  for (int h = 0; h < r_; ++h) {
    m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahkl1old_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<float>()*m_Rjl_).array()+b_-1)/(m_vtkrl+RealMin);
    m3_Alphahkl1_[h] = m3_Alphahkl_[h];
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h].array()+RealMin).log();
  }

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

float CategoricalLBModel::ICLCriteriaValue(){
  float criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
      -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
      +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(r_*b_)-r_*lgamma(b_))
      -lgamma(Mparam_.nbrowdata_+Mparam_.nbrowclust_*a_)
      -lgamma(Mparam_.nbcoldata_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    criteria+= lgamma(a_+ (v_Zi_.array()== k).count());
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    criteria+= lgamma(a_+ (v_Wj_.array()==l).count());
  }

  Eigen::ArrayXXi temp(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int h = 0; h < r_; ++h) {
    temp = (m_Zik_.transpose()*m3_Yhij_[h].cast<int>()*m_Wjl_).array()+b_;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        criteria+=lgamma(temp(k,l));
      }
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      criteria-= lgamma(((v_Zi_.array()== k).count())*((v_Wj_.array()==l).count())+r_*b_);
    }
  }

  return criteria;
}

void CategoricalLBModel::FinalizeOutput(){
  CommonFinalizeOutput();
}

void CategoricalLBModel::ConsoleOut(){
#ifndef RPACKAGE
  std::cout<<"Output Model parameters\n";
  std::cout<<"\npie_k:"<<v_Piek_<<"\nrho_l:"<<v_Rhol_<<"\n";
  for (int h = 0; h < r_; ++h) {
    std::cout<<"Alpha_kl for category "<<h<<"\n";
    std::cout<<m3_Alphahkl_[h]<<"\n";
  }
  std::cout<<"\nICL value:"<<ICLCriteriaValue()<<"\n";
#endif
}

const MatrixInteger& CategoricalLBModel::GetArrangedDataClusters()
{
  ArrangedDataCluster<MatrixInteger>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

