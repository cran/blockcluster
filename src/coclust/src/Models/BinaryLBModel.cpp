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


/** @file BinaryLBModel.cpp
 *  @brief Implements concrete model class BinaryLBModel derived from ICoClustModel.
 **/

#include <limits.h>
#include <math.h>
#include "BinaryLBModel.h"
#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModel::BinaryLBModel( MatrixBinary const& m_Dataij,ModelParameters const& Mparam,int a,int b)
                           : ICoClustModel(Mparam) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  dimprod_ = nbSample_*nbVar_;
};

BinaryLBModel::BinaryLBModel(MatrixBinary const& m_Dataij,VectorInteger const & rowlabels,
                             VectorInteger const & collabels,ModelParameters const& Mparam,int a,int b)
                            : ICoClustModel(Mparam,rowlabels,collabels) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  dimprod_ = nbSample_*nbVar_;
}


bool BinaryLBModel::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (InitCEMRows()) {
    if (InitCEMCols()) {
      m_Alphakl1_ = m_Alphakl_;
      m_Alphakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Alphaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_akl_ = MatrixBinary::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
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

void BinaryLBModel::FinalizeOutput()
{
  CommonFinalizeOutput();
  // Calculate summary Matrix (Class mean)
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l) {
      if(m_Alphakl_(k,l)>=.5)
        m_akl_(k,l) = 1;
      else {
        m_akl_(k,l) = 0;
      }
    }
  }
  // Calculate probability for Summary Matrix
  m_epsilonkl_ = m_akl_.cast<float>().array()*(1-m_Alphakl_.array()) + (1-m_akl_.cast<float>().array())*m_Alphakl_.array();

}

void BinaryLBModel::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\nakl:\n"<<m_akl_<<"\nepsilonkl:\n"<<m_epsilonkl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<"\n";

  std::cout<<"ICL: "<<ICLCriteriaValue()<<"\n";
#endif
}


//Compute Bernoulli log-sum for all rows
void BinaryLBModel::LogSumRows(MatrixReal & m_sum)
{
  m_sum = m_Uil_*(((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix()).transpose())+
      MatrixReal::Ones(nbSample_,1)*((v_logPiek_+((1-m_Alphakl_.array()+RealMin).log()).matrix()*v_Rl_).transpose());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModel::LogSumCols(MatrixReal & m_sum)
{
  m_sum = m_Vjk_*((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix())+
      MatrixReal::Ones(nbVar_,1)*((v_logRhol_+(((1-m_Alphakl_.array()+RealMin).log()).matrix()).transpose()*v_Tk_).transpose());
}

//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModel::EMRows()
{
  //Initializations
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ERows()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    MStepRows();
    //Termination check
    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}


//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModel::CEMRows()
{
  //Initializations
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    //CE-step
    if(!CERows()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    MStepRows();

    //Termination check
    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModel::EMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    if(!ECols()) return false;
     //M-step
    m_Alphaklold_ = m_Alphakl_;
    MStepCols();

    //Termination check
    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.epsilon_int_){
      break;
    }
  }

  // Update Alpha for outer loop
  m_Alphakl1old_ = m_Alphakl1_;
  m_Alphakl1_ = m_Alphakl_;
  return true;
}

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModel::CEMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //CE-step
    if(!CECols()) return false;
    //M-step
    m_Alphaklold_ = m_Alphakl_;
    MStepCols();
    //Termination check
    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.epsilon_int_){
      break;
    }
  }
  // Update Alpha for outer loop
  m_Alphakl1old_ = m_Alphakl1_;
  m_Alphakl1_ = m_Alphakl_;
  return true;
}

bool BinaryLBModel::SEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;

  if(!SERows()) return false;

  //M-step : update row proportions and model parameters
  MStepRows();
  return true;
}

bool BinaryLBModel::SEMCols()
{
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;

  if(!SECols()) return false;

  //M-step: Update column proportions and model parameters
  MStepCols();

  return true;
}

float BinaryLBModel::EstimateLikelihood()
{
	Likelihood_ = (v_Tk_.transpose()*
	          (m_Alphakl_.array()*((m_Alphakl_.array()+RealMin).log()) + (1-m_Alphakl_.array())*((1-m_Alphakl_.array()+RealMin).log())).matrix()*v_Rl_
	          +v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
	          -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
	          -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum())/dimprod_;

	return Likelihood_;
}
//Computer Fuzzy clustering criteria and set the terminate variable accordingly
void BinaryLBModel::likelihoodStopCriteria()
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

//Compute change in Alpha and set the terminate variable accordingly
void BinaryLBModel::ParameterStopCriteria()
{
  float relativechange = (((m_Alphakl1_.array()-m_Alphakl1old_.array()).abs()/(m_Alphakl1_.array()+RealMin)).sum());
  if(relativechange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

const MatrixBinary& BinaryLBModel::GetArrangedDataClusters()
{
  ArrangedDataCluster<MatrixBinary>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}
#ifndef RPACKAGE
void BinaryLBModel::DisplayCluster()
{
  CImg<unsigned char>  cluster(nbVar_,nbSample_,1,1,0);
  CImg<unsigned char>  data(nbVar_,nbSample_,1,1,0);
  MatrixBinary m_ClusterDataij_ = GetArrangedDataClusters();

  // Assign value to images
  for ( int i = 0; i < nbSample_; ++i) {
    for ( int j = 0; j < nbVar_; ++j) {
      if(m_ClusterDataij_(i,j) == 0)
        cluster(j,i) = (unsigned char)(0);
      else
        cluster(j,i) = (unsigned char)(255);

      if(m_Dataij_(i,j) == 0)
        data(j,i) = (unsigned char)(0);
      else
        data(j,i) = (unsigned char)(255);

    }
  }

  //Display data and cluster
  CImgDisplay data_disp(data,"Original Data"), cluster_disp(cluster,"Co-Clusters");
  while (!data_disp.is_closed() && !cluster_disp.is_closed()) {
    data_disp.wait();
    cluster_disp.wait();
  }

  data.save("data.jpg");
  cluster.save("cluster.jpg");
}
#endif
//Compute Bernoulli log-sum for all rows
void BinaryLBModel::InitBernoulliLogsumRows(MatrixReal & _m_sum)
{
  //Array2DReal m_Alphakltmp = m_Alphakl_.array()+RealMin;
  _m_sum = m_Uil_*(((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix()).transpose())+
      VectorReal::Ones(Mparam_.nbrowdata_)*((((1-m_Alphakl_.array()+RealMin).log()).matrix()*v_Rl_).transpose());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModel::InitBernoulliLogsumCols(MatrixReal & _m_sum)
{
  _m_sum = m_Vjk_*(((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix()))+
      VectorReal::Ones(Mparam_.nbcoldata_)*(v_Tk_.transpose()*(((1-m_Alphakl_.array()+RealMin).log()).matrix()));
}

bool BinaryLBModel::InitCEMRows()
{
  // Initialization of various parameters
  int cols=std::min(100,int(nbVar_));
  m_Uil_ = MatrixReal::Zero(nbSample_,cols);
  SelectRandomColsfromdata(m_Uil_,cols);
  m_Vjk_ = MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Rl_ = VectorReal::Ones(cols);
  v_Tk_ = VectorReal::Zero(Mparam_.nbrowclust_);
  m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
  m_Alphakl_ = MatrixReal::Zero(Mparam_.nbrowclust_,cols);
  GenerateRandomBernoulliParameterRows(m_Alphakl_,cols);
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
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

  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
    InitBernoulliLogsumRows(m_sumik);
#ifdef RANGEBASEDFORLOOP
    for ( int i : UnknownLabelsRows_) {
      m_sumik.row(i).maxCoeff(&maxIndex);
      m_Tik_.row(i).setZero();
      m_Tik_(i,maxIndex)=1;
    }
#else
    for ( int i=0; i<UnknownLabelsRows_.size();i++) {
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
    }
    else
    {empty_cluster = false;}

    // M-step
    m_Alphaklold_ = m_Alphakl_;
    m_Alphakl_ = (m_Tik_.transpose()*m_Uil_).array()/(v_Tk_*(v_Rl_.transpose())).array();

    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.initepsilon_) {
      break;
    }
  }
  return true;
}

bool BinaryLBModel::InitCEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  m_Alphakl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  v_Rl_.resize(Mparam_.nbcolclust_);
  GenerateRandomBernoulliParameterCols(m_Alphakl_);
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
    InitBernoulliLogsumCols(m_sumjl);
#ifdef RANGEBASEDFORLOOP
    for ( int j : UnknownLabelsCols_) {
      m_sumjl.row(j).maxCoeff(&maxIndex);
      m_Rjl_.row(j).setZero();
      m_Rjl_(j,maxIndex)=1;
    }
#else
    for ( int j=0;j< UnknownLabelsCols_.size();j++) {
      m_sumjl.row(UnknownLabelsCols_[j]).maxCoeff(&maxIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZero();
      m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
    }
#endif

    //check for empty cluster
    v_Rl_ = m_Rjl_.colwise().sum();
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

    // M-step
    m_Alphaklold_ = m_Alphakl_;
    m_Alphakl_ = (m_Vjk_.transpose()*m_Rjl_).array()/ (v_Tk_*(v_Rl_.transpose())).array();
    if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/(m_Alphakl_.array()+RealMin)).sum())<Mparam_.initepsilon_){
      break;
    }
  }
  return true;
}

void BinaryLBModel::SelectRandomColsfromdata(MatrixReal& _m_il,int cols)
{
  if(cols==Mparam_.nbcoldata_)
    _m_il=m_Dataij_.cast<float>();
  else{
    //random shuffle Algorithm
    VectorInteger _v_temp = RandSample(nbVar_,cols);

    for ( int l = 0; l < cols; ++l){
      _m_il.col(l)=m_Dataij_.cast<float>().col(_v_temp(l));
      //_m_il.col(l)=m_Dataij_.cast<float>().col(l);
    }
  }
}

void BinaryLBModel::GenerateRandomBernoulliParameterRows(MatrixReal& _m_kl,int cols)
{
  int index;
  float epsilon = 0.1;
  VectorInteger _v_temp = RandSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k){
    index=_v_temp(k);
    //index=k;
    for ( int l = 0; l < cols; ++l) {
      _m_kl(k,l)=epsilon*(1.0-m_Uil_(index,l))+m_Uil_(index,l)*(1.0-epsilon);
    }
  }
}

void BinaryLBModel::GenerateRandomBernoulliParameterCols(MatrixReal& _m_kl)
{
  int index;
  VectorInteger _v_temp = RandSample(nbVar_,Mparam_.nbcolclust_);

  for ( int l = 0; l < Mparam_.nbcolclust_; ++l){
    index=_v_temp(l);
    //index=l;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
      _m_kl(k,l)=m_Vjk_(index,k)/v_Tk_(k);
    }
  }
}

void BinaryLBModel::Modify_theta_start()
{
  m_Alphaklstart_ = m_Alphakl_;
  m_epsilonklstart_ = m_epsilonkl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
  m_Tikstart_ = m_Tik_;
}

void BinaryLBModel::Copy_theta_start()
{
  m_Alphakl_ = m_Alphaklstart_;
  m_epsilonkl_ = m_epsilonklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;
  m_Rjl_ = m_Rjlstart_;
  m_Tik_ = m_Tikstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
  m_Alphakl1_ = m_Alphakl_;
}

void BinaryLBModel::Copy_theta_max()
{
  m_Alphakl_ = m_Alphaklmax_;
  m_epsilonkl_ = m_epsilonklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;
  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

void BinaryLBModel::Modify_theta_max()
{
  m_Alphaklmax_ = m_Alphakl_;
  m_epsilonklmax_ = m_epsilonkl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;
  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}


void BinaryLBModel::MStepFull()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  m_Alphakl_ = (m_Tik_.transpose()*m_Dataij_.cast<float>()*m_Rjl_).array()/(v_Tk_*v_Rl_.transpose()).array();
}


void BinaryLBModel::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=((v_Tk_.array()+a_-1)/(nbSample_+Mparam_.nbrowclust_*(a_-1))).log();
  }

  m_Alphakl_ = (((m_Tik_.transpose())*m_Uil_).array()+b_-1)/((v_Tk_*(v_Rl_.transpose())).array()+2*(b_-1));
}

void BinaryLBModel::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=((v_Rl_.array()+a_-1)/(nbVar_+Mparam_.nbcolclust_*(a_-1))).log();
  }

  m_Alphakl_ = ((m_Vjk_.transpose()*m_Rjl_).array()+b_-1)/((v_Tk_*v_Rl_.transpose()).array()+2*(b_-1));
}


float BinaryLBModel::ICLCriteriaValue(){
  float criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
      -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
      +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(2*b_)-2*lgamma(b_))
      -lgamma(Mparam_.nbrowdata_+Mparam_.nbrowclust_*a_)
      -lgamma(Mparam_.nbcoldata_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    criteria+= lgamma(a_+ (v_Zi_.array()== k).count());
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    criteria+= lgamma(a_+ (v_Wj_.array()==l).count());
  }

  Eigen::ArrayXXi temp0(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  Eigen::ArrayXXi temp1(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixBinary m_tempdata = (m_Dataij_.array()==0);
  temp0 = (m_Zik_.transpose()*m_tempdata.cast<int>()*m_Wjl_).array()+b_;
  temp1 = (m_Zik_.transpose()*m_Dataij_.cast<int>()*m_Wjl_).array()+b_;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      criteria+=lgamma(temp0(k,l))+lgamma(temp1(k,l));
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      criteria-= lgamma(((v_Zi_.array()== k).count())*((v_Wj_.array()==l).count())+2*b_);
    }
  }

  return criteria;
}
