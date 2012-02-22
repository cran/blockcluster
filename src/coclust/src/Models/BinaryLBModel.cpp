/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2011  Parmeet Singh Bhatia

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
 * Project:  coclust
 * created on: Nov 25, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file BinaryLBModel.cpp
 *  @brief Implements concrete model class BinaryLBModel derived from ICoClustModel.
 **/

#include <limits.h>
#include <math.h>
#include "BinaryLBModel.h"
#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModel::BinaryLBModel( MatrixBinary const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  nbSample_ = m_Dataij.rows();
  nbVar_ = m_Dataij.cols();
  //std::cout<<nbSample_<<" "<<nbVar_;
  Mparam_.nbrowdata_ = nbSample_;
  Mparam_.nbcoldata_ = nbVar_;
  dimprod_ = nbSample_*nbVar_;
};

bool BinaryLBModel::Estep()
{
  if (EMRows()) {
    if (EMCols()) {
      return true;
    }
  }
  return false;
}

bool BinaryLBModel::CEstep()
{
  if (CEMRows()) {
    if (CEMCols()) {
      return true;
    }
  }
  return false;
}

void BinaryLBModel::Mstep()
{
  // Update Alpha for outer loop
  m_Alphakl1old_ = m_Alphakl1_;
  m_Alphakl1_ = m_Alphakl_;
}

bool BinaryLBModel::FuzzyCEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Fuzzy CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Fuzzy CEM initialization is not valid for this model.";
#endif
  return false;
}
bool BinaryLBModel::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Random initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Random initialization is not valid for this model.";
#endif
  return false;
}

bool BinaryLBModel::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (InitCEMRows()) {
    if (InitCEMCols()) {
      v_Rl_ = m_Rjl_.colwise().sum();
      m_Alphakl1_ = m_Alphakl_;
      m_Alphakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Alphaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_akl_ = MatrixBinary::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Uil_ = MatrixReal::Zero(nbSample_,Mparam_.nbcolclust_);
      m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
      v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
      m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
      v_Zi_ = MatrixInteger::Zero(nbSample_,1);
      v_Wj_ = MatrixInteger::Zero(nbVar_,1);
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
  // Calculate row and column proportions
  if(!Mparam_.fixedproportions_){
    v_Piek_ = v_logPiek_.array().exp();
    v_Rhol_ = v_logRhol_.array().exp();
  }else
  {
    v_Piek_ = (1.0/Mparam_.nbrowclust_)*MatrixReal::Ones(Mparam_.nbrowclust_,1);
    v_Rhol_ = (1.0/Mparam_.nbcolclust_)*MatrixReal::Ones(Mparam_.nbcolclust_,1);;
  }

  // Calculate summary Matrix (Class mean)
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
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
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}


//Compute Bernoulli log-sum for all rows
void BinaryLBModel::BernoulliLogsumRows(MatrixReal & m_sum)
{
  m_sum = m_Uil_*(((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix()).transpose())+
      MatrixReal::Ones(nbSample_,1)*((v_logPiek_+((1-m_Alphakl_.array()+RealMin).log()).matrix()*v_Rl_).transpose());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModel::BernoulliLogsumCols(MatrixReal & m_sum)
{
  m_sum = m_Vjk_*((((m_Alphakl_.array()+RealMin)/(1-m_Alphakl_.array()+RealMin)).log()).matrix())+
      MatrixReal::Ones(nbVar_,1)*((v_logRhol_+(((1-m_Alphakl_.array()+RealMin).log()).matrix()).transpose()*v_Tk_).transpose());
}

//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModel::EMRows()
{
  //Initializations
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;

  //Temporary variables
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      BernoulliLogsumRows(m_sumik);
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
      m_Alphaklold_ = m_Alphakl_;
      m_Alphakl_ = ((m_Tik_.transpose())*m_Uil_).array()/(v_Tk_*(v_Rl_.transpose())).array();

      //Termination check
      if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.epsilon_int_) {
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

//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModel::CEMRows()
{
  //Initializations
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;

  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  //CEMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
        BernoulliLogsumRows(m_sumik);
        m_Tik_.setZero();
        for (int i = 0; i < nbSample_; ++i) {
          m_sumik.row(i).maxCoeff(&maxIndex);
          m_Tik_(i,maxIndex)=1;
        }
        v_Tk_ = m_Tik_.colwise().sum();
        if((v_Tk_.array()<.00001).any()){
          empty_cluster = true;
          throw "Row clustering failed while running model.";
        }
        else
        {empty_cluster = false;}

        //M-step
        if(!Mparam_.fixedproportions_) {
          v_logPiek_=(v_Tk_.array()/nbSample_).log();
        }
        m_Alphaklold_ = m_Alphakl_;
        m_Alphakl_ = ((m_Tik_.transpose())*m_Uil_).array()/(v_Tk_*(v_Rl_.transpose())).array();

        if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.epsilon_int_) {
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

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModel::EMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;

  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_),Wsuml(Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      BernoulliLogsumCols(m_sumjl);
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
      m_Alphaklold_ = m_Alphakl_;
      m_Alphakl_ = (m_Vjk_.transpose()*m_Rjl_).array()/(v_Tk_*v_Rl_.transpose()).array();
      if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.epsilon_int_){
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

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModel::CEMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  //CEMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
        BernoulliLogsumCols(m_sumjl);
        m_Rjl_.setZero();
        for (int j = 0; j < nbVar_; ++j) {
          m_sumjl.row(j).maxCoeff(&maxIndex);
          m_Rjl_(j,maxIndex)=1;
        }
        v_Rl_ = m_Rjl_.colwise().sum();
        if((v_Rl_.array()<.00001).any()){
          empty_cluster = true;
          throw "Column clustering failed while runnin model.";
        }
        else
        {empty_cluster = false;}

        if(!Mparam_.fixedproportions_) {
          v_logRhol_=(v_Rl_.array()/nbVar_).log();
        }

        m_Alphaklold_ = m_Alphakl_;
        m_Alphakl_ = (m_Vjk_.transpose()*m_Rjl_).array()/(v_Tk_*v_Rl_.transpose()).array();

        if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.epsilon_int_){
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

float BinaryLBModel::EstimateLikelihood()
{
	Likelihood_ = (v_Tk_.transpose()*
	          (m_Alphakl_.array()*(m_Alphakl_.array().log()) + (1-m_Alphakl_.array())*((1-m_Alphakl_.array()).log())).matrix()*v_Rl_
	          +v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
	          -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
	          -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum())/dimprod_;

	return Likelihood_;
}
//Computer Fuzzy clustering criteria and set the terminate variable accordingly
void BinaryLBModel::likelihoodStopCriteria()
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

//Compute change in Alpha and set the terminate variable accordingly
void BinaryLBModel::ParameterStopCriteria()
{
  float relativechange = (((m_Alphakl1_.array()-m_Alphakl1old_.array()).abs()/m_Alphakl1_.array()).sum());
  if(relativechange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

MatrixBinary BinaryLBModel::GetArrangedDataClusters()
{
  Eigen::ArrayXi v_Zi = (GetRowClassificationVector()).array();
  Eigen::ArrayXi v_Wj = (GetColumnClassificationVector()).array();
  MatrixBinary m_clusterij = MatrixBinary::Zero(nbSample_,nbVar_);

  //Rearrange data into clusters

  VectorInteger rowincrement = MatrixInteger::Zero(Mparam_.nbrowclust_,1);

  VectorInteger nbindrows = MatrixInteger::Zero(Mparam_.nbrowclust_+1,1);
  for (int k = 1; k < Mparam_.nbrowclust_; ++k) {
    nbindrows(k) = (v_Zi==(k-1)).count()+nbindrows(k-1);
  }

  VectorInteger colincrement = MatrixInteger::Zero(Mparam_.nbcolclust_,1);

  VectorInteger nbindcols = MatrixInteger::Zero(Mparam_.nbcolclust_+1,1);
  for (int l = 1; l < Mparam_.nbcolclust_; ++l) {
    nbindcols(l)=(v_Wj==(l-1)).count()+nbindcols(l-1);
  }

  for (int j = 0; j < nbVar_; ++j) {
    m_clusterij.col(colincrement(v_Wj(j)) + nbindcols(v_Wj(j))) = m_Dataij_.col(j);
    colincrement(v_Wj(j))+=1;
  }
  MatrixBinary temp = m_clusterij;

  for (int i = 0; i < nbSample_; ++i) {
    m_clusterij.row( rowincrement(v_Zi(i)) + nbindrows(v_Zi(i))) = temp.row(i);
    rowincrement(v_Zi(i))+=1;
  }

  /*// calculate mean of each cocluster in image

  //nbindrows.resize(Mparam_.nbrowclust_+1);
  //nbindcols.resize(Mparam_.nbcolclust_+1);
  nbindrows(Mparam_.nbrowclust_) = Mparam_.nbrowdata_;
  nbindcols(Mparam_.nbcolclust_) = Mparam_.nbcoldata_;
  MatrixBinary mean(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int i = 0; i < Mparam_.nbrowclust_; ++i) {
    for (int j = 0; j < Mparam_.nbcolclust_; ++j) {
      MatrixBinary temp = m_clusterij.block(nbindrows(i),nbindcols(j),nbindrows(i+1)-nbindrows(i),nbindcols(j+1)-nbindcols(j));
      //std::cout<<temp.cast<int>().array().sum()<<" "<<temp.rows()<<" "<<temp.cols()<<"\n";
      mean(i,j) = (temp.cast<float>().array().sum()/(temp.rows()*temp.cols()))>=0.5?1:0;
    }
  }
  std::cout<<"cluster mean:"<<mean;*/
  return m_clusterij;

}
#ifndef RPACKAGE
void BinaryLBModel::DisplayCluster()
{
  CImg<unsigned char>  cluster(nbVar_,nbSample_,1,1,0);
  CImg<unsigned char>  data(nbVar_,nbSample_,1,1,0);
  MatrixBinary m_clusterij = GetArrangedDataClusters();

  // Assign value to images
  for (int i = 0; i < nbSample_; ++i) {
    for (int j = 0; j < nbVar_; ++j) {
      if(m_clusterij(i,j) == 0)
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
  Array2DReal m_Alphakltmp = m_Alphakl_.array()+RealMin;
  _m_sum = m_Uil_*(((((m_Alphakltmp)/(1-m_Alphakltmp)).log()).matrix()).transpose())+
      VectorReal::Ones(Mparam_.nbrowdata_)*((((1-m_Alphakltmp).log()).matrix()*v_Rl_).transpose());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModel::InitBernoulliLogsumCols(MatrixReal & _m_sum)
{
  Array2DReal m_Alphakltmp = m_Alphakl_.array()+RealMin;
  _m_sum = m_Vjk_*(((((m_Alphakltmp)/(1-m_Alphakltmp)).log()).matrix()))+
      VectorReal::Ones(Mparam_.nbcoldata_)*(v_Tk_.transpose()*(((1-m_Alphakltmp).log()).matrix()));
}

bool BinaryLBModel::InitCEMRows()
{
  // Initialization of various parameters
  int cols=std::min(100,int(nbVar_));
  m_Uil_ = MatrixReal::Zero(nbSample_,cols);
  SelectRandomColsfromdata(m_Uil_,cols);
  m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Rl_ = MatrixReal::Ones(cols,1);
  v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
  m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
  m_Alphakl_ = MatrixReal::Zero(Mparam_.nbrowclust_,cols);
  GenerateRandomBernoulliParameterRows(m_Alphakl_,cols);
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal _m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index _maxIndex;
  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      InitBernoulliLogsumRows(_m_sumik);
      m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
      for (int i = 0; i < nbSample_; ++i) {
        _m_sumik.row(i).maxCoeff(&_maxIndex);
        m_Tik_(i,_maxIndex)=1;
      }

      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        v_Tk_(k) = (m_Tik_.col(k)).sum();
        if(v_Tk_(k)<.00001){
          empty_cluster = true;
          throw "Row clustering failed while initialization.";
        }
        else
        {empty_cluster = false;}

      }
      // M-step
      m_Alphaklold_ = m_Alphakl_;
      m_Alphakl_ = (m_Tik_.transpose()*m_Uil_).array()/(v_Tk_*(v_Rl_.transpose())).array();

      if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.initepsilon_) {
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

bool BinaryLBModel::InitCEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  VectorReal::Index _maxIndex;
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  m_Alphakl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  GenerateRandomBernoulliParameterCols(m_Alphakl_);

  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      // CE-step
      InitBernoulliLogsumCols(m_sumjl);
      m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
      for (int j = 0; j < nbVar_; ++j) {
        m_sumjl.row(j).maxCoeff(&_maxIndex);
        m_Rjl_(j,_maxIndex)=1;
      }
      v_Rl_.resize(Mparam_.nbcolclust_);
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        v_Rl_(l)= (m_Rjl_.col(l)).sum();
       if(v_Rl_(l)<.00001){
         empty_cluster = true;
         throw "Column cluster failed while initialization.";
       }
       else
       {empty_cluster = false;}

       }

      // M-step
      m_Alphaklold_ = m_Alphakl_;
      m_Alphakl_ = (m_Vjk_.transpose()*m_Rjl_).array()/ (v_Tk_*(v_Rl_.transpose())).array();
      if((((m_Alphakl_.array()-m_Alphaklold_.array()).abs()/m_Alphakl_.array()).sum())<Mparam_.initepsilon_){
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

void BinaryLBModel::SelectRandomColsfromdata(MatrixReal& _m_il,int cols)
{
  if(cols==Mparam_.nbcoldata_)
    _m_il=m_Dataij_.cast<float>();
  else{
    //random shuffle Algorithm
    VectorInteger _v_temp = RandSample(nbVar_,cols);

    for (int l = 0; l < cols; ++l){
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
  for (int k = 0; k < Mparam_.nbrowclust_; ++k){
    index=_v_temp(k);
    //index=k;
    for (int l = 0; l < cols; ++l) {
      _m_kl(k,l)=epsilon*(1.0-m_Uil_(index,l))+m_Uil_(index,l)*(1.0-epsilon);
    }
  }
}

void BinaryLBModel::GenerateRandomBernoulliParameterCols(MatrixReal& _m_kl)
{
  int index;
  VectorInteger _v_temp = RandSample(nbVar_,Mparam_.nbcolclust_);

  for (int l = 0; l < Mparam_.nbcolclust_; ++l){
    index=_v_temp(l);
    //index=l;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      _m_kl(k,l)=m_Vjk_(index,k)/v_Tk_(k);
    }
  }
}
