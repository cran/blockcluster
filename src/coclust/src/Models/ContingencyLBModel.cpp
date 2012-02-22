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

/** @file ContingencyLBModel.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#include "ContingencyLBModel.h"

ContingencyLBModel::ContingencyLBModel(MatrixInteger const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  nbSample_ = m_Dataij.rows();
  nbVar_ = m_Dataij.cols();
  Mparam_.nbrowdata_ = nbSample_;
  Mparam_.nbcoldata_ = nbVar_;
  DataSum_ = m_Dataij.array().sum();
};

bool ContingencyLBModel::Estep()
{
  if (EMRows()) {
    if (EMCols()) {
      return true;
    }
  }
  return false;
}

void ContingencyLBModel::Mstep()
{
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
}

bool ContingencyLBModel::CEstep()
{
  if (CEMRows()) {
    if (CEMCols()) {
      return true;
    }
  }
  return false;
}

bool ContingencyLBModel::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (InitCEMRows()) {
    if (InitCEMCols()) {
      v_Rl_ = m_Rjl_.colwise().sum().transpose();
      m_Gammakl1_ = m_Gammakl_;
      m_Gammakl1old_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
      m_Gammaklold_ =MatrixReal::Zero(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
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

bool ContingencyLBModel::FuzzyCEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Fuzzy CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Fuzzy CEM initialization is not valid for this model.";
#endif
  return false;
}
bool ContingencyLBModel::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Random initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Random initialization is not valid for this model.";
#endif
  return false;
}

bool ContingencyLBModel::EMRows()
{
  //Initialization
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_Yl_ = m_Uil_.colwise().sum().transpose();

  //Temporary Variables
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
      m_Gammakl_ = m_Ykl_.array()/(m_Ykl_.rowwise().sum()*v_Yl_.transpose()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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


bool ContingencyLBModel::CEMRows()
{
  m_Uil_ = m_Dataij_.cast<float>()*m_Rjl_;
  v_Yl_ = m_Uil_.colwise().sum().transpose();

  //Temporary Variables
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
      m_Gammakl_ = m_Ykl_.array()/(m_Ykl_.rowwise().sum()*v_Yl_.transpose()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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

bool ContingencyLBModel::EMCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_Yk_ = m_Vjk_.colwise().sum().transpose();

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
      m_Gammakl_ = m_Ykl_.array()/(v_Yk_*m_Ykl_.colwise().sum()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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

bool ContingencyLBModel::CEMCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.cast<float>().transpose()*m_Tik_;
  v_Yk_ = m_Vjk_.colwise().sum().transpose();

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
      }else{empty_cluster = false;}

      if(!Mparam_.fixedproportions_) {
        v_logRhol_=(v_Rl_.array()/nbVar_).log();
      }
      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Gammaklold_ = m_Gammakl_;
      m_Gammakl_ = m_Ykl_.array()/(v_Yk_*m_Ykl_.colwise().sum()).array();

      if((((m_Gammakl_.array()-m_Gammaklold_.array())/m_Gammakl_.array()).abs().sum())<Mparam_.epsilon_int_) {
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
    for (int j = 0; j < Mparam_.nbcoldata_; ++j) {
      _v_temp(j)=j;
    }
    for (int l = 0; l < cols; ++l){
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
  MatrixReal _m_sumik(Mparam_.nbrowdata_,Mparam_.nbrowclust_);
  VectorReal::Index _maxIndex;

  //Model parameters
  m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Rl_ = MatrixReal::Ones(cols,1);
  v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
  m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);

  // Initializations

  m_Uil_ = MatrixReal::Zero(nbSample_,cols);
  SelectRandomColsfromdata(m_Uil_,cols);
  v_Ui_ = m_Uil_.rowwise().sum();
  GenerateRandomPoissonParameterRows(m_Akl,cols);
  //Determine row partition using CEM algorithm with equal proportions
  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      _m_sumik = m_Uil_*(m_Akl.transpose());
      m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
      for (int i = 0; i < nbSample_; ++i) {
        _m_sumik.row(i).maxCoeff(&_maxIndex);
        m_Tik_(i,_maxIndex)=1;
      }
        v_Tk_ = m_Tik_.colwise().sum();
        if((v_Tk_.array()<.00001).any()){
          empty_cluster = true;
          throw "Row clustering failed  while initialization.";
        }else{empty_cluster = false;}

      // M-step
      m_Aklold = m_Akl;
      m_Akl = ((m_Tik_.transpose()*m_Uil_).array()/((m_Tik_.transpose()*v_Ui_)*MatrixReal::Ones(1,cols)).array()+RealMin).log();
      if((((m_Akl.array()-m_Aklold.array()).abs()/m_Akl.array()).sum())<Mparam_.initepsilon_) {
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

bool ContingencyLBModel::InitCEMCols()
{
  //Temporary variables
  MatrixReal m_Alk(Mparam_.nbcolclust_,Mparam_.nbrowclust_) , m_Alkold(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  VectorReal::Index _maxIndex;

  //Initializations
  m_Vjk_ = (m_Dataij_.cast<float>().transpose())*m_Tik_;
  v_Vj_ = m_Vjk_.rowwise().sum();
  GenerateRandomPoissonParameterCols(m_Alk);

  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      // CE-step
      m_sumjl = m_Vjk_*(m_Alk.transpose());
      m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
      for (int j = 0; j < nbVar_; ++j) {
        m_sumjl.row(j).maxCoeff(&_maxIndex);
        m_Rjl_(j,_maxIndex)=1;
      }
      v_Rl_ = m_Rjl_.colwise().sum();
      if((v_Rl_.array()<.00001).any()){
        empty_cluster = true;
        throw "Column clustering failed while initialization.";
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

void ContingencyLBModel::PoissonLogsumRows(MatrixReal & m_ik)
{
  m_ik = VectorReal::Ones(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.array().log()).matrix().transpose());
}

void ContingencyLBModel::PoissonLogsumCols(MatrixReal & m_jl)
{
  m_jl = VectorReal::Ones(nbVar_)*v_logRhol_.transpose() +
      m_Vjk_*((m_Gammakl_.array().log()).matrix());
}

void ContingencyLBModel::GenerateRandomPoissonParameterRows(MatrixReal& m_kl,int cols)
{
  int index;
  VectorInteger _v_temp = RandSample(nbSample_,Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k){
    index=_v_temp(k);
    //index=k;
    for (int l = 0; l < cols; ++l) {
      m_kl(k,l)=std::log(m_Uil_(index,l)/v_Ui_(index)+RealMin);
    }
  }
}

void ContingencyLBModel::GenerateRandomPoissonParameterCols(MatrixReal& m_lk)
{
  int index;
  VectorInteger _v_temp = RandSample(nbVar_,Mparam_.nbcolclust_);
  for (int l = 0; l < Mparam_.nbcolclust_; ++l){
    index=_v_temp(l);
    //index=l;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      m_lk(l,k) = std::log(m_Vjk_(index,k)/v_Vj_(index)+RealMin);
    }
  }
}

void ContingencyLBModel::FinalizeOutput()
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

void ContingencyLBModel::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

MatrixInteger ContingencyLBModel::GetArrangedDataClusters()
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
