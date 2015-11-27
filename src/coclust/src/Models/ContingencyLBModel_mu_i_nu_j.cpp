/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

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

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j(MatrixReal const& m_Dataij,VectorReal const& v_Mui,
                                                           VectorReal const& v_Nuj,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij),v_Mui_(v_Mui),v_Nuj_(v_Nuj)
{
  DataSum_ = m_Dataij.sum();
};

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,VectorInteger const & collabels,
                                                           VectorReal const& v_Mui,VectorReal const& v_Nuj,ModelParameters const& Mparam)
                             : ICoClustModel(Mparam,rowlabels,collabels)
                             , m_Dataij_(m_Dataij),v_Mui_(v_Mui),v_Nuj_(v_Nuj)
{
  DataSum_ = m_Dataij.sum();
};

bool ContingencyLBModel_mu_i_nu_j::emRows()
{
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::semRows()
{
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::cemRows()
{
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_nul_ = m_Rjl_.transpose()*v_Nuj_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {

    if(!ceStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::emCols()
{
  //Initialization
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!eStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();

    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::semCols()
{
  //Initialization
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;

  if(!seStepCols()) return false;
  //M-step
  mStepCols();

  return true;
}

bool ContingencyLBModel_mu_i_nu_j::cemCols()
{
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_muk_ = m_Tik_.transpose()*v_Mui_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    if(!ceStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();

    //Termination check
    if((((m_Gammakl_-m_Gammaklold_).abs()/m_Gammakl_).sum())<Mparam_.epsilon_int_) {
      break;
    }
  }

  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

void ContingencyLBModel_mu_i_nu_j::logSumRows(MatrixReal & m_ik)
{
  m_ik = STK::Const::Vector<STK::Real>(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.log()).transpose())
      -v_Mui_*(m_Gammakl_*v_nul_).transpose();
}

void ContingencyLBModel_mu_i_nu_j::logSumCols(MatrixReal & m_jl)
{
  m_jl = STK::Const::Vector<STK::Real>(nbVar_)*v_logRhol_.transpose() + m_Vjk_*((m_Gammakl_.log()))
      -v_Nuj_*(v_muk_.transpose()*m_Gammakl_);
}

STK::Real ContingencyLBModel_mu_i_nu_j::estimateLikelihood()
{
  likelihood_ = (m_Ykl_.prod(m_Gammakl_.log()) ).sum() - DataSum_ + v_Tk_.transpose()*v_logPiek_ + v_Rl_.transpose()*v_logRhol_
            -(m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
            -(m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum();

  return likelihood_;
}

void ContingencyLBModel_mu_i_nu_j::parameterStopCriteria()
{
  STK::Real relativechange = (((m_Gammakl1_-m_Gammakl1old_).abs()/m_Gammakl1_).sum());
  if(relativechange<Mparam_.epsilon_)
    stopAlgo_ = true;
  else
    stopAlgo_ = false;
}

bool ContingencyLBModel_mu_i_nu_j::randomInit()
{
#ifdef COVERBOSE
  std::cout<<"Entering ContingencyLBModel_mu_i_nu_j::randomInit()."<<"\n";
#endif
  //Initialize random row and column partition
  VectorReal probarows = (1.0/Mparam_.nbrowclust_)*STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_);
  VectorReal probacols = (1.0/Mparam_.nbcolclust_)*STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_);
  v_Zi_ = partRnd(nbSample_,probarows);
  v_Wj_ = partRnd(nbVar_,probacols);
  m_Zik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Wjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  std::pair<int,int> Label_pair;
  for(int i =0;i< (int)knownLabelsRows_.size();i++)
  {
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1.0;
  }

  for ( int i =0;i< (int)UnknownLabelsRows_.size();i++)
  { m_Tik_(i,v_Zi_[UnknownLabelsRows_[i]]-1) = 1.0;}

  for ( int j=0;j< (int)knownLabelsCols_.size();j++)
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1.0;
  }
  for ( int j =0;j< (int)UnknownLabelsCols_.size();j++)
  { m_Rjl_(j,v_Wj_[UnknownLabelsCols_[j]]-1) = 1.0;}
  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ContingencyLBModel_mu_i_nu_j::randomInit(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_*m_Rjl_)/(m_Tik_.transpose()*v_Mui_*v_Nuj_.transpose()*m_Rjl_);
  //Initializing model parameters
  m_Gammakl1_ = m_Gammakl_;
  m_Gammakl1old_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
  m_Gammaklold_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
  m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;
  v_Zi_.resize(nbSample_) = 0;
  v_Wj_.resize(nbVar_) = 0;
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_));

#ifdef COVERBOSE
  std::cout<<"Terminating ContingencyLBModel_mu_i_nu_j::randomInit()."<<"\n";
#endif

  return true;
}

void ContingencyLBModel_mu_i_nu_j::finalizeOutput()
{ commonFinalizeOutput();}

void ContingencyLBModel_mu_i_nu_j::consoleOut()
{
#ifndef COVERBOSE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

MatrixReal const& ContingencyLBModel_mu_i_nu_j::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContingencyLBModel_mu_i_nu_j::modifyThetaStart()
{
  m_Gammaklstart_ = m_Gammakl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

void ContingencyLBModel_mu_i_nu_j::copyThetaStart()
{
  m_Gammakl_ = m_Gammaklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = STK::sum(m_Rjl_);
  m_Gammakl1_ = m_Gammakl_;
}

void ContingencyLBModel_mu_i_nu_j::copyThetaMax()
{
  m_Gammakl_ = m_Gammaklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  likelihood_ = Lmax_;
}

void ContingencyLBModel_mu_i_nu_j::modifyThetaMax()
{
  m_Gammaklmax_ = m_Gammakl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = likelihood_;
}

void ContingencyLBModel_mu_i_nu_j::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/nbVar_).log();
    v_logPiek_=(v_Tk_/nbSample_).log();
  }
  m_Ykl_     = m_Tik_.transpose()*m_Dataij_*m_Rjl_;
  m_Gammakl_ = m_Ykl_/((m_Tik_.transpose()*v_Mui_)*(v_Nuj_.transpose()*m_Rjl_));
}
