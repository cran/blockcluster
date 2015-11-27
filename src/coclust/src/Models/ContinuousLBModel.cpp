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

void ContinuousLBModel::logSumRows(MatrixReal& m_sumik)
{
//  VectorReal tmp1(Mparam_.nbrowclust_);
//  tmp1 = ((0.5*((m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_)*v_Rl_)) - v_logPiek_);
  m_sumik = STK::Const::Vector<STK::Real>(nbSample_)
            * ( ( v_logPiek_ - 0.5*( (m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_)*v_Rl_))).transpose()
            - 0.5     * (m_Uil2_ * m_Sigma2kl_.inverse().transpose())
            + m_Uil1_ * (m_Mukl_/m_Sigma2kl_).transpose();
}

void ContinuousLBModel::logSumCols(MatrixReal & m_sumjl)
{
//  VectorReal tmp1(Mparam_.nbcolclust_);
//  tmp1 = ((0.5*(v_Tk_.transpose()*(m_Sigma2kl_.log()+m_Mukl_.square()/m_Sigma2kl_))) - v_logRhol_);
  m_sumjl = STK::Const::Vector<STK::Real>(nbVar_)
            * ( v_logRhol_.transpose() - 0.5 * v_Tk_.transpose() * (m_Sigma2kl_.log() + m_Mukl_.square()/m_Sigma2kl_))
          - 0.5    * (m_Vjk2_ * m_Sigma2kl_.inverse())
          + m_Vjk1_* (m_Mukl_/m_Sigma2kl_);
}

bool ContinuousLBModel::emRows()
{
  //Initialization
  //m_Mukl2_ = m_Mukl_.square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepRows();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContinuousLBModel::semRows()
{
  //Initialization
  //m_Mukl2_ = m_Mukl_.square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

//  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContinuousLBModel::emCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;
  //  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepCols();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_) {
      break;
    }
  }
  return true;
}

bool ContinuousLBModel::semCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

//  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}


bool ContinuousLBModel::cemRows()
{
  //Initialization
 //m_Mukl2_ = m_Mukl_.square();
  m_Uil1_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;

//  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepRows()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepRows();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}


bool ContinuousLBModel::cemCols()
{
  //Initialization
  m_Vjk1_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;

//  MatrixReal m_trkl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepCols()) return false;
    //M-step
    m_Muklold2_ = m_Mukl_;
    mStepCols();
    //Termination check
    if((((m_Mukl_-m_Muklold2_)/m_Mukl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

void ContinuousLBModel::parameterStopCriteria()
{
  if((((m_Mukl_-m_Muklold1_)/m_Mukl_).abs().sum())<Mparam_.epsilon_)
  {
    stopAlgo_ = true;
  }
  else
  {
    stopAlgo_ = false;
    m_Muklold1_ =m_Mukl_;
  }
}

STK::Real ContinuousLBModel::estimateLikelihood()
{
  likelihood_ = (-0.5*(dimprod_+v_Tk_.dot(m_Sigma2kl_.log()*v_Rl_+v_logPiek_))
//              + v_Tk_.dot(v_logPiek_)
              + v_Rl_.dot(v_logRhol_)
              - (m_Tik_.prod( (RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod( (RealMin + m_Rjl_).log()) ).sum())/dimprod_;
  return likelihood_;
}

// Initialization using CEM algo
bool ContinuousLBModel::cemInitStep()
{
#ifdef COVERBOSE
  std::cout<<"Running Initialization.."<<"\n";
#endif
  if (initCEMCols())
  {
    v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_));
    v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_));
    m_Dataij2_ = m_Dataij_.square();
    m_Uil1_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
    m_Uil2_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
    m_Vjk1_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
    m_Vjk2_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
    v_Tk_.resize(Mparam_.nbrowclust_) = 0;
    m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
    m_Muklold1_ = m_Mukl_;
    //m_Mukl2_ = m_Mukl_.square();
    m_Muklold2_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
#ifdef COVERBOSE
    std::cout<<"Initialization over."<<"\n";
#endif
    return true;
  }
    return false;
}

// Private initialization  functions

bool ContinuousLBModel::initCEMCols()
{
  //Temporary variables
  MatrixReal m_Vjk(nbVar_,Mparam_.nbrowclust_);
  MatrixReal m_Vkj(Mparam_.nbrowclust_,nbVar_);
  MatrixReal m_Djl(nbVar_,Mparam_.nbcolclust_);
  MatrixReal m_Mulk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
 // PointReal  p_Mul2(Mparam_.nbcolclust_);
  MatrixReal m_Rlk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  int minIndex;
  STK::Real W1 = RealMax , W1_old;

  //Private members initialization
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_);
  m_Rjl_ = 0;
  std::pair<int,int> Label_pair;
  for ( int j=0;j<knownLabelsCols_.size();j++)
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
  m_Sigma2kl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
  //Initializations
  selectRandomRowsFromData(m_Vkj);
  m_Vjk = m_Vkj.transpose();
  VectorReal v_Vj2 = STK::sumByRow(m_Vjk.square());
  generateRandomMean(m_Vjk,m_Mulk);

  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    m_Djl  = v_Vj2*STK::Const::Point<STK::Real>(Mparam_.nbcolclust_)
           + STK::Const::Vector<STK::Real>(nbVar_)*STK::sumByRow(m_Mulk.square()).transpose()
           - 2*(m_Vjk*m_Mulk.transpose());
    for ( int j =0;j< UnknownLabelsCols_.size();j++)
    {
      m_Djl.row(UnknownLabelsCols_[j]).minElt(minIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
      m_Rjl_(UnknownLabelsCols_[j],minIndex)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepCols()) )
    {
      Error_msg_  = "In ContinuousLBModel::initCEMCols(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
#endif
      return false;
    }
    //Check for termination
    W1_old = W1;
    W1 = (m_Rjl_.prod(m_Djl)).sum();
    //Initialization of model parameters
    m_Mukl_ = ( (m_Rjl_.transpose()*m_Vjk)
               /(v_Rl_*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_))).transpose();
    m_Sigma2kl_ = W1/(Mparam_.nbrowclust_*Mparam_.nbcolclust_);
    if (std::abs((W1-W1_old)/W1)<Mparam_.initepsilon_)
    { break;}
  }

  return true;
}

void ContinuousLBModel::finalizeOutput()
{
  commonFinalizeOutput();
}

void ContinuousLBModel::consoleOut()
{
#ifndef COVERBOSE
    std::cout<<"Output Model parameter:"<<"\nBlock Mean:\n"<<m_Mukl_<<"\nBlock Sigma:\n"<<m_Sigma2kl_<<"\npiek: "<<
        v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}

void ContinuousLBModel::selectRandomRowsFromData(MatrixReal& _m_kj)
{
  VectorInteger v_temp= randSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    //_m_kj.row(k) = m_Dataij_.row(v_temp(k));
    _m_kj.row(k) = m_Dataij_.row(k);
  }
}

void ContinuousLBModel::generateRandomMean(const MatrixReal & m_jk, MatrixReal & mean_lk)
{
  VectorInteger v_temp= randSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    //mean_lk.row(l) = m_jk.row(v_temp[l]);
    mean_lk.row(l) = m_jk.row(l);
  }
}

MatrixReal const& ContinuousLBModel::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContinuousLBModel::modifyThetaStart()
{
  m_Muklstart_ = m_Mukl_;
  m_Sigma2klstart_ = m_Sigma2kl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

void ContinuousLBModel::copyThetaStart()
{
  m_Mukl_ = m_Muklstart_;
  m_Sigma2kl_ = m_Sigma2klstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  // initializations
  //m_Mukl2_ = m_Mukl_.square();
  m_Muklold1_ = m_Mukl_;
  v_Rl_ = STK::sum(m_Rjl_);
}

void ContinuousLBModel::copyThetaMax()
{
  m_Mukl_ = m_Muklmax_;
  m_Sigma2kl_ = m_Sigma2klmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  likelihood_ = Lmax_;
}

void ContinuousLBModel::modifyThetaMax()
{
  m_Muklmax_ = m_Mukl_;
  m_Sigma2klmax_ = m_Sigma2kl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = likelihood_;
}

void ContinuousLBModel::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/nbVar_).log();
    v_logPiek_=(v_Tk_/nbSample_).log();
  }

  MatrixReal m_trkl = (v_Tk_*v_Rl_.transpose());
  if (nbSample_<nbVar_)
  { m_Mukl_ = ((m_Tik_.transpose()*m_Dataij_)*m_Rjl_)/m_trkl;}
  else
  { m_Mukl_ = ( m_Tik_.transpose()*(m_Dataij_*m_Rjl_))/m_trkl;}
  //m_Mukl2_    = m_Mukl_.square();
  if (nbSample_<nbVar_)
  { m_Sigma2kl_ = ( (m_Tik_.transpose()*m_Dataij2_)*m_Rjl_)/m_trkl - m_Mukl_.square();}
  else
  { m_Sigma2kl_ = ( m_Tik_.transpose()*(m_Dataij2_*m_Rjl_))/m_trkl - m_Mukl_.square();}
}

