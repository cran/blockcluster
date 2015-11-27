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

/** @file ContingencyLBModel.cpp
 *  @brief Implements concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#include "ContingencyLBModel.h"

ContingencyLBModel::ContingencyLBModel(MatrixReal const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  DataSum_ = m_Dataij.sum();
};

ContingencyLBModel::ContingencyLBModel(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,
                                       VectorInteger const & collabels,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam,rowlabels,collabels)
                           , m_Dataij_(m_Dataij)
{
  DataSum_ = m_Dataij.sum();
};


bool ContingencyLBModel::cemInitStep()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (initCEMRows())
  {
    if (initCEMCols())
    {
      m_Gammakl1_ = m_Gammakl_;
      m_Gammakl1old_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
      m_Gammaklold_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
      m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
      m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_));
#ifdef COVERBOSE
  std::cout<<"Initialization over."<<std::endl;
#endif
      return true;
    }
  }
  return false;
}

bool ContingencyLBModel::emRows()
{
  //Initialization
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_Yl_ = STK::sum(m_Uil_).transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();

    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContingencyLBModel::semRows()
{
  //Initialization
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_Yl_ = STK::sum(m_Uil_).transpose();

  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool ContingencyLBModel::cemRows()
{
  m_Uil_ = m_Dataij_*m_Rjl_;
  v_Yl_ = STK::sum(m_Uil_).transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepRows()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepRows();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  return true;
}

bool ContingencyLBModel::emCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_Yk_  = STK::sum(m_Vjk_).transpose();
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!eStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

bool ContingencyLBModel::semCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_Yk_  = STK::sum(m_Vjk_).transpose();

  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool ContingencyLBModel::cemCols()
{
  //Initializations
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  v_Yk_  = STK::sum(m_Vjk_).transpose();

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    if(!ceStepCols()) return false;
    //M-step
    m_Gammaklold_ = m_Gammakl_;
    mStepCols();
    //Termination check
    if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.epsilon_int_)
    { break;}
  }
  m_Gammakl1old_ = m_Gammakl1_;
  m_Gammakl1_ = m_Gammakl_;
  return true;
}

STK::Real ContingencyLBModel::estimateLikelihood()
{
  likelihood_ = (m_Ykl_.prod((m_Gammakl_).log()) ).sum()
              - DataSum_ + v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
              -(m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              -(m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum();
  return likelihood_;
}

void ContingencyLBModel::parameterStopCriteria()
{
  STK::Real relativechange = (((m_Gammakl1_-m_Gammakl1old_)/m_Gammakl1_).abs().sum());
  if(relativechange<Mparam_.epsilon_)
    stopAlgo_ = true;
  else
    stopAlgo_ = false;
}

void ContingencyLBModel::selectRandomColsFromData(MatrixReal& _m_il,int cols)
{
  if(cols==nbVar_)
    _m_il = m_Dataij_;
  else{
    //random shuffle Algorithm
    int random,index,temp;
    VectorInteger _v_temp(Mparam_.nbcoldata_);
    for ( int j = 0; j < Mparam_.nbcoldata_; ++j)
    { _v_temp[j]=j;}
    for ( int l = 0; l < cols; ++l)
    {
      random=std::rand()%(Mparam_.nbcoldata_-l);
      index=_v_temp[random];
      _m_il.col(l)=m_Dataij_.col(index);
      //swap elements
      temp=_v_temp[Mparam_.nbcoldata_-l-1];
      _v_temp[Mparam_.nbcoldata_-l-1]=_v_temp[random];
      _v_temp[random]=temp;
    }
  }
}

bool ContingencyLBModel::initCEMRows()
{
  //Temporary variables
  int cols=std::min(100,int(nbVar_));
  MatrixReal m_Akl(Mparam_.nbrowclust_,cols),m_Aklold(Mparam_.nbrowclust_,cols);
  MatrixReal m_sumik(Mparam_.nbrowdata_,Mparam_.nbrowclust_);
  int maxIndex;

  //Model parameters
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Rl_ = STK::Const::Vector<STK::Real>(cols);
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  // Initializations
  m_Uil_.resize(nbSample_,cols) = 0;
  selectRandomColsFromData(m_Uil_,cols);
  v_Ui_ = STK::sumByRow(m_Uil_);
  randomPoissonParameterRows(m_Akl,cols);
  std::pair<int,int> Label_pair;
  for(int i=0;i<knownLabelsRows_.size();i++)
  {
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }
  //Determine row partition using CEM algorithm with equal proportions
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    m_sumik = m_Uil_*(m_Akl.transpose());
    for ( int i =0;i< UnknownLabelsRows_.size();i++)
    {
      m_sumik.row(UnknownLabelsRows_[i]).maxElt(maxIndex);
      m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
      m_Tik_(UnknownLabelsRows_[i],maxIndex)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In ContingencyLBModel::initCEMRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
#endif
      return false;
    }
    // M-step
    m_Aklold = m_Akl;
    m_Akl = ( (m_Tik_.transpose()*m_Uil_)
             /( (m_Tik_.transpose()*v_Ui_)*STK::Const::Point<STK::Real>(cols))+RealMin).log();

    if((((m_Akl-m_Aklold).abs()/m_Akl).sum())<Mparam_.initepsilon_)
    { break;}
  }
  return true;
}

bool ContingencyLBModel::initCEMCols()
{
  //Temporary variables
  MatrixReal m_Alk(Mparam_.nbcolclust_,Mparam_.nbrowclust_) , m_Alkold(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  int maxIndex;

  //Initializations
  m_Vjk_ = (m_Dataij_.transpose())*m_Tik_;
  v_Vj_  = STK::sumByRow(m_Vjk_);
  randomPoissonParameterCols(m_Alk);
  std::pair<int,int> Label_pair;
  for ( int j=0;j<knownLabelsCols_.size();j++)
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    // CE-step
    m_sumjl = m_Vjk_*(m_Alk.transpose());
    for ( int j =0;j< UnknownLabelsCols_.size();j++)
    {
      m_sumjl.row(UnknownLabelsCols_[j]).maxElt(maxIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
      m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepCols()) )
    {
      Error_msg_  = "In ContingencyLBModel::initCEMCols(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
#endif
      return false;
    }
    // M-step
    m_Alkold = m_Alk;
    m_Alk = ( (m_Rjl_.transpose()*m_Vjk_)
             /( (m_Rjl_.transpose()*v_Vj_)*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_))+RealMin).log();
    if((((m_Alk-m_Alkold).abs()/m_Alk).sum())<Mparam_.initepsilon_)
    { break;}
  }
  // try some optimization
  if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
  { m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
  else
  { m_Gammakl_ = (m_Tik_.transpose()* (m_Dataij_)*m_Rjl_);}
  m_Gammakl_ /= STK::sumByRow(m_Gammakl_)*STK::sum(m_Gammakl_);
  return true;
}

void ContingencyLBModel::logSumRows(MatrixReal & m_ik)
{
  m_ik = STK::Const::Vector<STK::Real>(nbSample_)*v_logPiek_.transpose()
       + m_Uil_* m_Gammakl_.log().transpose();
}

void ContingencyLBModel::logSumCols(MatrixReal & m_jl)
{
  m_jl = STK::Const::Vector<STK::Real>(nbVar_)*v_logRhol_.transpose()
       + m_Vjk_* m_Gammakl_.log();
}

void ContingencyLBModel::randomPoissonParameterRows(MatrixReal& m_kl,int cols)
{
  int index;
  VectorInteger _v_temp = randSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    index=_v_temp[k];
    //index=k;
    for ( int l = 0; l < cols; ++l)
    {
      m_kl(k,l)=std::log(m_Uil_(index,l)/v_Ui_[index]+RealMin);
    }
  }
}

void ContingencyLBModel::randomPoissonParameterCols(MatrixReal& m_lk)
{
  int index;
  VectorInteger _v_temp = randSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    index=_v_temp[l];
    //index=l;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      m_lk(l,k) = std::log(m_Vjk_(index,k)/v_Vj_[index]+RealMin);
    }
  }
}

void ContingencyLBModel::finalizeOutput()
{
  commonFinalizeOutput();
}

void ContingencyLBModel::consoleOut()
{
#ifndef COVERBOSE
  std::cout<<"Output Model parameter:"<<"\ngammakl:\n"<<m_Gammakl_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}
MatrixReal const& ContingencyLBModel::arrangedDataClusters()
{
  arrangedDataCluster<MatrixReal>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

void ContingencyLBModel::modifyThetaStart()
{
  m_Gammaklstart_ = m_Gammakl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
}

void ContingencyLBModel::copyThetaStart()
{
  m_Gammakl_ = m_Gammaklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = STK::sum(m_Rjl_);
  m_Gammakl1_ = m_Gammakl_;
}

void ContingencyLBModel::copyThetaMax()
{
  m_Gammakl_ = m_Gammaklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  likelihood_ = Lmax_;
}

void ContingencyLBModel::modifyThetaMax()
{
  m_Gammaklmax_ = m_Gammakl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = likelihood_;
}

void ContingencyLBModel::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/nbVar_).log();
    v_logPiek_=(v_Tk_/nbSample_).log();
  }
  // try some optimization
  if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
  { m_Ykl_     = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
  else
  { m_Ykl_     = m_Tik_.transpose()*(m_Dataij_*m_Rjl_);}
  m_Gammakl_ = m_Ykl_/( STK::sumByRow(m_Ykl_)* STK::sum(m_Dataij_*m_Rjl_));
}
