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

#include <math.h>
#include "CategoricalLBModel.h"

CategoricalLBModel::CategoricalLBModel( MatrixInteger const& m_Dataij
                                      , ModelParameters const& Mparam
                                      , int a, int b)
                                      : ICoClustModel(Mparam) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  int maxr =  m_Dataij_.maxElt();
  int minr =  m_Dataij_.minElt();
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

  for (int i = 0; i < Mparam.nbrowdata_; ++i)
  { m3_Yijh_[i].resize(Mparam_.nbcoldata_,r_);}

  for (int j = 0; j < Mparam.nbcoldata_; ++j)
  { m3_Yjih_[j].resize(Mparam_.nbrowdata_,r_);}

  for (int h = 0; h < r_; ++h)
  {
    m3_Yhij_[h] = (m_Dataij_ == minr+h);
    for (int i = 0; i < Mparam.nbrowdata_; ++i)
    {
      for (int j = 0; j < Mparam.nbcoldata_; ++j)
      {
        m3_Yijh_[i](j,h) = m3_Yhij_[h](i,j);
        m3_Yjih_[j](i,h) = m3_Yhij_[h](i,j);
      }
    }
  }
}

CategoricalLBModel::CategoricalLBModel( MatrixInteger const& m_Dataij
                                      , VectorInteger const & rowlabels
                                      , VectorInteger const & collabels
                                      , ModelParameters const& Mparam
                                      , int a, int b)
                                      : ICoClustModel(Mparam,rowlabels,collabels) , m_Dataij_(m_Dataij)
{
  a_ = a;
  b_ = b;
  int maxr =  m_Dataij_.maxElt();
  int minr =  m_Dataij_.minElt();
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
  for (int i = 0; i < Mparam.nbrowdata_; ++i)
  { m3_Yijh_[i].resize(Mparam_.nbcoldata_,r_);}
  for (int j = 0; j < Mparam.nbcoldata_; ++j)
  { m3_Yjih_[j].resize(Mparam_.nbrowdata_,r_);}
  for (int h = 0; h < r_; ++h)
  {
    m3_Yhij_[h] = (m_Dataij_ == minr+h);
    for (int i = 0; i < Mparam.nbrowdata_; ++i)
    {
      for (int j = 0; j < Mparam.nbcoldata_; ++j)
      {
        m3_Yijh_[i](j,h) = m3_Yhij_[h](i,j);
        m3_Yjih_[j](i,h) = m3_Yhij_[h](i,j);
      }
    }
  }
}

CategoricalLBModel::~CategoricalLBModel() {}

void CategoricalLBModel::logSumRows(MatrixReal & m_sum)
{
  std::vector<MatrixReal> v_mlogtemphl(Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    v_mlogtemphl[k].resize(r_,Mparam_.nbcolclust_);
    for (int h = 0; h < r_; ++h)
    {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l)
      { v_mlogtemphl[k](h,l) = m3_logAlhphahkl_[h](k,l);}
    }
  }
  for (int i = 0; i < Mparam_.nbrowdata_; ++i)
  {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k)
    { m_sum(i,k) = v_logPiek_[k] + (m3_Yijh_[i].cast<STK::Real>()*v_mlogtemphl[k]).sum();}
  }
}

void CategoricalLBModel::logSumCols(MatrixReal & m_sum)
{
  std::vector<MatrixReal> v_mlogtemphk(Mparam_.nbcolclust_);
  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    v_mlogtemphk[l].resize(r_,Mparam_.nbrowclust_);
    for (int h = 0; h < r_; ++h)
    {
      for (int k = 0; k < Mparam_.nbrowclust_; ++k)
      { v_mlogtemphk[l](h,k) = m3_logAlhphahkl_[h](k,l);}
    }
  }
  for (int j = 0; j < Mparam_.nbcoldata_; ++j)
  {
    for (int l = 0;l < Mparam_.nbcolclust_; ++l)
    { m_sum(j,l) = v_logRhol_[l] + (m3_Yjih_[j].cast<STK::Real>()*v_mlogtemphk[l]).sum();}
  }
}

void CategoricalLBModel::mStepRows()
{
  if(!Mparam_.fixedproportions_)
  { v_logPiek_=((v_Tk_+a_-1)/(nbSample_+Mparam_.nbrowclust_*(a_-1))).log();}

  Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
}

void CategoricalLBModel::mStepCols()
{
  if(!Mparam_.fixedproportions_)
  { v_logRhol_=((v_Rl_+a_-1)/(nbVar_+Mparam_.nbcolclust_*(a_-1))).log();}

  Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
}

void CategoricalLBModel::modifyThetaStart()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::modifyThetaStart().\n";
#endif
  m3_Alphahklstart_ = m3_Alphahkl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
  m_Tikstart_ = m_Tik_;
}

void CategoricalLBModel::copyThetaStart()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::copyThetaStart().\n";
#endif
  m3_Alphahkl_ = m3_Alphahklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;
  m_Rjl_ = m_Rjlstart_;
  m_Tik_ = m_Tikstart_;
  //initialization
  v_Rl_ = STK::sum(m_Rjl_);
}

void CategoricalLBModel::copyThetaMax()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::copyThetaMax().\n";
#endif
  m3_Alphahkl_ = m3_Alphahklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;
  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  likelihood_ = Lmax_;
}

void CategoricalLBModel::modifyThetaMax()
{
#ifdef COVERBOSE
  std::cout<<"Entering CategoricalLBModel::modifyThetaMax().\n";
#endif
  m3_Alphahklmax_ = m3_Alphahkl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;
  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = likelihood_;
}


bool CategoricalLBModel::emRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepRows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepRows();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    {
      netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();
    }
    netchange/=r_;
    //Termination check
    if (netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::cemRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!ceStepRows()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepRows();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if (netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::semRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepRows()) return false;
  mStepRows();
  return true;
}

bool CategoricalLBModel::emCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepCols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepCols();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_ = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::cemCols(){
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!ceStepCols()) return false;
    //M-step
    m3_Alphahklold_ = m3_Alphahkl_;
    mStepCols();
    STK::Real netchange = 0.0;
    for (int h = 0; h < r_; ++h)
    { netchange+= ((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum();}
    netchange/=r_;
    //Termination check
    if(netchange<Mparam_.epsilon_int_) break;
  }
  // Update Alpha for outer loop
  m3_Alphahkl1old_ = m3_Alphahkl1_;
  m3_Alphahkl1_    = m3_Alphahkl_;
  return true;
}

bool CategoricalLBModel::semCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepCols()) return false;
  mStepCols();
  return true;
}

STK::Real CategoricalLBModel::estimateLikelihood()
{
  Array2DReal m_TbyRkl = v_Tk_*v_Rl_.transpose();
  STK::Real tempsum = -(m_TbyRkl.prod(m_TbyRkl+RealMin).log()).sum();
  for (int h = 0; h < r_; ++h)
  {
    Array2DReal m_Ukl = m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>()*m_Rjl_;
    tempsum+= ( m_Ukl.prod( (m_Ukl+RealMin).log() )).sum()+(b_-1)*((m3_Alphahkl_[h]+RealMin).log()).sum();
  }
  likelihood_ = tempsum
              + v_Tk_.dot(v_logPiek_) - Mparam_.nbrowdata_*log(STK::Real(Mparam_.nbrowdata_))
              + v_Rl_.dot(v_logRhol_) - Mparam_.nbcoldata_*log(STK::Real(Mparam_.nbcoldata_))
              - (m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum()
              + (a_-1)*(v_logPiek_.sum()+v_logRhol_.sum());
  return likelihood_;
}

void CategoricalLBModel::parameterStopCriteria()
{
  STK::Real netchange = 0.0;
  for (int h = 0; h < r_; ++h)
  { netchange+= ((m3_Alphahkl1_[h]-m3_Alphahkl1old_[h]).abs()/(m3_Alphahkl1_[h]+RealMin)).sum();}
  netchange/=r_;
  stopAlgo_ = (netchange<Mparam_.epsilon_);
}

bool CategoricalLBModel::randomInit()
{
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
  { Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1.0;}

  for ( int j =0;j< (int)UnknownLabelsCols_.size();j++)
  { m_Rjl_(j,v_Wj_[UnknownLabelsCols_[j]]-1) = 1.0;}

  v_Tk_ = STK::sum(m_Tik_).transpose();
  v_Rl_ = STK::sum(m_Rjl_).transpose();
  //Check for empty cluster
  if( (empty_cluster_ = (v_Rl_<3).any()) )
  {
    Error_msg_  = "CategoricalLBModel::randomInit(). Empty column while running initialization.\n";
#ifdef COVERBOSE
  std::cout<<Error_msg_;
#endif
    return false;
  }
  if((v_Tk_<3).any())
  {
    Error_msg_  = "CategoricalLBModel::randomInit(). Empty row while running initialization.\n";
#ifdef COVERBOSE
  std::cout<<Error_msg_;
#endif
    return false;
  }
  //Initializing model parameters
  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);
  m3_Alphahkl1_.resize(r_);
  m3_Alphahklold_.resize(r_);
  m3_Alphahkl1old_.resize(r_);
  Array2DReal m_vtkrl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahkl1old_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>()*m_Rjl_)+b_-1)/(m_vtkrl+RealMin);
    m3_Alphahkl1_[h] = m3_Alphahkl_[h];
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;
  v_Zi_.resize(nbSample_) = 0;
  v_Wj_.resize(nbVar_) = 0;
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_));
  return true;
}

STK::Real CategoricalLBModel::iclCriteriaValue()
{
  STK::Real criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
          -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
          +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(r_*b_)-r_*lgamma(b_))
          -lgamma(Mparam_.nbrowdata_+Mparam_.nbrowclust_*a_)
          -lgamma(Mparam_.nbcoldata_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  { criteria+= lgamma(a_+ (v_Zi_== k).count());}
  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  { criteria+= lgamma(a_+ (v_Wj_==l).count());}

  STK::ArrayXXi temp(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int h = 0; h < r_; ++h)
  {
    temp = ((m_Zik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Wjl_)+b_;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l)
      { criteria+=lgamma(temp(k,l));}
    }
  }
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    { criteria-= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+r_*b_);}
  }
  return criteria;
}

void CategoricalLBModel::finalizeOutput()
{ commonFinalizeOutput();}

void CategoricalLBModel::consoleOut()
{
#ifdef COVERBOSE
  std::cout<<"Output Model parameters\n";
  std::cout<<"\npie_k:"<<v_Piek_<<"\nrho_l:"<<v_Rhol_<<"\n";
  for (int h = 0; h < r_; ++h)
  {
    std::cout<<"Alpha_kl for category "<<h<<"\n";
    std::cout<<m3_Alphahkl_[h]<<"\n";
  }
  std::cout<<"likelihood value:"<<likelihood()<<"\n";
  std::cout<<"ICL value:"<<iclCriteriaValue()<<"\n";
#endif
}

const MatrixInteger& CategoricalLBModel::arrangedDataClusters()
{
  arrangedDataCluster(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

