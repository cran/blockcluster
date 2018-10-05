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

CategoricalLBModel::CategoricalLBModel( MatrixInt const& m_Dataij
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

CategoricalLBModel::CategoricalLBModel( MatrixInt const& m_Dataij
                                      , VectorInt const & rowlabels
                                      , VectorInt const & collabels
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
  m_sum = STK::Const::VectorX(nbSample_)*v_logPiek_.transpose();
  for (int h = 0; h < r_; ++h) {
    m_sum +=  (m3_Yhij_[h].cast<STK::Real>()*m_Rjl_*m3_logAlhphahkl_[h].transpose());
  }
}

bool CategoricalLBModel::cemInitStep()
{
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
  if (initCEMRows())
  {
#ifdef COVERBOSE
    std::cout << "CategoricalLBModel::initCEMRows done with success."<<std::endl;
    consoleOut();
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    if (initCEMCols())
    {
      m3_Alphahkl1_ = m3_Alphahkl_;
      m3_Alphahkl1old_.resize(r_);
      m3_Alphahklold_.resize(r_);
      m3_logAlhphahkl_.resize(r_);
      for (int h = 0; h < r_; ++h)
      {
        m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        m3_Alphahkl1old_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        m3_Alphahkl1_[h] = m3_Alphahkl_[h];
        m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
      }
      m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
      m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
#ifdef COVERBOSE
      std::cout<<"ContingencyLBModel::cemInitStep. Initialization done with success."<<std::endl;
      consoleOut();
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return true;
    }
  }
  return false;
}

bool CategoricalLBModel::emInitStep()
{
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
  if (initEMRows())
  {
#ifdef COVERBOSE
    std::cout << "CategoricalLBModel::initEMRows done with success."<<std::endl;
    consoleOut();
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    if (initEMCols())
    {
      m3_Alphahkl1_ = m3_Alphahkl_;
      m3_Alphahkl1old_.resize(r_);
      m3_Alphahklold_.resize(r_);
      m3_logAlhphahkl_.resize(r_);
      for (int h = 0; h < r_; ++h)
      {
        m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        m3_Alphahkl1old_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
        m3_Alphahkl1_[h] = m3_Alphahkl_[h];
        m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
      }
      m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
      m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
#ifdef COVERBOSE
      std::cout<<"ContingencyLBModel::emInitStep. Initialization done with success."<<std::endl;
      consoleOut();
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return true;
    }
  }
  return false;
}

void CategoricalLBModel::logSumCols(MatrixReal & m_sum)
{
  m_sum = STK::Const::VectorX(nbVar_)*v_logRhol_.transpose();
  for (int h = 0; h < r_; ++h) {
    m_sum +=  (m3_Yhij_[h].transpose().cast<STK::Real>()*m_Tik_*m3_logAlhphahkl_[h]);
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

void CategoricalLBModel::mGibbsStepRows()
{
  v_logPiek_=(v_Tk_+a_);

  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbrowclust_);
  STK::Real sumRng = 0.0;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_randgamma[k] = STK::Law::Gamma::rand(v_logPiek_[k],1);
    sumRng += v_randgamma[k];
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_logPiek_[k] = v_randgamma[k]/sumRng;
  }
  v_logPiek_ = (v_logPiek_+RealMin).log();

  std::vector<MatrixReal> m_randgamma;
  std::vector<VectorReal> v_sumRng(Mparam_.nbrowclust_);
  m_randgamma.resize(r_);
  v_sumRng.resize(r_);
  for (int h = 0; h < r_; ++h) {
    m_randgamma[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    v_sumRng[h].resize(Mparam_.nbrowclust_);
  }
  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m_randgamma[h](k,l) = STK::Law::Gamma::rand(m3_Alphahkl_[h](k,l),1);
        v_sumRng[h][k] += m_randgamma[h](k,l);
      }
    }
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m3_Alphahkl_[h](k,l) = m_randgamma[h](k,l)/v_sumRng[h][k];
      }
    }
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

}

void CategoricalLBModel::mGibbsStepCols()
{
  v_logRhol_=(v_Rl_+a_);

  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }

  //generate random numbers
  VectorReal v_randgamma(Mparam_.nbrowclust_);
  STK::Real sumRng = 0.0;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_randgamma[k] = STK::Law::Gamma::rand(v_logRhol_[k],1);
    sumRng += v_randgamma[k];
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    v_logRhol_[k] = v_randgamma[k]/sumRng;
  }
  v_logRhol_ = (v_logRhol_+RealMin).log();

  std::vector<MatrixReal> m_randgamma;
  std::vector<VectorReal> v_sumRng(Mparam_.nbcolclust_);
  m_randgamma.resize(r_);
  v_sumRng.resize(r_);
  for (int h = 0; h < r_; ++h) {
    m_randgamma[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    v_sumRng[h].resize(Mparam_.nbrowclust_);
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m_randgamma[h](k,l) = STK::Law::Gamma::rand(m3_Alphahkl_[h](k,l),1);
        v_sumRng[h][k] += m_randgamma[h](k,l);
      }
    }
  }

  for (int h = 0; h < r_; ++h) {
    for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
        m3_Alphahkl_[h](k,l) = m_randgamma[h](k,l)/v_sumRng[h][k];
      }
    }
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

bool CategoricalLBModel::GibbsRows()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepRows()) return false;
  mGibbsStepRows();
  return true;
}

bool CategoricalLBModel::GibbsCols()
{
  //Initializations
  for (int h = 0; h < r_; ++h)
  { m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();}
  if(!seStepCols()) return false;
  mGibbsStepCols();
  return true;
}

STK::Real CategoricalLBModel::estimateLikelihood()
{
  Array2DReal m_Ukl = v_Tk_*v_Rl_.transpose();
  STK::Real tempsum = -(m_Ukl.prod(m_Ukl+RealMin).log()).sum();
  // Compute \sum_h \sum_{i,j,k,l} t_{ik} Y_{ij}^h r_{jl} \log(\alpha_{k,l}^h
  for (int h = 0; h < r_; ++h)
  {
    m_Ukl    = m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>()*m_Rjl_;
    // BUGFIXE (Serge Iovleff) Seems log is wrong
    //    tempsum += ( m_Ukl.prod( (m_Ukl+RealMin).log() )).sum()
    //    		 + (b_-1)*((m3_Alphahkl_[h]+RealMin).log()).sum();
    tempsum += ( m_Ukl.prod( (m3_Alphahkl_[h]+RealMin).log() )).sum();
//    		 + (b_-1)*((m3_Alphahkl_[h]+RealMin).log()).sum();
  }
  likelihood_ = tempsum
              + v_Tk_.dot(v_logPiek_) // - nbSample_*log(STK::Real(nbSample_))
              + v_Rl_.dot(v_logRhol_) // - nbVar_   *log(STK::Real(nbVar_))
              - ( m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              - ( m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum()
			  //              + (a_-1)*(v_logPiek_.sum()+v_logRhol_.sum())
              ;
  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int CategoricalLBModel::nbFreeParameters() const
{ return Mparam_.nbcolclust_ * Mparam_.nbrowclust_ * (r_-1);}

void CategoricalLBModel::parameterStopCriteria()
{
  STK::Real netchange = 0.0;
  for (int h = 0; h < r_; ++h)
  { netchange+= ((m3_Alphahkl1_[h]-m3_Alphahkl1old_[h]).abs()/(m3_Alphahkl1_[h]+RealMin)).sum();}
  netchange/=r_;
  stopAlgo_ = (netchange<Mparam_.epsilon_);
}

//Compute Bernoulli log-sum for all rows
void CategoricalLBModel::initBernoulliLogSumRows(MatrixReal & _m_sum)
{
//std::vector<MatrixReal> v_mlogtemphl(Mparam_.nbrowclust_);
//for (int k = 0; k < Mparam_.nbrowclust_; ++k)
//{
//  v_mlogtemphl[k].resize(r_,Mparam_.nbcolclust_);
//  for (int h = 0; h < r_; ++h)
//  {
//    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
//    { v_mlogtemphl[k](h,l) = m3_logAlhphahkl_[h](k,l);}
//  }
//}
//for (int i = 0; i < Mparam_.nbrowdata_; ++i)
//{
//  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
//  { _m_sum(i,k) = v_logPiek_[k] + (m3_Yijh_[i].cast<STK::Real>()*v_mlogtemphl[k]).sum();}
//}
  _m_sum = STK::Const::VectorX(nbSample_)*v_logPiek_.transpose();
  for (int h = 0; h < r_; ++h) {
    _m_sum +=  (m3_Yhij_[h].cast<STK::Real>()*m_Rjl_*m3_logAlhphahkl_[h].transpose());
  }
}

//Compute Bernoulli log-sum for all columns
void CategoricalLBModel::initBernoulliLogSumCols(MatrixReal & _m_sum)
{
//std::vector<MatrixReal> v_mlogtemphk(Mparam_.nbcolclust_);
//for (int l = 0; l < Mparam_.nbcolclust_; ++l)
//{
//  v_mlogtemphk[l].resize(r_,Mparam_.nbrowclust_);
//  for (int h = 0; h < r_; ++h)
//  {
//    for (int k = 0; k < Mparam_.nbrowclust_; ++k)
//    { v_mlogtemphk[l](h,k) = m3_logAlhphahkl_[h](k,l);}
//  }
//}
//for (int j = 0; j < Mparam_.nbcoldata_; ++j)
//{
//  for (int l = 0;l < Mparam_.nbcolclust_; ++l)
//  { _m_sum(j,l) = v_logRhol_[l] + (m3_Yjih_[j].cast<STK::Real>()*v_mlogtemphk[l]).sum();}
//}
  _m_sum = STK::Const::VectorX(nbVar_)*v_logRhol_.transpose();
  for (int h = 0; h < r_; ++h) {
    _m_sum += (m3_Yhij_[h].transpose().cast<STK::Real>()*m_Tik_*m3_logAlhphahkl_[h]);
  }
}

bool CategoricalLBModel::initCEMRows()
{
  // Initialization of various parameters
//int cols=std::min(100, nbVar_);
  int cols = nbVar_;
  m_Uil_.resize(nbSample_,cols) = 0;
  selectRandomColsFromData(m_Uil_,cols);
  v_Ui_ = STK::sumByRow(m_Uil_);

  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Rl_  = STK::Const::VectorX(cols);

  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;

  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);

  for (int h = 0; h < r_; ++h)
  {
//  m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,cols) = 0;
//  m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,cols);
//  randomParameterRows(m3_Alphahkl_[h],cols);
    m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
    m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    randomParameterRows(m3_Alphahkl_[h],Mparam_.nbcolclust_);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  std::pair<int,int> Label_pair;
  for(int i=0;i<knownLabelsRows_.size();i++)
  {
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    initBernoulliLogSumRows(m_sumik);
    for (int i=0; i<UnknownLabelsRows_.size();i++)
    {
   	  int maxIndex;
      m_sumik.row(UnknownLabelsRows_[i]).maxElt(maxIndex);
      m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
      m_Tik_(UnknownLabelsRows_[i],maxIndex)=1;
    }
    // compute v_Tk_ and check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In CategoricalLBModel::initCEMRows(). Class size too small while estimating model.\n";
  #ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
  #endif
      return false;
    }
    // M-step
    v_Rl_ = STK::sum(m_Rjl_);
    m3_Alphahklold_ = m3_Alphahkl_;
    Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
    for (int h = 0; h < r_; ++h)
    {
      m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
//    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_);
      if((((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum())<Mparam_.initepsilon_)
      { break;}
    }
  }
  return true;
}

bool CategoricalLBModel::initCEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;
  v_Vj_  = STK::sumByRow(m_Vjk_);

  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    randomParameterCols(m3_Alphahkl_[h]);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
  std::pair<int,int> Label_pair;
  for (int j=0;j<knownLabelsCols_.size();j++)
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    // CE-step
    initBernoulliLogSumCols(m_sumjl);
    for (int j=0;j< UnknownLabelsCols_.size();j++)
    {
      int maxIndex;
      m_sumjl.row(UnknownLabelsCols_[j]).maxElt(maxIndex);
      m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
      m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
    }
    // compute v_Rl_ and check empty class
    if( (empty_cluster_ = finalizeStepCols()) )
    {
      Error_msg_  = "In CategoricalLBModel::initCEMCols(). Class size too small while estimating model.\n";
  #ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
  #endif
      return false;
    }
    // M-step
    v_Tk_ = STK::sum(m_Tik_);
    m3_Alphahklold_ = m3_Alphahkl_;
    Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
    for (int h = 0; h < r_; ++h)
    {
      m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
      if((((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum())<Mparam_.initepsilon_)
      { break;}
    }
  }
  return true;
}

bool CategoricalLBModel::initEMRows()
{
  // Initialization of various parameters
//int cols=std::min(100, nbVar_);
  int cols = /*Mparam_.nbcolclust_;*/ nbVar_;
  m_Uil_.resize(nbSample_,cols) = 0;
  selectRandomColsFromData(m_Uil_,cols);
  v_Ui_ = STK::sumByRow(m_Uil_);

  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Rl_  = STK::Const::VectorX(cols);

  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;

  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);

  for (int h = 0; h < r_; ++h)
  {
//  m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,cols) = 0;
//  m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,cols);
//  randomParameterRows(m3_Alphahkl_[h],cols);
    m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
    m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    randomParameterRows(m3_Alphahkl_[h],Mparam_.nbcolclust_);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
  //Determine row partition using EM algorithm with equal proportions
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  std::pair<int,int> Label_pair;
  for(int i=0;i<knownLabelsRows_.size();i++)
  {
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    initBernoulliLogSumRows(m_sumik);
    m_Tik_  = (m_sumik-STK::maxByRow(m_sumik)*STK::Const::PointX(Mparam_.nbrowclust_)).exp();
    m_Tik_ /= STK::sumByRow(m_Tik_)*STK::Const::PointX(Mparam_.nbrowclust_);
    // reinitialize known labels
    for ( int i=0;i< (int)knownLabelsRows_.size();i++)
    {
      m_Tik_.row(knownLabelsRows_[i].first).setZeros();
      m_Tik_(knownLabelsRows_[i].first, knownLabelsRows_[i].second)=1;
    }


/*///////TEST S STEP
  //take cumulative sum of probabilities
  MatrixReal m_Tiktemp = m_Tik_;
  for ( int k = 1; k < m_Tiktemp.sizeCols(); ++k)
  { m_Tiktemp.col(k) += m_Tiktemp.col(k-1);}

  //generate random numbers
  std::vector<STK::Real> randnumbers(nbSample_);
  for ( int i = 0; i < nbSample_; ++i)
  {
#ifdef _OPENMP
#pragma omp critical
#endif
    randnumbers[i] = STK::Law::Uniform::rand(0,1);
  }

  m_Zik_.setZeros();
  //chose randomly the row class using generated random numbers
  for ( int i = 0; i < nbSample_; ++i)
  {
    for ( int k = 0; k < m_Tiktemp.sizeCols(); ++k)
    {
      if(randnumbers[i]< m_Tiktemp(i,k))
      {
        m_Zik_(i,k) = 1;
        break;
      }
    }
  }
  m_Tik_ = m_Zik_.cast<STK::Real>();
//////FIN TEST*/


    // check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In CategoricalLBModel::InitEMRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return false;
    }

    // M-step
    v_Rl_ = STK::sum(m_Rjl_);
    Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
    m3_Alphahklold_ = m3_Alphahkl_;
    for (int h = 0; h < r_; ++h)
    {
      m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
//    m3_Alphahkl_[h] = ((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_);
      if((((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum())<Mparam_.initepsilon_)
      { break;}
    }
  }
  return true;
}

bool CategoricalLBModel::initEMCols()
{
  //Determine row partition using EM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  m3_Alphahkl_.resize(r_);
  m3_logAlhphahkl_.resize(r_);
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;
  v_Vj_  = STK::sumByRow(m_Vjk_);
  for (int h = 0; h < r_; ++h)
  {
    m3_Alphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_logAlhphahkl_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    m3_Alphahklold_[h].resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
    randomParameterCols(m3_Alphahkl_[h]);
    m3_logAlhphahkl_[h] = (m3_Alphahkl_[h]+RealMin).log();
  }
  consoleOut();
  std::pair<int,int> Label_pair;
  for (int j=0;j<knownLabelsCols_.size();j++)
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    //E-step
    initBernoulliLogSumCols(m_sumjl);
    m_Rjl_ = (m_sumjl-STK::maxByRow(m_sumjl)*STK::Const::PointX(Mparam_.nbcolclust_)).exp();
    m_Rjl_ /= (STK::sumByRow(m_Rjl_)*STK::Const::PointX(Mparam_.nbcolclust_));
    //
    for ( int j=0;j< (int)knownLabelsCols_.size();j++)
    {
      m_Rjl_.row(knownLabelsCols_[j].first).setZeros();
      m_Rjl_(knownLabelsCols_[j].first,knownLabelsCols_[j].second)=1;
    }


/*///////TEST S STEP
  //take cumulative sum of probabilities
  MatrixReal m_Rjltemp = m_Rjl_;
  for ( int l = 1; l < m_Rjltemp.sizeCols(); ++l)
  { m_Rjltemp.col(l) += m_Rjltemp.col(l-1);}

  //generate random numbers
  std::vector<STK::Real> randnumbers(nbSample_);
  for ( int j = 0; j < nbVar_; ++j)
  {
    //std::srand(j);
//  randnumbers[j] = STK::Real(std::rand())/STK::Real(RAND_MAX);
#ifdef _OPENMP
#pragma omp critical
#endif
    randnumbers[j] = STK::Law::Uniform::rand(0,1);
  }
  m_Wjl_.setZeros();
  //chose randomly the row class using generated random numbers
  for ( int j = 0; j < nbVar_; ++j)
  {
    for ( int l = 0; l < m_Rjltemp.sizeCols(); ++l)
    {
      if(randnumbers[j]<m_Rjltemp(j,l))
      {
        m_Wjl_(j,l) = 1;
        break;
      }
    }
  }

  m_Rjl_ = m_Wjl_.cast<STK::Real>();
///////FIN TEST*/


    // check empty class
    if( (empty_cluster_ = finalizeStepCols()) )
    {
      Error_msg_  = "In CategoricalLBModel::InitEMCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return false;
    }

    // M-step
    v_Tk_ = STK::sum(m_Tik_);
    Array2DReal m_TbyRkl = (v_Tk_*v_Rl_.transpose())+r_*(b_-1);
    m3_Alphahklold_ = m3_Alphahkl_;
    for (int h = 0; h < r_; ++h)
    {
      m3_Alphahkl_[h] = (((m_Tik_.transpose()*m3_Yhij_[h].cast<STK::Real>())*m_Rjl_)+b_-1)/(m_TbyRkl+RealMin);
      if((((m3_Alphahkl_[h]-m3_Alphahklold_[h]).abs()/(m3_Alphahkl_[h]+RealMin)).sum())<Mparam_.initepsilon_)
      { break;}
    }
  }
  return true;
}

void CategoricalLBModel::selectRandomColsFromData(MatrixReal& _m_il,int cols)
{
  if(cols==Mparam_.nbcoldata_)
  { _m_il=m_Dataij_.cast<STK::Real>();}
  else
  {
    //random shuffle Algorithm
    VectorInt _v_temp = randSample(nbVar_,cols);
    for ( int l = 0; l < cols; ++l)
    { _m_il.col(l)=m_Dataij_.cast<STK::Real>().col(_v_temp[l]);}
  }
}

void CategoricalLBModel::randomParameterRows(MatrixReal& _m_kl,int cols)
{
  int index;
  STK::Real epsilon = 0.1;
  VectorInt _v_temp = randSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    index=_v_temp[k];
    //index=k;
    for ( int l = 0; l < cols; ++l)
    { _m_kl(k,l)=epsilon*(1.0-m_Uil_(index,l))+m_Uil_(index,l)*(1.0-epsilon);}
  }
}

void CategoricalLBModel::randomParameterCols(MatrixReal& _m_kl)
{
  int index;
  VectorInt _v_temp = randSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    index=_v_temp[l];
    //index=l;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    { _m_kl(k,l)=m_Vjk_(index,k)/v_Tk_[k];}
  }
}

bool CategoricalLBModel::randomInit()
{
  std::vector<MatrixReal> test(Mparam_.nbrowclust_);
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters..with random"<<std::endl;
#endif
  //Initialize random row and column partition
  VectorReal probarows = (1.0/Mparam_.nbrowclust_)*STK::Const::VectorX(Mparam_.nbrowclust_);
  VectorReal probacols = (1.0/Mparam_.nbcolclust_)*STK::Const::VectorX(Mparam_.nbcolclust_);
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
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
  return true;
}

STK::Real CategoricalLBModel::iclCriteriaValue()
{
  STK::Real criteria = lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
                     - (Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
                     + Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(r_*b_)-r_*lgamma(b_))
                     - lgamma(Mparam_.nbrowdata_+Mparam_.nbrowclust_*a_)
                     - lgamma(Mparam_.nbcoldata_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  { criteria += lgamma(a_+ (v_Zi_== k).count());}
  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  { criteria += lgamma(a_+ (v_Wj_==l).count());}

  STK::ArrayXXi temp(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int h = 0; h < r_; ++h)
  {
    temp = ((m_Zik_.transpose()*m3_Yhij_[h].cast<int>())*m_Wjl_)+b_;
    for (int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for (int l = 0; l < Mparam_.nbcolclust_; ++l)
      { criteria += lgamma(temp(k,l));}
    }
  }
  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    { criteria -= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+r_*b_);}
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
  std::cout << "v_Tk_= " << v_Tk_.transpose();
  std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
}

const MatrixInt& CategoricalLBModel::arrangedDataClusters()
{
  arrangedDataCluster(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

