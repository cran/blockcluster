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


/** @file BinaryLBModelequalepsilon.cpp
 *  @brief Implements concrete BinaryLBModelequalepsilon model class derived from ICoClustModel.
 **/

#include "BinaryLBModelequalepsilon.h"

#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModelequalepsilon::BinaryLBModelequalepsilon( MatrixBinary const&  m_Dataij
                                                    , ModelParameters const& Mparam
                                                    ,int a,int b)
                                                    : ICoClustModel(Mparam)
                                                    , a_(a), b_(b)
                                                    , m_Dataij_(m_Dataij)
                                                    , Epsilon_(0),Epsilonstart_(0),Epsilonmax_(0)
                                                    , W1_(0),W1_old_(0)
{
  m_Xjl_ = STK::sum(m_Dataij.cast<STK::Real>()).transpose()*STK::Const::Point<STK::Real>(Mparam_.nbcolclust_);
  m_Xik_ = STK::sumByRow(m_Dataij.cast<STK::Real>())*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_);
};

BinaryLBModelequalepsilon::BinaryLBModelequalepsilon( MatrixBinary const&  m_Dataij
                                                    , VectorInt const& rowlabels
                                                    , VectorInt const& collabels
                                                    , ModelParameters const& Mparam
                                                    , int a,int b)
                                                    : ICoClustModel(Mparam,rowlabels,collabels)
                                                    , a_(a), b_(b)
                                                    , m_Dataij_(m_Dataij)
                                                    , Epsilon_(0),Epsilonstart_(0),Epsilonmax_(0)
                                                    , W1_(0),W1_old_(0)
{
  m_Xjl_ = STK::sum(m_Dataij.cast<STK::Real>()).transpose()*STK::Const::Point<STK::Real>(Mparam_.nbcolclust_);
  m_Xik_ = STK::sum(m_Dataij.cast<STK::Real>())*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_);
};

bool BinaryLBModelequalepsilon::cemInitStep()
{
  // initialize parameters
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  v_Rl_.resize(Mparam_.nbcolclust_) = 0;
  if (initCEMRows())
  {
#ifdef COVERBOSE
  std::cout << "BinaryLBModelequalepsilon::initCEMRows done with success."<<std::endl;
  consoleOut();
#endif
    if (initCEMCols())
    {
      m_Ykl_old1_ = m_Ykl_;
      m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
      m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
      v_Zi_.resize(nbSample_) = 0;
      v_Wj_.resize(nbVar_) = 0;
      m_Zik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
      m_Wjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
#ifdef COVERBOSE
      std::cout<<"BinaryLBModelequalepsilon::cemInitStep. Initialization done with success."<<std::endl;
      consoleOut();
#endif

      return true;
    }
  }
  return false;
}

bool BinaryLBModelequalepsilon::emInitStep()
{
  // initialize parameters
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;
  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  v_Rl_.resize(Mparam_.nbcolclust_) = 0;
  if (initEMRows())
  {
#ifdef COVERBOSE
  std::cout << "BinaryLBModelequalepsilon::initEMRows done with success."<<std::endl;
  consoleOut();
#endif
    if (initEMCols())
    {
      m_Ykl_old1_ = m_Ykl_;
      m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
      m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
      v_Zi_.resize(nbSample_) = 0;
      v_Wj_.resize(nbVar_) = 0;
      m_Zik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
      m_Wjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
      v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
      v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
#ifdef COVERBOSE
      std::cout<<"BinaryLBModelequalepsilon::emInitStep. Initialization done with success."<<std::endl;
      consoleOut();
#endif

      return true;
    }
  }
  return false;
}


void BinaryLBModelequalepsilon::finalizeOutput()
{
  commonFinalizeOutput();
}

void BinaryLBModelequalepsilon::consoleOut()
{
#ifdef COVERBOSE
  std::cout<<"Output Model parameter:"<<"\nakl:\n"<<m_Akl_<<"\nepsilon:\n"<<Epsilon_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
  std::cout<<"ICL: "<<iclCriteriaValue()<<"\n";
#endif
}

//Compute Bernoulli log-sum for all rows
void BinaryLBModelequalepsilon::logSumRows(MatrixReal & m_sum)
{
  STK::Real logepsilon = log(Epsilon_/(-Epsilon_+1));
  m_sum = STK::Const::VectorX(nbSample_)*(v_logPiek_+logepsilon*m_Akl_.cast<STK::Real>()*v_Rl_).transpose()
          -logepsilon*(2*m_Uil_*m_Akl_.cast<STK::Real>().transpose() + m_Xik_.cast<STK::Real>());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModelequalepsilon::logSumCols(MatrixReal & m_sum)
{
  STK::Real logepsilon = log(Epsilon_/(-Epsilon_+1));
  m_sum = STK::Const::VectorX(nbVar_)*(v_logRhol_.transpose()+logepsilon*v_Tk_.transpose()*m_Akl_.cast<STK::Real>())
          -logepsilon*(2*m_Vjk_*m_Akl_.cast<STK::Real>() + m_Xjl_.cast<STK::Real>());
}

void BinaryLBModelequalepsilon::initBernoulliLogSumRows(MatrixReal & m_sum)
{
  int cols = m_Uil_.sizeCols();
  m_sum = STK::Const::Array<STK::Real>(nbSample_,cols)*m_Akl_.cast<STK::Real>().transpose()
        - 2*m_Uil_*m_Akl_.cast<STK::Real>().transpose()
        + v_Ui_*STK::Const::Point<STK::Real>(Mparam_.nbrowclust_);
}

void BinaryLBModelequalepsilon::initBernoulliLogSumCols(MatrixReal & m_sum)
{
  m_sum = STK::Const::VectorX(nbVar_)* (v_Tk_.transpose()*m_Akl_.cast<STK::Real>())
        - 2*m_Vjk_*m_Akl_.cast<STK::Real>()
        + m_Xjl_;
}

//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::emRows()
{
  //Initializations
  m_Ykl_ = (m_Akl_.cast<STK::Real>()*(-Epsilon_+1)+(-m_Akl_.cast<STK::Real>()+1)*Epsilon_)*dimprod_;
  m_Uil_ = m_Dataij_.cast<STK::Real>()*m_Rjl_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr)
  {
    //E-step
    if(!eStepRows()) return false;
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepRows();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }
  return true;
}

//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::cemRows()
{
  //Initializations
  m_Ykl_ = (m_Akl_.cast<STK::Real>()*(-Epsilon_+1)+(-m_Akl_.cast<STK::Real>()+1)*Epsilon_)*dimprod_;
  m_Uil_=m_Dataij_.cast<STK::Real>()*m_Rjl_;
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_);
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ceStepRows()) return false;
    //M-step
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepRows();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
       break;
    }
  }
  return true;
}

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::emCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;

  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!eStepCols()) return false;
     //M-step

    if(!Mparam_.fixedproportions_) {
      v_logRhol_=(v_Rl_/nbVar_).log();
    }
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepCols();
    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }
  return true;
}

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::cemCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  for ( int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
    //E-step
    if(!ceStepCols()) return false;
    //M-step
    m_Ykl_old2_ = m_Ykl_;
    mStepCols();

    if(((m_Ykl_-m_Ykl_old2_)/m_Ykl_).abs().sum()<Mparam_.epsilon_int_)
    {
      break;
    }
  }

  return true;
}

bool BinaryLBModelequalepsilon::semRows()
{
  //Initializations
  m_Ykl_ = (m_Akl_.cast<STK::Real>()*(-Epsilon_+1)+(-m_Akl_.cast<STK::Real>()+1)*Epsilon_)*dimprod_;
  m_Uil_=m_Dataij_.cast<STK::Real>()*m_Rjl_;

  if(!seStepRows()) return false;
  //M-step
  mStepRows();
  return true;
}

bool BinaryLBModelequalepsilon::semCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;

  if(!seStepCols()) return false;
  //M-step
  mStepCols();
  return true;
}

bool BinaryLBModelequalepsilon::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool BinaryLBModelequalepsilon::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

STK::Real BinaryLBModelequalepsilon::estimateLikelihood()
{
  likelihood_ = (dimprod_*(Epsilon_*std::log(Epsilon_/(-Epsilon_+1)) + std::log(-Epsilon_+1))
              + v_Tk_.dot(v_logPiek_)
			  + v_Rl_.dot(v_logRhol_)
              - (m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
              - (m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum()
			    );
  return likelihood_;
}
/** @return the number of free parameters of the distribution of a block.*/
int BinaryLBModelequalepsilon::nbFreeParameters() const
{ return 1;}

//Compute change in Alpha and set the terminate variable accordingly
void BinaryLBModelequalepsilon::parameterStopCriteria()
{
  STK::Real relativechange = ((m_Ykl_-m_Ykl_old1_)/m_Ykl_old1_).abs().sum();
  if(relativechange<Mparam_.epsilon_)
    stopAlgo_ = true;
  else
    stopAlgo_ = false;

  // Update Ykl for outer loop
  m_Ykl_old1_ = m_Ykl_;
}

MatrixBinary const& BinaryLBModelequalepsilon::arrangedDataClusters()
{
  arrangedDataCluster<MatrixBinary>(m_ClusterDataij_,m_Dataij_);
  return m_ClusterDataij_;
}

#ifndef RPACKAGE
void BinaryLBModelequalepsilon::displayCluster()
{
  CImg<unsigned char>  cluster(nbVar_,nbSample_,1,1,0);
  CImg<unsigned char>  data(nbVar_,nbSample_,1,1,0);

  MatrixBinary m_ClusterDataij_ = arrangedDataClusters();

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

bool BinaryLBModelequalepsilon::initCEMRows()
{
  // Initialization of various parameters
//int cols=std::min(100, nbVar_);
  int cols =  nbVar_;
  m_Uil_.resize(nbSample_,cols) = 0;
  selectRandomColsFromData(m_Uil_,cols);
  v_Ui_ = STK::sumByRow(m_Uil_);
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  m_Akl_.resize(Mparam_.nbrowclust_,cols) = 0;
  generateRandomMean(m_Akl_);
  W1_ = RealMax;
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  int minIndex;
  std::pair<int,int> Label_pair;
  for(int i = 0;i< (int)knownLabelsRows_.size();i++){
    Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }

  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    initBernoulliLogSumRows(m_sumik);
    for ( int i =0;i< (int)UnknownLabelsRows_.size();i++) {
      m_sumik.row(UnknownLabelsRows_[i]).minElt(minIndex);
      m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
      m_Tik_(UnknownLabelsRows_[i],minIndex)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In BinaryLBModelequalepsilon::initCEMRows(). Class size too small while initializing model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
#endif
      return false;
    }
    // M-step
    W1_old_ = W1_;
    W1_ = (m_Tik_.prod(m_sumik)).sum();
    m_Ukl_ = m_Tik_.transpose()*m_Uil_;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for ( int l = 0; l < cols; ++l)
      {
        m_Akl_(k,l)  = (m_Ukl_(k,l)>=v_Tk_[k]/2) ? 1 : 0;
      }
    }
    //m_Akl_ = m_Ukl_ >=(v_Tk_*MatrixReal(1,cols)/2);
    if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_)
    { break;}
  }
  return true;
}

bool BinaryLBModelequalepsilon::initCEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  MatrixReal m_Tk_Rl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixReal m_Ulk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  int minIndex;
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;
  SelectRandomRows(m_Ulk);
  m_Ukl_ = m_Ulk.transpose();
  m_Akl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      m_Akl_(k,l) = (m_Ukl_(k,l)>=v_Tk_[k]/2) ? 1: 0;
    }
  }

  std::pair<int,int> Label_pair;
  for ( int j=0;j< (int)knownLabelsCols_.size();j++ )
  {
    Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
    for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
    {
      // CE-step
      initBernoulliLogSumCols(m_sumjl);
      for ( int j=0 ;j< (int)UnknownLabelsCols_.size();j++)
      {
        m_sumjl.row(UnknownLabelsCols_[j]).minElt(minIndex);
        m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
        m_Rjl_(UnknownLabelsCols_[j],minIndex)=1;
      }
      // check empty class
      if( (empty_cluster_ = finalizeStepCols()) )
      {
        Error_msg_  = "In BinaryLBModelequalepsilon::initCEMCols(). Class size too small while initializing model.\n";
  #ifdef COVERBOSE
        std::cout << Error_msg_;
  #endif
        return false;
      }
      // M-step
      W1_old_ = W1_;
      W1_ = (m_Rjl_.prod(m_sumjl)).sum();
      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Tk_Rl = v_Tk_*v_Rl_.transpose()/2.0;
      for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
      {
        for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
        {
          m_Akl_(k,l) = (m_Ykl_(k,l)>=m_Tk_Rl(k,l)) ? 1 : 0;
        }
      }
      //m_Akl_ = m_Ykl_>=(v_Tk_*v_Rl_.transpose()/2);
      if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_)
      {
        Epsilon_ = W1_/dimprod_;
        break;
      }
    }
    Epsilon_ = W1_/dimprod_;
  return true;
}

bool BinaryLBModelequalepsilon::initEMRows()
{
  // Initialization of various parameters
//int cols=std::min(100, nbVar_);
  int cols =  nbVar_;
  m_Uil_.resize(nbSample_,cols) = 0;
  selectRandomColsFromData(m_Uil_,cols);
  v_Ui_ = STK::sumByRow(m_Uil_);
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  m_Akl_.resize(Mparam_.nbrowclust_,cols) = 0;
  generateRandomMean(m_Akl_);
  W1_ = RealMax;
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  std::pair<int,int> Label_pair;
  for(int i = 0;i< (int)knownLabelsRows_.size();i++){
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
    // check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In BinaryLBModelequalepsilon::InitEMRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return false;
    }

    // M-step
    W1_old_ = W1_;
    W1_ = (m_Tik_.prod(m_sumik)).sum();
    m_Ukl_ = m_Tik_.transpose()*m_Uil_;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for ( int l = 0; l < cols; ++l)
      {
        m_Akl_(k,l)  = (m_Ukl_(k,l)>=v_Tk_[k]/2) ? 1 : 0;
      }
    }
    //m_Akl_ = m_Ukl_ >=(v_Tk_*MatrixReal(1,cols)/2);
    if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_)
    { break;}
  }
  return true;
}

bool BinaryLBModelequalepsilon::initEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  MatrixReal m_Tk_Rl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixReal m_Ulk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  m_Vjk_=m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;
  SelectRandomRows(m_Ulk);
  m_Ukl_ = m_Ulk.transpose();
  m_Akl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      m_Akl_(k,l) = (m_Ukl_(k,l)>=v_Tk_[k]/2) ? 1: 0;
    }
  }

  for ( int j=0;j< (int)knownLabelsCols_.size();j++ )
  {
    std::pair<int,int> Label_pair = knownLabelsCols_[j];
    m_Rjl_(Label_pair.first,Label_pair.second)=1;
  }
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    // E-step
    initBernoulliLogSumCols(m_sumjl);
    m_Rjl_ = (m_sumjl-STK::maxByRow(m_sumjl)*STK::Const::PointX(Mparam_.nbcolclust_)).exp();
    m_Rjl_ /= (STK::sumByRow(m_Rjl_)*STK::Const::PointX(Mparam_.nbcolclust_));
    //
    for ( int j=0;j< (int)knownLabelsCols_.size();j++)
    {
      m_Rjl_.row(knownLabelsCols_[j].first).setZeros();
      m_Rjl_(knownLabelsCols_[j].first,knownLabelsCols_[j].second)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepCols()) )
    {
      Error_msg_  = "In BinaryLBModelequalepsilon::InitEMCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return false;
    }

    // M-step
    W1_old_ = W1_;
    W1_ = (m_Rjl_.prod(m_sumjl)).sum();
    m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
    m_Tk_Rl = v_Tk_*v_Rl_.transpose()/2.0;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
      {
        m_Akl_(k,l) = (m_Ykl_(k,l)>=m_Tk_Rl(k,l)) ? 1 : 0;
      }
    }
    //m_Akl_ = m_Ykl_>=(v_Tk_*v_Rl_.transpose()/2);
    if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_)
    {
      Epsilon_ = W1_/dimprod_;
      break;
    }
  }
  Epsilon_ = W1_/dimprod_;
  return true;
}

void BinaryLBModelequalepsilon::selectRandomColsFromData(MatrixReal& _m_il,int cols)
{
  if(cols==Mparam_.nbcoldata_)
    _m_il=m_Dataij_.cast<STK::Real>();
  else{
    //random shuffle Algorithm
    VectorInt _v_temp = randSample(nbVar_,cols);

    for ( int l = 0; l < cols; ++l)
    {
      _m_il.col(l)=m_Dataij_.cast<STK::Real>().col(_v_temp[l]);
      //_m_il.col(l)=m_Dataij_.cast<STK::Real>().col(l);
    }
  }
}

void BinaryLBModelequalepsilon::SelectRandomRows(MatrixReal& m_lk)
{
  VectorInt v_temp = randSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    m_lk.row(l) = m_Vjk_.row(v_temp[l]);
    //m_lk.row(l) = m_Vjk_.row(l);
  }
}

void BinaryLBModelequalepsilon::generateRandomMean(MatrixBinary& m_kl)
{
  VectorInt _v_temp = randSample(nbSample_,Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    m_kl.row(k) = m_Uil_.cast<bool>().row(_v_temp[k]);
    //m_kl.row(k) = m_Uil_.cast<bool>().row(k);
  }
}

void BinaryLBModelequalepsilon::modifyThetaStart()
{
  m_Aklstart_ = m_Akl_;
  Epsilonstart_ = Epsilon_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;
  m_Rjlstart_ = m_Rjl_;
}

void BinaryLBModelequalepsilon::copyThetaStart()
{
  m_Akl_ = m_Aklstart_;
  Epsilon_ = Epsilonstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;
  m_Rjl_ = m_Rjlstart_;
  //initialization
  v_Rl_ = STK::sum(m_Rjl_);
  m_Ykl_old1_ = m_Ykl_;
}

void BinaryLBModelequalepsilon::copyThetaMax()
{
  m_Akl_ = m_Aklmax_;
  Epsilon_ = Epsilonmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;
  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  likelihood_ = Lmax_;
}

void BinaryLBModelequalepsilon::modifyThetaMax()
{
  m_Aklmax_ = m_Akl_;
  Epsilonmax_ = Epsilon_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;
  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = likelihood_;
}

void BinaryLBModelequalepsilon::mStepFull()
{
  if(!Mparam_.fixedproportions_)
  {
    v_logRhol_=(v_Rl_/nbVar_).log();
    v_logPiek_=(v_Tk_/nbSample_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Dataij_.cast<STK::Real>()*m_Rjl_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      m_Akl_(k,l) = (m_Ykl_(k,l)>=m_Tk_Rl_(k,l)) ? 1 : 0;
    }
  }
  Epsilon_= (m_Ykl_-(v_Tk_*v_Rl_.transpose()).prod(m_Akl_.cast<STK::Real>()) ).abs().sum()/dimprod_;

}

STK::Real BinaryLBModelequalepsilon::iclCriteriaValue(){
  STK::Real criteria = 0.0;

  criteria+= lgamma(Mparam_.nbrowclust_*a_)+lgamma(Mparam_.nbcolclust_*a_)
      -(Mparam_.nbrowclust_+Mparam_.nbcolclust_)*lgamma(a_)
      +Mparam_.nbrowclust_*Mparam_.nbcolclust_*(lgamma(2*b_)-2*lgamma(b_))
      -lgamma(Mparam_.nbrowdata_+Mparam_.nbrowclust_*a_)
      -lgamma(Mparam_.nbcoldata_+Mparam_.nbcolclust_*a_);

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    criteria+= lgamma(a_+ (v_Zi_== k).count());
  }

  for (int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    criteria+= lgamma(a_+ (v_Wj_==l).count());
  }

  STK::ArrayXXi temp0(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  STK::ArrayXXi temp1(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixBinary m_tempdata = (m_Dataij_==0);
  temp0 = (m_Zik_.transpose()*m_tempdata.cast<int>()*m_Wjl_)+b_;
  temp1 = (m_Zik_.transpose()*m_Dataij_.cast<int>()*m_Wjl_)+b_;
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      criteria+=lgamma(temp0(k,l))+lgamma(temp1(k,l));
    }
  }

  for (int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      criteria-= lgamma(((v_Zi_== k).count())*((v_Wj_==l).count())+2*b_);
    }
  }
  return criteria;
}

