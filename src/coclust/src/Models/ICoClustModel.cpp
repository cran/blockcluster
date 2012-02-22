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


/** @file ICoClustModel.cpp
 *  @brief This file only initializes the static members of ICoClustModel.
 **/

#include "ICoClustModel.h"


ICoClustModel::ICoClustModel(ModelParameters const& Mparam):Mparam_(Mparam)
{
  nbSample_ = Mparam.nbrowdata_;
  nbVar_ = Mparam.nbcoldata_;
  v_nbRowClusterMembers_ = VectorInteger::Zero(Mparam.nbrowclust_);
  v_nbColClusterMembers_ = VectorInteger::Zero(Mparam.nbcolclust_);
  v_Zi_ = VectorInteger::Zero(nbSample_);
  v_Wj_ = VectorInteger::Zero(nbVar_);
  m_Zik_ = MatrixInteger::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Wjl_ = MatrixInteger::Zero(nbVar_,Mparam_.nbcolclust_);
  for ( int i = 0; i < Mparam.nbrowdata_; ++i) {
    UnknownLabelsRows_.push_back(i);
  }

  for ( int j = 0; j < Mparam.nbcoldata_; ++j) {
    UnknownLabelsCols_.push_back(j);
  }
}

ICoClustModel::ICoClustModel(ModelParameters const& Mparam,VectorInteger const & rowlabels,
              VectorInteger const & collabels):  Mparam_(Mparam)
{
  nbSample_ = Mparam.nbrowdata_;
  nbVar_ = Mparam.nbcoldata_;
  v_nbRowClusterMembers_ = VectorInteger::Zero(Mparam.nbrowclust_);
  v_nbColClusterMembers_ = VectorInteger::Zero(Mparam.nbcolclust_);
  v_Zi_ = VectorInteger::Zero(nbSample_);
  v_Wj_ = VectorInteger::Zero(nbVar_);
  m_Zik_ = MatrixInteger::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Wjl_ = MatrixInteger::Zero(nbVar_,Mparam_.nbcolclust_);
  SetRowLabels(rowlabels);
  SetColLabels(collabels);
}

VectorInteger ICoClustModel::RandSample(int n,int k)
{
  //random shuffle Algorithm
  //srand(float(clock()));
  int random,temp;
  VectorInteger v_temp(n),v_randint(k);
  for ( int j = 0; j < n; ++j) {
    v_temp(j)=j;
  }

  for ( int l = 0; l < k; ++l){
    random=std::rand()%(n-l);
    v_randint(l) = v_temp(random);
    //swap elements
    temp = v_temp(n-l-1);
    v_temp(n-l-1)=v_temp(random);
    v_temp(random)=temp;
  }
  return v_randint;
}

bool ICoClustModel::ERows()
{
  //Temporary variables
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_);
  VectorInteger Zsumk(VectorInteger::Zero(Mparam_.nbrowclust_));
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);
  //E-step
  LogSumRows(m_sumik);
  for ( int i = 0; i < nbSample_; ++i) {
    v_sumikmax(i) = m_sumik.row(i).maxCoeff(&maxIndex);
    Zsumk(maxIndex)+=1;
  }

  if((Zsumk.array()+v_nbRowClusterMembers_.array()<0.00001).any())
  {
    Error_msg_  = "Row clustering failed while running model.";
    empty_cluster = true;
#ifdef COVERBOSE
    std::cout<<Error_msg_<<"\n";
#endif
    return false;
  }
  else
  {
    empty_cluster = false;
  }

  m_prodik=(m_sumik-v_sumikmax*Onesk.transpose()).array().exp();
  v_sumi = m_prodik.rowwise().sum();
  m_Tik_ = m_prodik.array()/(v_sumi*Onesk.transpose()).array();
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for (LaLabel_pair : knownLabelsRows_) {
    m_Tik_.row(LaLabel_pair.first).setZero();
    m_Tik_(LabLabel_pair.first,LabLabel_pair.second)=1;
  }
#else
  for ( int i=0;i<knownLabelsRows_.size();i++) {
    Label_pair = knownLabelsRows_[i];
    m_Tik_.row(Label_pair.first).setZero();
    m_Tik_(Label_pair.first,Label_pair.second)=1;
  }
#endif
  v_Tk_ = m_Tik_.colwise().sum();
  return true;
}

bool ICoClustModel::ECols()
{
  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_);
  VectorInteger Wsuml(VectorInteger::Zero(Mparam_.nbcolclust_));
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);
  LogSumCols(m_sumjl);
  Wsuml.setZero();
   for ( int j = 0; j < nbVar_; ++j) {
     v_sumjlmax(j) = m_sumjl.row(j).maxCoeff(&maxIndex);
     Wsuml(maxIndex)+=1;
   }

   //Check for empty cluster
   if((Wsuml.array()+v_nbColClusterMembers_.array()<.00001).any()){
     Error_msg_  = "Column clustering failed while running model.";
#ifdef COVERBOSE
    std::cout<<Error_msg_<<"\n";
#endif
    empty_cluster = true;
    return false;
   }
   else
   {empty_cluster = false;}


   m_prodjl=(m_sumjl-v_sumjlmax*Onesl.transpose()).array().exp();
   v_sumj = m_prodjl.rowwise().sum();
   m_Rjl_ = m_prodjl.array()/(v_sumj*Onesl.transpose()).array();
   std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
   for (Label_pair : knownLabelsCols_) {
     m_Rjl_.row(Label_pair.first).setZero();
     m_Rjl_(Label_pair.first,Label_pair.second)=1;
   }
#else
   for ( int j=0;j<knownLabelsCols_.size();j++) {
     Label_pair = knownLabelsCols_[j];
     m_Rjl_.row(Label_pair.first).setZero();
     m_Rjl_(Label_pair.first,Label_pair.second)=1;
   }
#endif
   v_Rl_ = m_Rjl_.colwise().sum();

   return true;
}

bool ICoClustModel::CERows()
{
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  LogSumRows(m_sumik);
#ifdef RANGEBASEDFORLOOP
  for ( int i : UnknownLabelsRows_) {
    m_sumik.row(i).maxCoeff(&maxIndex);
    m_Tik_.row(i).setZero();
    m_Tik_(i,maxIndex)=1;
  }
#else
  for ( int i =0; i<UnknownLabelsRows_.size();i++) {
    m_sumik.row(UnknownLabelsRows_[i]).maxCoeff(&maxIndex);
    m_Tik_.row(UnknownLabelsRows_[i]).setZero();
    m_Tik_(UnknownLabelsRows_[i],maxIndex)=1;
  }
#endif
  v_Tk_ = m_Tik_.colwise().sum();
  if((v_Tk_.array()<.00001).any()){
    Error_msg_  = "Row clustering failed while running model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  empty_cluster = true;
  return false;
  }
  else
  {empty_cluster = false;}

  return true;
}

bool ICoClustModel::CECols()
{
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  LogSumCols(m_sumjl);
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
  v_Rl_ = m_Rjl_.colwise().sum();
  if((v_Rl_.array()<.00001).any()){
    Error_msg_  = "Column clustering failed while running model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  empty_cluster = true;
  return false;
  }
  else
  {empty_cluster = false;}

  return true;
}

bool ICoClustModel::SERows()
{
  //Temporary variables
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);

  //E-step: Calculate conditional row class probabilities
  LogSumRows(m_sumik);
  for ( int i = 0; i < nbSample_; ++i) {
    v_sumikmax(i) = m_sumik.row(i).maxCoeff(&maxIndex);
  }

  m_prodik = (m_sumik-v_sumikmax*Onesk.transpose()).array().exp();
  v_sumi = m_prodik.rowwise().sum();
  m_Tik_ = m_prodik.array()/(v_sumi*Onesk.transpose()).array();
  std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
  for (Label_pair : knownLabelsRows_) {
    m_Tik_.row(Label_pair.first).setZero();
    m_Tik_(Label_pair.first,Label_pair.second)=1;
  }
#else
  for ( int i=0;i<knownLabelsRows_.size();i++) {
    Label_pair = knownLabelsRows_[i];
    m_Tik_.row(Label_pair.first).setZero();
    m_Tik_(Label_pair.first,Label_pair.second)=1;
  }
#endif
  //S-step : generate Row class matrix m_Zik_ using m_Tik_ and copy back to m_Tik_
  RowClassMatrixdraw();
  m_Tik_ = m_Zik_.cast<float>();
  v_Tk_ = m_Tik_.colwise().sum();

  if((v_Tk_.array()<0.00001).any())
  {
    Error_msg_  = "Row clustering failed while running model.";
#ifdef COVERBOSE
    std::cout<<Error_msg_<<"\n";
#endif
    empty_cluster = true;
    return false;
  }
  else
  {
    empty_cluster = false;
  }

  return true;
}

bool ICoClustModel::SECols()
{
  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_);
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);

  //E-step: Calculation of conditional column class probabilities

  LogSumCols(m_sumjl);
   for ( int j = 0; j < nbVar_; ++j) {
     v_sumjlmax(j) = m_sumjl.row(j).maxCoeff(&maxIndex);
   }

   m_prodjl=(m_sumjl-v_sumjlmax*Onesl.transpose()).array().exp();
   v_sumj = m_prodjl.rowwise().sum();
   m_Rjl_ = m_prodjl.array()/(v_sumj*Onesl.transpose()).array();
   std::pair<int,int> Label_pair;
#ifdef RANGEBASEDFORLOOP
   for (Label_pair : knownLabelsCols_) {
     m_Rjl_.row(Label_pair.first).setZero();
     m_Rjl_(Label_pair.first,Label_pair.second)=1;
   }
#else
   for ( int j=0;j<knownLabelsCols_.size();j++) {
     Label_pair = knownLabelsCols_[j];
     m_Rjl_.row(Label_pair.first).setZero();
     m_Rjl_(Label_pair.first,Label_pair.second)=1;
   }
#endif
   //S-step: Draw column class matrix m_Wjl_ using m_Rjl_ and copy back to m_Rjl_
   ColClassMatrixdraw();
   m_Rjl_ = m_Wjl_.cast<float>();

   v_Rl_ = m_Rjl_.colwise().sum();

   if((v_Rl_.array()<.00001).any()){
     Error_msg_  = "Column clustering failed while running model.";
#ifdef COVERBOSE
   std::cout<<Error_msg_<<"\n";
#endif
   empty_cluster = true;
   return false;
   }
   else
   {
     empty_cluster = false;
   }

   return true;
}

void ICoClustModel::RowClassMatrixdraw()
{
  //take cumulative sum of probabilities
  MatrixReal m_Tiktemp = m_Tik_;
  for ( int k = 1; k < m_Tiktemp.cols(); ++k) {
    m_Tiktemp.col(k) += m_Tiktemp.col(k-1);
  }

  //generate random numbers
  std::vector<float> randnumbers(nbSample_);
  for ( int i = 0; i < nbSample_; ++i) {
    //std::srand(i);
    randnumbers[i] = float(std::rand())/float(RAND_MAX);
  }

  m_Zik_.setZero();

  //chose randomly the row class using generated random numbers
  for ( int i = 0; i < nbSample_; ++i) {
    for ( int k = 0; k < m_Tiktemp.cols(); ++k) {
      if(randnumbers[i]< m_Tiktemp(i,k))
      {
        m_Zik_(i,k) = 1;
        break;
      }
    }
  }
}

void ICoClustModel::ColClassMatrixdraw()
{
  //take cumulative sum of probabilities
  MatrixReal m_Rjltemp = m_Rjl_;
  for ( int l = 1; l < m_Rjltemp.cols(); ++l) {
    m_Rjltemp.col(l) += m_Rjltemp.col(l-1);
  }

  //generate random numbers
  std::vector<float> randnumbers(nbSample_);
  for ( int j = 0; j < nbVar_; ++j) {
    //std::srand(j);
    randnumbers[j] = float(std::rand())/float(RAND_MAX);
  }

  m_Wjl_.setZero();

  //chose randomly the row class using generated random numbers
  for ( int j = 0; j < nbVar_; ++j) {
    for ( int l = 0; l < m_Rjltemp.cols(); ++l) {
      if(randnumbers[j]<m_Rjltemp(j,l))
      {
        m_Wjl_(j,l) = 1;
        break;
      }
    }
  }
}

bool ICoClustModel::CEMInit()
{

  Error_msg_ = "CEM initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ICoClustModel::FuzzyCEMInit()
{
  Error_msg_ = "Fuzzy CEM initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ICoClustModel::RandomInit()
{
  Error_msg_ = "Random initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

float ICoClustModel::ICLCriteriaValue(){
  Error_msg_ = "ICL creteria is not yet implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return -RealMax;
}
VectorInteger ICoClustModel::PartRnd(int n,VectorReal proba)
{
  int clusters = proba.rows();
  VectorInteger v_Z = VectorInteger::Zero(n);
  VectorInteger v_Randperm = RandSample(n,n);
  VectorInteger remainingclusters(n-clusters);

  for ( int ind = 0; ind < clusters; ++ind) {
    v_Z(v_Randperm(ind)) = ind+1;
  }
  remainingclusters = (clusters+1)*MatrixInteger::Ones(n-clusters,1) - ((MatrixReal::Ones(clusters,1)*(Unifrnd(0,1,1,n-clusters))).array() < (Cumsum(proba)*MatrixReal::Ones(1,n-clusters)).array()).matrix().cast<int>().colwise().sum().transpose();

  for ( int ind = clusters; ind < n; ++ind) {
    v_Z(v_Randperm(ind)) = remainingclusters(ind-clusters)==(clusters+1)?clusters:remainingclusters(ind-clusters);
  }
  return v_Z;

}

VectorReal ICoClustModel::Cumsum(VectorReal proba)
{
  int size = proba.rows();
  VectorReal v_temp = VectorReal::Zero(size);
  v_temp(0) = proba(0);
  for ( int itr = 1; itr < size; ++itr) {
    v_temp(itr) = v_temp(itr-1) + proba(itr);
  }
  return v_temp;
}

MatrixReal ICoClustModel::Unifrnd(float a,float b, int row, int col)
{
  MatrixReal m_temp(row,col);
  for ( int r = 0; r < row; ++r) {
    for ( int c = 0; c < col; ++c) {
      m_temp(r,c) = (b-a)*(std::rand()/float(RAND_MAX)) + a;
    }
  }
  return m_temp;
}

