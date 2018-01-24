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


/** @file ICoClustModel.cpp
 *  @brief This file only initializes the static members of ICoClustModel.
 **/

#include "ICoClustModel.h"

#define MINSIZE 3

ICoClustModel::ICoClustModel( ModelParameters const& Mparam)
                            : nbSample_(Mparam.nbrowdata_)
                            , nbVar_(Mparam.nbcoldata_)
                            , Mparam_(Mparam)
                            , likelihood_(-STK::Arithmetic<STK::Real>::infinity())
                            , Lmax_(-STK::Arithmetic<STK::Real>::infinity())
                            , empty_cluster_(false)
                            , stopAlgo_(false)
							, dimprod_(nbSample_*nbVar_)
{
  v_nbRowClusterMembers_.resize(Mparam.nbrowclust_) = 0;
  v_nbColClusterMembers_.resize(Mparam.nbcolclust_) = 0;
  v_Zi_.resize(nbSample_) = 0;
  v_Wj_.resize(nbVar_) = 0;
  m_Zik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Wjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  for ( int i = 0; i < Mparam.nbrowdata_; ++i)
  { UnknownLabelsRows_.push_back(i);}

  for ( int j = 0; j < Mparam.nbcoldata_; ++j)
  { UnknownLabelsCols_.push_back(j);}
}

ICoClustModel::ICoClustModel( ModelParameters const& Mparam
                            , VectorInteger const& rowlabels
                            , VectorInteger const& collabels)
                            : nbSample_(Mparam.nbrowdata_)
                            , nbVar_(Mparam.nbcoldata_)
                            , Mparam_(Mparam)
                            , likelihood_(-STK::Arithmetic<STK::Real>::infinity())
                            , empty_cluster_(false)
                            , stopAlgo_(false)
							, dimprod_(nbSample_*nbVar_)
{
  v_nbRowClusterMembers_.resize(Mparam.nbrowclust_) = 0;
  v_nbColClusterMembers_.resize(Mparam.nbcolclust_) = 0;
  v_Zi_.resize(nbSample_) = 0;
  v_Wj_.resize(nbVar_) = 0;
  m_Zik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Wjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  setRowLabels(rowlabels);
  setColLabels(collabels);
}

/* sample k integer from {0,1,...,n-1} */
VectorInteger ICoClustModel::randSample(int n,int k)
{
  VectorInteger v_temp(n), v_randint(k);
  for ( int j = 0; j < n; ++j) { v_temp[j]=j;}
  for ( int lfirst = 0, lend = n; lfirst < k; ++lfirst)
  {
    // get ramdom index
    int irand;
#ifdef _OPENMP
#pragma omp critical
#endif
    irand = STK::Law::UniformDiscrete::rand(0,lend - 1);
    v_randint[lfirst] = v_temp[irand];
    //swap elements
    --lend;
    std::swap(v_temp[lend], v_temp[irand]);
  }
  return v_randint;
}

/*
 * Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
 * ICoClustModel::stopAlgo_ parameter to either true or false depending on whether the change
 * in Likelihood is less than Mparam_.epsilon_ or not respectively.
 */
void ICoClustModel::likelihoodStopCriteria()
{
  STK::Real likelihood_old = likelihood_;
  estimateLikelihood();
  if(std::abs(likelihood_-likelihood_old)<std::abs(likelihood_)*Mparam_.epsilon_)
  { stopAlgo_ = true;}
  else
  { stopAlgo_ = false;}
}

/* compute the vector v_Tk_ and check if the size block is not too small
 *  @return false if the size block is under the threshold, true otherwise
 **/
bool ICoClustModel::finalizeStepRows()
{
  // compute size of the blocks by column
  v_Tk_ = STK::sum(m_Tik_);
  // check empty class
  return (v_Tk_ * v_Rl_.transpose() < MINSIZE).any();
}
/* compute the vector v_Rl_ and check if the size block is not too small
 *  @return false if the size block is under the threshold, true otherwise
 **/
bool ICoClustModel::finalizeStepCols()
{
  // compute size of the blocks by column
  v_Rl_ = STK::sum(m_Rjl_);
  // check empty class
  return (v_Tk_ * v_Rl_.transpose() < MINSIZE).any();
}

bool ICoClustModel::eStepRows()
{
  //E-step, compute sumik = log(\psi)
  MatrixReal  m_sumik(nbSample_, Mparam_.nbrowclust_);
  logSumRows(m_sumik);
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
    Error_msg_  = "In ICoClustModel::eStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::eStepCols()
{
  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);
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
    Error_msg_  = "In ICoClustModel::eStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
    std::cout << "v_Tk_= " << v_Tk_.transpose();
    std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::ceStepRows()
{
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_);
  logSumRows(m_sumik);

  for ( int i =0; i< (int)UnknownLabelsRows_.size();i++)
  {
    int maxIndex;
    m_sumik.row(UnknownLabelsRows_[i]).maxElt(maxIndex);
    m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
    m_Tik_(UnknownLabelsRows_[i],maxIndex)=1;
  }
  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ICoClustModel::ceStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::ceStepCols()
{
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);
  // adjust
  for ( int j=0;j< (int)UnknownLabelsCols_.size();j++)
  {
    int maxIndex;
    m_sumjl.row(UnknownLabelsCols_[j]).maxElt(maxIndex);
    m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
    m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
  }
  // check empty class
  if( (empty_cluster_ = finalizeStepCols()) )
  {
    Error_msg_  = "In ICoClustModel::ceStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::seStepRows()
{
  //E-step: Calculate conditional row class probabilities
  MatrixReal m_sumik(nbSample_,Mparam_.nbrowclust_);
  logSumRows(m_sumik);

  m_Tik_ = (m_sumik-STK::maxByRow(m_sumik)*STK::Const::PointX(Mparam_.nbrowclust_)).exp();
  m_Tik_ /=(STK::sumByRow(m_Tik_)*STK::Const::PointX(Mparam_.nbrowclust_));
  //
  for ( int i=0;i< (int)knownLabelsRows_.size();i++)
  {
    m_Tik_.row(knownLabelsRows_[i].first).setZeros();
    m_Tik_(knownLabelsRows_[i].first,knownLabelsRows_[i].second)=1;
  }
  //S-step : generate Row class matrix m_Zik_ using m_Tik_ and copy back to m_Tik_
  rowClassMatrixDraw();
  m_Tik_ = m_Zik_.cast<STK::Real>();
  // check empty class
  if( (empty_cluster_ = finalizeStepRows()) )
  {
    Error_msg_  = "In ICoClustModel::seStepRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

bool ICoClustModel::seStepCols()
{
  //E-step: Calculation of conditional column class probabilities
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_);
  logSumCols(m_sumjl);
  m_Rjl_ =(m_sumjl-STK::maxByRow(m_sumjl)*STK::Const::PointX(Mparam_.nbcolclust_)).exp();
  m_Rjl_ /=(STK::sumByRow(m_Rjl_)*STK::Const::PointX(Mparam_.nbcolclust_));
  //
  for ( int j=0;j< (int)knownLabelsCols_.size();j++)
  {
    m_Rjl_.row(knownLabelsCols_[j].first).setZeros();
    m_Rjl_(knownLabelsCols_[j].first,knownLabelsCols_[j].second)=1;
  }
  //S-step: Draw column class matrix m_Wjl_ using m_Rjl_ and copy back to m_Rjl_
  colClassMatrixDraw();
  m_Rjl_ = m_Wjl_.cast<STK::Real>();
  // check empty class
  if( (empty_cluster_ = finalizeStepCols()) )
  {
    Error_msg_  = "In ICoClustModel::seStepCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
    std::cout << Error_msg_;
#endif
    return false;
  }
  return true;
}

void ICoClustModel::rowClassMatrixDraw()
{
  //take cumulative sum of probabilities
  MatrixReal m_Tiktemp = m_Tik_;
  for ( int k = 1; k < m_Tiktemp.sizeCols(); ++k)
  { m_Tiktemp.col(k) += m_Tiktemp.col(k-1);}

  //generate random numbers
  std::vector<STK::Real> randnumbers(nbSample_);
  for ( int i = 0; i < nbSample_; ++i)
  {
    //std::srand(i);
//  randnumbers[i] = STK::Real(std::rand())/STK::Real(RAND_MAX);
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
}

void ICoClustModel::colClassMatrixDraw()
{
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
}

bool ICoClustModel::cemInitStep()
{
  Error_msg_ = "CEM initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ICoClustModel::emInitStep()
{
  Error_msg_ = "EM initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ICoClustModel::randomInit()
{
  Error_msg_ = "Random initialization is not valid for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

STK::Real ICoClustModel::iclCriteriaValue()
{

  STK::Real criteria = 0.0;
  //  Error_msg_ = "ICL creteria is not yet implemented for this model.";
  criteria =  likelihood_ - std::log(Mparam_.nbrowdata_)*(Mparam_.nbrowclust_-1.)/2.
              - std::log(Mparam_.nbcoldata_)*(Mparam_.nbcolclust_-1.)/2.
              - std::log(Mparam_.nbcoldata_*Mparam_.nbrowdata_) * nbFreeParameters()/2.;
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return criteria;
}

/* generate random tik_ */
int ICoClustModel::randomFuzzyTik()
{
  m_Tik_.randUnif();
  v_Tk_ = sum(m_Tik_);
  m_Tik_ /= STK::Const::VectorX(nbSample_) * v_Tk_.transpose() ;
  return v_Tk_.minElt();
}

/* generate random Rjl_ */
int ICoClustModel::randomFuzzyRjl()
{
  m_Rjl_.randUnif();
  v_Rl_ = sum(m_Rjl_);
  m_Rjl_ /= STK::Const::VectorX(nbVar_) * v_Rl_.transpose();
  return v_Rl_.minElt();
}


VectorInteger ICoClustModel::partRnd(int n, VectorReal proba)
{
  int clusters = proba.sizeRows();
  VectorInteger v_Z(n, 0);
  VectorInteger v_Randperm = randSample(n,n);
  VectorInteger remainingclusters(n-clusters);

  for ( int ind = 0; ind < clusters; ++ind)
  { v_Z[v_Randperm[ind]] = ind+1;}
  //   remainingclusters = (clusters+1)*MatrixInteger::Ones(n-clusters,1)
  //- ((STK::Const::Array<STK::Real>(clusters,1)*(unifRnd(0,1,1,n-clusters))).array() < (cumSum(proba)*STK::Const::Array<STK::Real>(1,n-clusters)).array()).matrix().cast<int>().colwise().sum().transpose();
  remainingclusters = (clusters+1)*STK::Const::VectorXi(n-clusters)
  - STK::count(  (STK::Const::VectorX(clusters)*(unifRnd(0,1,n-clusters)))
              < (cumSum(proba)*STK::Const::PointX(n-clusters))
              );

  for ( int ind = clusters; ind < n; ++ind)
  {
    v_Z[v_Randperm[ind]] = (remainingclusters[ind-clusters]==(clusters+1))
                         ? clusters:remainingclusters[ind-clusters];
  }
  return v_Z;
}

VectorReal ICoClustModel::cumSum(VectorReal proba)
{
  int size = proba.sizeRows();
  VectorReal v_temp(size, 0);
  v_temp[0] = proba[0];
  for ( int itr = 1; itr < size; ++itr)
  { v_temp[itr] = v_temp[itr-1] + proba[itr];}
  return v_temp;
}

PointReal ICoClustModel::unifRnd(STK::Real a,STK::Real b, int col)
{
  PointReal m_temp(col);
  for ( int c = 0; c < col; ++c)
  {
//   m_temp[c] = (b-a)*(std::rand()/STK::Real(RAND_MAX)) + a;
#ifdef _OPENMP
#pragma omp critical
#endif
     m_temp[c] = (b-a)*STK::Law::Uniform::rand(0,1) + a;
  }
  return m_temp;
}

void ICoClustModel::setRowLabels(VectorInteger const& rowlabels)
{
  int cluster, length = rowlabels.size();
  for ( int i = 0; i < length; ++i)
  {
    cluster = rowlabels[i];
    if (cluster<0)
    {
      UnknownLabelsRows_.push_back(i);
    }
    else
    {
      knownLabelsRows_.push_back(std::pair<int,int>(i,cluster));
      v_nbRowClusterMembers_[cluster]+=v_nbRowClusterMembers_[cluster];
      v_Zi_[i] = 1;
      m_Zik_(i,cluster) = 1;
    }
  }
}

void ICoClustModel::setColLabels(VectorInteger const& collabels)
{
  int cluster , length = collabels.size();
  for ( int j = 0; j < length; ++j)
  {
    cluster = collabels[j];
    if (cluster<0)
    {
      UnknownLabelsCols_.push_back(j);
    }
    else
    {
      knownLabelsCols_.push_back(std::pair<int,int>(j,cluster));
      v_nbColClusterMembers_[cluster]+=v_nbColClusterMembers_[cluster];
      v_Wj_[j] = 1;
      m_Wjl_(j,cluster) = 1;
    }
  }
}

void ICoClustModel::commonFinalizeOutput()
{
  // Calculate row and column proportions
  if(!Mparam_.fixedproportions_)
  {
    v_Piek_ = v_logPiek_.exp();
    v_Rhol_ = v_logRhol_.exp();
  }
  else
  {
    v_Piek_ = (1.0/Mparam_.nbrowclust_)*STK::Const::Vector<STK::Real>(Mparam_.nbrowclust_);
    v_Rhol_ = (1.0/Mparam_.nbcolclust_)*STK::Const::Vector<STK::Real>(Mparam_.nbcolclust_);
  }

  int maxIndex;
  // Calculate classification vector for rows
  m_Zik_.setZeros();
  m_Wjl_.setZeros();
  for ( int i = 0; i < nbSample_; ++i)
  {
    m_Tik_.row(i).maxElt(maxIndex);
    v_Zi_[i] = maxIndex;
    m_Zik_(i,maxIndex) = 1;
  }

  // Calculate classification vector for columns
  for ( int j = 0; j < nbVar_; ++j)
  {
    m_Rjl_.row(j).maxElt(maxIndex);
    v_Wj_[j] = maxIndex;
    m_Wjl_(j,maxIndex) = 1;
  }
}

#undef MINSIZE
