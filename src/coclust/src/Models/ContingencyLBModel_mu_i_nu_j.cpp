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

ContingencyLBModel_mu_i_nu_j::ContingencyLBModel_mu_i_nu_j(MatrixReal const& m_Dataij,VectorInt const & rowlabels,VectorInt const & collabels,
                                                           VectorReal const& v_Mui,VectorReal const& v_Nuj,ModelParameters const& Mparam)
                             : ICoClustModel(Mparam,rowlabels,collabels)
                             , m_Dataij_(m_Dataij),v_Mui_(v_Mui),v_Nuj_(v_Nuj)
{
  DataSum_ = m_Dataij.sum();
};

bool ContingencyLBModel_mu_i_nu_j::cemInitStep()
{
#ifdef COVERBOSE
  std::cout<<"ContingencyLBModel_mu_i_nu_j::cemInitStep. Initializing Model Parameters."<<std::endl;
#endif
  // log proportions in rows and columns
//v_logPiek_ = -std::log(Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
//v_logRhol_ = -std::log(Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));

//// block posterior probabilities for individuals and variables
//m_Tik_.resize(nbSample_,Mparam_.nbrowclust_);
//m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_);

//// initialize Tik and Rjl
//randomFuzzyTik();
//randomFuzzyRjl();
//// create initial parameters
//if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
//{ m_Ykl_   = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
//else
//{ m_Ykl_   = m_Tik_.transpose()*(m_Dataij_*m_Rjl_);}
//m_Gammakl_ = m_Ykl_/(v_Tk_* v_Rl_.transpose());
//#ifdef COVERBOSE
//consoleOut();
//#endif
//bool fixedprop = Mparam_.fixedproportions_;
//Mparam_.fixedproportions_ = true;

    // start cem iterations
    if (initCEMRows())
    {
  #ifdef COVERBOSE
      std::cout << "ContingencyLBModel_mu_i_nu_j::initCEMRows done with success."<<std::endl;
      consoleOut();
  #endif
      if (initCEMCols())
      {
          m_Gammakl1_ = m_Gammakl_;
          m_Gammakl1old_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
          m_Gammaklold_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
          m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
          m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
          v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
          v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
    
  #ifdef COVERBOSE
    std::cout<<"ContingencyLBModel_mu_i_nu_j::cemInitStep. Initialization done with success."<<std::endl;
    consoleOut();
  #endif
//      Mparam_.fixedproportions_ = fixedprop;
        return true;
      }
    }

//Mparam_.fixedproportions_ = fixedprop;
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::initCEMRows()
{
////Initialization
//m_Uil_ = m_Dataij_*m_Rjl_;
////Determine row partition using EM algorithm with equal proportions
//for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
//{
//  if(!eStepRows()) return false;
//  //M-step
//  m_Gammaklold_ = m_Gammakl_;
//  mStepRows();
//  //Termination check
//  if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.initepsilon_)
//  { break;}
//}
//return true;


//Temporary variables, reduce the number of column
//int cols=std::min(100, nbVar_);
  int cols = nbVar_;

  MatrixReal m_Akl(Mparam_.nbrowclust_,cols), m_Aklold(Mparam_.nbrowclust_,cols);
  MatrixReal m_sumik(Mparam_.nbrowdata_,Mparam_.nbrowclust_);

  //Model parameters
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;

  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  v_Rl_ = STK::Const::VectorX(cols);

  // Initializations. m_Uil_ will contain cols columns of the original data set
  m_Uil_.resize(nbSample_, cols) = 0;
  selectRandomColsFromData(m_Uil_, cols);
  v_Ui_ = STK::sumByRow(m_Uil_);
  //
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
    m_sumik = m_Uil_*m_Akl.transpose();
    for ( int i =0;i< UnknownLabelsRows_.size();i++)
    {
      int maxIndex;
      m_sumik.row(UnknownLabelsRows_[i]).maxElt(maxIndex);
      m_Tik_.row(UnknownLabelsRows_[i]).setZeros();
      m_Tik_(UnknownLabelsRows_[i], maxIndex)=1;
    }
    // check empty class
    if( (empty_cluster_ = finalizeStepRows()) )
    {
      Error_msg_  = "In ContingencyLBModel_mu_i_nu_j::initCEMRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
#endif
      return false;
    }
    // M-step
    m_Aklold = m_Akl;
    m_Akl = ( (m_Tik_.transpose()*m_Uil_)
        /( (m_Tik_.transpose()*v_Ui_)*STK::Const::PointX(cols))+RealMin).log(); // @suppress("No return")
    // check convergence
    if((((m_Akl-m_Aklold).abs()/m_Akl).sum())<Mparam_.initepsilon_) { break;}
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::initCEMCols()
{
  //Initializations
//m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
////Determine column partition using EM algorithm with equal proportions
//for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
//{
//  if(!eStepCols()) return false;
//  //M-step
//  m_Gammaklold_ = m_Gammakl_;
//  mStepCols();
//  //Termination check
//  if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.initepsilon_)
//  { break;}
//}
//m_Gammakl1old_ = m_Gammakl1_;
//m_Gammakl1_    = m_Gammakl_;
//return true;

    //Temporary variables
    MatrixReal m_Alk(Mparam_.nbcolclust_,Mparam_.nbrowclust_)
             , m_Alkold(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
    MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  
    //Initializations
    m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
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
        int maxIndex;
        m_sumjl.row(UnknownLabelsCols_[j]).maxElt(maxIndex);
        m_Rjl_.row(UnknownLabelsCols_[j]).setZeros();
        m_Rjl_(UnknownLabelsCols_[j],maxIndex)=1;
      }
      // compute v_Rl_ and check empty class
      if( (empty_cluster_ = finalizeStepCols()) )
      {
        Error_msg_  = "In ContingencyLBModel_mu_i_nu_j::initCEMCols(). Class size too small while estimating model.\n";
  #ifdef COVERBOSE
        std::cout << Error_msg_;
  #endif
        return false;
      }
      // M-step
      m_Alkold = m_Alk;
      m_Alk = ( (m_Rjl_.transpose()*m_Vjk_)
               /( ( m_Rjl_.transpose()*v_Vj_)*STK::Const::PointX(Mparam_.nbrowclust_))
                  + RealMin).log();
      if((((m_Alk-m_Alkold).abs()/m_Alk).sum())<Mparam_.initepsilon_)
      { break;}
    }
    // try some optimization
    if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
    { m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
    else
    { m_Gammakl_ = (m_Tik_.transpose()* (m_Dataij_*m_Rjl_));}
    m_Gammakl_ /= STK::sumByRow(m_Gammakl_)*STK::sum(m_Gammakl_);
    return true;
}

bool ContingencyLBModel_mu_i_nu_j::emInitStep()
{
#ifdef COVERBOSE
  std::cout<<"ContingencyLBModel_mu_i_nu_j::emInitStep. Initializing Model Parameters."<<std::endl;
#endif
  // log proportions in rows and columns
//v_logPiek_ = -std::log(Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
//v_logRhol_ = -std::log(Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));

//// block posterior probabilities for individuals and variables
//m_Tik_.resize(nbSample_,Mparam_.nbrowclust_);
//m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_);

//// initialize Tik and Rjl
//randomFuzzyTik();
//randomFuzzyRjl();
//// create initial parameters
//if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
//{ m_Ykl_   = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
//else
//{ m_Ykl_   = m_Tik_.transpose()*(m_Dataij_*m_Rjl_);}
//m_Gammakl_ = m_Ykl_/(v_Tk_* v_Rl_.transpose());
//#ifdef COVERBOSE
//consoleOut();
//#endif
//bool fixedprop = Mparam_.fixedproportions_;
//Mparam_.fixedproportions_ = true;

    // start cem iterations
    if (initEMRows())
    {
  #ifdef COVERBOSE
      std::cout << "ContingencyLBModel_mu_i_nu_j::initEMRows done with success."<<std::endl;
      consoleOut();
  #endif
      if (initEMCols())
      {
          m_Gammakl1_ = m_Gammakl_;
          m_Gammakl1old_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
          m_Gammaklold_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_) = 0;
          m_Uil_.resize(nbSample_,Mparam_.nbcolclust_) = 0;
          m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
          v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
          v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));
    
  #ifdef COVERBOSE
    std::cout<<"ContingencyLBModel_mu_i_nu_j::emInitStep. Initialization done with success."<<std::endl;
    consoleOut();
  #endif
//      Mparam_.fixedproportions_ = fixedprop;
        return true;
      }
    }

//Mparam_.fixedproportions_ = fixedprop;
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::initEMRows()
{
////Initialization
//m_Uil_ = m_Dataij_*m_Rjl_;
////Determine row partition using EM algorithm with equal proportions
//for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
//{
//  if(!eStepRows()) return false;
//  //M-step
//  m_Gammaklold_ = m_Gammakl_;
//  mStepRows();
//  //Termination check
//  if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.initepsilon_)
//  { break;}
//}
//return true;


//Temporary variables, reduce the number of column
//int cols=std::min(100, nbVar_);
  int cols = nbVar_;

  MatrixReal m_Akl(Mparam_.nbrowclust_,cols), m_Aklold(Mparam_.nbrowclust_,cols);
  MatrixReal m_sumik(Mparam_.nbrowdata_,Mparam_.nbrowclust_);

  //Model parameters
  m_Tik_.resize(nbSample_,Mparam_.nbrowclust_) = 0;
  m_Vjk_.resize(nbVar_,Mparam_.nbrowclust_) = 0;
  v_Tk_.resize(Mparam_.nbrowclust_) = 0;

  m_Rjl_.resize(nbVar_,Mparam_.nbcolclust_) = 0;
  v_Rl_ = STK::Const::VectorX(cols);

  // Initializations. m_Uil_ will contain cols columns of the original data set
  m_Uil_.resize(nbSample_, cols) = 0;
  selectRandomColsFromData(m_Uil_, cols);
  v_Ui_ = STK::sumByRow(m_Uil_);
  //
  randomPoissonParameterRows(m_Akl,cols);

  for(int i=0;i<knownLabelsRows_.size();i++)
  {
    std::pair<int,int> Label_pair = knownLabelsRows_[i];
    m_Tik_(Label_pair.first,Label_pair.second) = 1;
  }

  //Determine row partition using EM algorithm with equal proportions
  for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
  {
    //E-Step
    m_sumik = m_Uil_*m_Akl.transpose();
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
      Error_msg_  = "In ContingencyLBModel_mu_i_nu_j::InitEMRows(). Class size too small while estimating model.\n";
#ifdef COVERBOSE
      std::cout << Error_msg_;
      std::cout << "v_Tk_= " << v_Tk_.transpose();
      std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
      return false;
    }
    // M-step
    m_Aklold = m_Akl;
    m_Akl = ( (m_Tik_.transpose()*m_Uil_)
        /( (m_Tik_.transpose()*v_Ui_)*STK::Const::PointX(cols))+RealMin).log();
    // check convergence
    if((((m_Akl-m_Aklold).abs()/m_Akl).sum())<Mparam_.initepsilon_) { break;}
  }
  return true;
}

bool ContingencyLBModel_mu_i_nu_j::initEMCols()
{
  //Initializations
//m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
////Determine column partition using EM algorithm with equal proportions
//for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
//{
//  if(!eStepCols()) return false;
//  //M-step
//  m_Gammaklold_ = m_Gammakl_;
//  mStepCols();
//  //Termination check
//  if((((m_Gammakl_-m_Gammaklold_)/m_Gammakl_).abs().sum())<Mparam_.initepsilon_)
//  { break;}
//}
//m_Gammakl1old_ = m_Gammakl1_;
//m_Gammakl1_    = m_Gammakl_;
//return true;

    //Temporary variables
    MatrixReal m_Alk(Mparam_.nbcolclust_,Mparam_.nbrowclust_)
             , m_Alkold(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
    MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  
    //Initializations
    m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
    v_Vj_  = STK::sumByRow(m_Vjk_);
    randomPoissonParameterCols(m_Alk);
    for ( int j = 0; j < knownLabelsCols_.size(); ++j)
    {
      std::pair<int,int> Label_pair = knownLabelsCols_[j];
      m_Rjl_(Label_pair.first,Label_pair.second)=1;
    }
    for ( int itr = 0; itr < Mparam_.nbinititerations_; ++itr)
    {
      //E-step
      m_sumjl = m_Vjk_*(m_Alk.transpose());
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
        Error_msg_  = "In ContingencyLBModel_mu_i_nu_j::InitEMCols(). Class size too small while running model.\n";
#ifdef COVERBOSE
        std::cout << Error_msg_;
        std::cout << "v_Tk_= " << v_Tk_.transpose();
        std::cout << "v_Rl_= " << v_Rl_.transpose();
#endif
        return false;
      }
      // M-step
      m_Alkold = m_Alk;
      m_Alk = ( (m_Rjl_.transpose()*m_Vjk_)
               /( ( m_Rjl_.transpose()*v_Vj_)*STK::Const::PointX(Mparam_.nbrowclust_))
                  + RealMin).log();
      if((((m_Alk-m_Alkold).abs()/m_Alk).sum())<Mparam_.initepsilon_)
      { break;}
    }
    // try some optimization
    if (m_Tik_.sizeCols() < m_Rjl_.sizeCols())
    { m_Gammakl_ = (m_Tik_.transpose()*m_Dataij_)*m_Rjl_;}
    else
    { m_Gammakl_ = (m_Tik_.transpose()* (m_Dataij_*m_Rjl_));}
    m_Gammakl_ /= STK::sumByRow(m_Gammakl_)*STK::sum(m_Gammakl_);
    return true;
}

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

bool ContingencyLBModel_mu_i_nu_j::GibbsRows()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

bool ContingencyLBModel_mu_i_nu_j::GibbsCols()
{
  Error_msg_ = "Gibbs is not implemented for this model.";
#ifdef COVERBOSE
  std::cout<<Error_msg_<<"\n";
#endif
  return false;
}

void ContingencyLBModel_mu_i_nu_j::selectRandomColsFromData(MatrixReal& _m_il,int cols)
{
  if(cols==nbVar_) {  _m_il = m_Dataij_;}
  else
  {
    //random shuffle Algorithm
    VectorInt _v_temp(Mparam_.nbcoldata_);
    for ( int j = 0; j < Mparam_.nbcoldata_; ++j) { _v_temp[j]=j;}
    for ( int l = 0; l < cols; ++l)
    {
//    int random=std::rand()%(Mparam_.nbcoldata_-l);	
		int random;
#ifdef _OPENMP
#pragma omp critical
#endif
      random = STK::Law::UniformDiscrete::rand(0, Mparam_.nbcoldata_-l - 1);
      int index  =_v_temp[random];
      _m_il.col(l)=m_Dataij_.col(index);
      //swap elements
      std::swap(_v_temp[Mparam_.nbcoldata_-l-1], _v_temp[random]);
    }
  }
}

void ContingencyLBModel_mu_i_nu_j::randomPoissonParameterRows(MatrixReal& m_kl,int cols)
{
  // sample nbrowclust_ integer form [1..nbSample]
  VectorInt _v_temp = randSample(nbSample_, Mparam_.nbrowclust_);
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    int index=_v_temp[k];
    for ( int l = 0; l < cols; ++l)
    {
      m_kl(k,l)=std::log(m_Uil_(index,l)/v_Ui_[index]+RealMin);
    }
  }
}

void ContingencyLBModel_mu_i_nu_j::randomPoissonParameterCols(MatrixReal& m_lk)
{
  VectorInt _v_temp = randSample(nbVar_,Mparam_.nbcolclust_);
  for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
  {
    int index=_v_temp[l];
    //index=l;
    for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
    {
      m_lk(l,k) = std::log(m_Vjk_(index,k)/v_Vj_[index]+RealMin);
    }
  }
}


void ContingencyLBModel_mu_i_nu_j::logSumRows(MatrixReal & m_ik)
{
  m_ik = STK::Const::VectorX(nbSample_)*v_logPiek_.transpose() +
      m_Uil_*((m_Gammakl_.log()).transpose())
      -v_Mui_*(m_Gammakl_*v_nul_).transpose();
}

void ContingencyLBModel_mu_i_nu_j::logSumCols(MatrixReal & m_jl)
{
  m_jl = STK::Const::VectorX(nbVar_)*v_logRhol_.transpose() + m_Vjk_*((m_Gammakl_.log()))
      -v_Nuj_*(v_muk_.transpose()*m_Gammakl_);
}

STK::Real ContingencyLBModel_mu_i_nu_j::estimateLikelihood()
{
  likelihood_ = (m_Ykl_.prod(m_Gammakl_.log()) ).sum() - DataSum_
		      + v_Tk_.transpose()*v_logPiek_
			  + v_Rl_.transpose()*v_logRhol_
	          - (m_Tik_.prod((RealMin + m_Tik_).log()) ).sum()
  		      - (m_Rjl_.prod((RealMin + m_Rjl_).log()) ).sum();

  return likelihood_;
}

/* @return the number of free parameters of the distribution of a block.*/
int ContingencyLBModel_mu_i_nu_j::nbFreeParameters() const
{return Mparam_.nbcolclust_*Mparam_.nbrowclust_;}

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
  v_logPiek_ = std::log(1.0/Mparam_.nbrowclust_)*(STK::Const::VectorX(Mparam_.nbrowclust_));
  v_logRhol_ = std::log(1.0/Mparam_.nbcolclust_)*(STK::Const::VectorX(Mparam_.nbcolclust_));

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
