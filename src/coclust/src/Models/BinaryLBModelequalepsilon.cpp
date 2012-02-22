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
 * created on: Mar 1, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file BinaryLBModelequalepsilon.cpp
 *  @brief Implements concrete BinaryLBModelequalepsilon model class derived from ICoClustModel.
 **/

#include "BinaryLBModelequalepsilon.h"

#ifndef RPACKAGE
using namespace cimg_library;
#endif
BinaryLBModelequalepsilon::BinaryLBModelequalepsilon( MatrixBinary const& m_Dataij,ModelParameters const& Mparam)
                           : ICoClustModel(Mparam)
                           , m_Dataij_(m_Dataij)
{
  nbSample_ = m_Dataij.rows();
  nbVar_ = m_Dataij.cols();
  Mparam_.nbrowdata_ = nbSample_;
  Mparam_.nbcoldata_ = nbVar_;
  dimprod_ = nbSample_*nbVar_;
  m_Xjl_ = m_Dataij.cast<float>().colwise().sum().transpose()*MatrixReal::Ones(1,Mparam_.nbcolclust_);
  m_Xik_ = m_Dataij.cast<float>().rowwise().sum()*MatrixReal::Ones(1,Mparam_.nbrowclust_);
};

bool BinaryLBModelequalepsilon::Estep()
{
  if (EMRows()) {
    if (EMCols()) {
      return true;
    }
  }
  return false;
}


bool BinaryLBModelequalepsilon::CEstep()
{
  if (CEMRows()) {
    if (CEMCols()) {
      return true;
    }
  }
  return false;
}
void BinaryLBModelequalepsilon::Mstep()
{
  // Update Ykl for outer loop
  m_Ykl_old1_ = m_Ykl_;
}

bool BinaryLBModelequalepsilon::FuzzyCEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Fuzzy CEM initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Fuzzy CEM initialization is not valid for this model.";
#endif
  return false;
}
bool BinaryLBModelequalepsilon::RandomInit()
{
#ifdef COVERBOSE
  std::cout<<"Random initialization is not valid for this model.\n";
#endif

#ifdef RPACKAGE
  R_errormsg = "Random initialization is not valid for this model.";
#endif
  return false;
}

bool BinaryLBModelequalepsilon::CEMInit()
{
#ifdef COVERBOSE
  std::cout<<"Initializing Model Parameters.."<<std::endl;
#endif

  if (InitCEMRows()) {
    if (InitCEMCols()) {
      v_Rl_ = m_Rjl_.colwise().sum();
      m_Ykl_old1_ = m_Ykl_;
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

void BinaryLBModelequalepsilon::FinalizeOutput()
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

void BinaryLBModelequalepsilon::ConsoleOut()
{
#ifndef RPACKAGE
  std::cout<<"Output Model parameter:"<<"\nakl:\n"<<m_Akl_<<"\nepsilon:\n"<<Epsilon_<<"\npiek: "<<
      v_Piek_.transpose()<<"\nRhol: "<<v_Rhol_.transpose()<<std::endl;
#endif
}


//Compute Bernoulli log-sum for all rows
void BinaryLBModelequalepsilon::BernoulliLogsumRows(MatrixReal & m_sum)
{
  float logepsilon = log(Epsilon_/(1-Epsilon_));
  m_sum = VectorReal::Ones(nbSample_)*(v_logPiek_+logepsilon*m_Akl_.cast<float>()*v_Rl_).transpose()
          -logepsilon*(2*m_Uil_*m_Akl_.cast<float>().transpose() + m_Xik_.cast<float>());
}

//Compute Bernoulli log-sum for all columns
void BinaryLBModelequalepsilon::BernoulliLogsumCols(MatrixReal & m_sum)
{
  float logepsilon = log(Epsilon_/(1-Epsilon_));
  m_sum = VectorReal::Ones(nbVar_)*(v_logRhol_.transpose()+logepsilon*v_Tk_.transpose()*m_Akl_.cast<float>())
          -logepsilon*(2*m_Vjk_*m_Akl_.cast<float>() + m_Xjl_.cast<float>());
}

//Run EM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::EMRows()
{
  //Initializations
  m_Ykl_ = (m_Akl_.cast<float>().array()*(1-Epsilon_)+(1-m_Akl_.cast<float>().array())*Epsilon_)*dimprod_;
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;

  //Temporary variables
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_),m_prodik(nbSample_,Mparam_.nbrowclust_);
  VectorReal v_sumikmax(nbSample_),v_sumi(nbSample_),Zsumk(Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesk = VectorReal::Ones(Mparam_.nbrowclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      //E-step
      BernoulliLogsumRows(m_sumik);
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
      m_Ykl_old2_ = m_Ykl_;
      m_Ykl_ = m_Tik_.transpose()*m_Uil_;
      m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
          if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
          m_Akl_(k,l) = 1;
          else
            m_Akl_(k,l) = 0;
        }
      }
      Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;
      if(((m_Ykl_-m_Ykl_old2_).array()/m_Ykl_.array()).abs().sum()<Mparam_.epsilon_int_)
      {
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

//Run CEM algorithm on data Matrix m_Uil_
bool BinaryLBModelequalepsilon::CEMRows()
{
  //Initializations
  m_Ykl_ = (m_Akl_.cast<float>().array()*(1-Epsilon_)+(1-m_Akl_.cast<float>().array())*Epsilon_)*dimprod_;
  m_Uil_=m_Dataij_.cast<float>()*m_Rjl_;
  MatrixReal  m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index maxIndex;
  //CEMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
        BernoulliLogsumRows(m_sumik);
        m_Tik_.setZero(nbSample_,Mparam_.nbrowclust_);
        for (int i = 0; i < nbSample_; ++i) {
          m_sumik.row(i).maxCoeff(&maxIndex);
          m_Tik_(i,maxIndex)=1;
        }
        v_Tk_ = m_Tik_.colwise().sum();
        if((v_Tk_.array()<.00001).any()){
          empty_cluster = true;
          throw "Row clustering failed while running model.";
        }
        else
        {empty_cluster = false;}
        //M-step
        if(!Mparam_.fixedproportions_) {
          v_logPiek_=(v_Tk_.array()/nbSample_).log();
        }
        m_Ykl_old2_ = m_Ykl_;
        m_Ykl_ = m_Tik_.transpose()*m_Uil_;
        m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
        for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
          for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
            if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
            m_Akl_(k,l) = 1;
            else
              m_Akl_(k,l) = 0;
          }
        }
        //m_Akl_ = m_Ykl_.array()>=(v_Tk_*v_Rl_.transpose()/2).array();
        Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;
        if(((m_Ykl_-m_Ykl_old2_).array()/m_Ykl_.array()).abs().sum()<Mparam_.epsilon_int_)
        {
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

// Run EM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::EMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;

  //Temporary variables
  MatrixReal m_sumjl(nbVar_,Mparam_.nbcolclust_),m_prodjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal v_sumjlmax(nbVar_),v_sumj(nbVar_),Wsuml(Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  VectorReal Onesl = VectorReal::Ones(Mparam_.nbcolclust_);
  //EMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
      BernoulliLogsumCols(m_sumjl);
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

      //M-step
      if(!Mparam_.fixedproportions_) {
        v_logPiek_=(v_Tk_.array()/nbSample_).log();
      }
      m_Ykl_old2_ = m_Ykl_;
      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
          if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
          m_Akl_(k,l) = 1;
          else
            m_Akl_(k,l) = 0;
        }
      }
      //m_Akl_ = m_Ykl_.array()>=(v_Tk_*v_Rl_.transpose()/2).array();
      Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;
      if(((m_Ykl_-m_Ykl_old2_).array()/m_Ykl_.array()).abs().sum()<Mparam_.epsilon_int_)
      {
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

// Run CEM algorithm on data matrix m_Vjk_
bool BinaryLBModelequalepsilon::CEMCols()
{
  //Initializations
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  MatrixReal  m_sumjl(nbVar_,Mparam_.nbcolclust_);
  VectorReal::Index maxIndex;
  //CEMAlgo begins
  try {
    for (int itr = 0; itr < Mparam_.nbiterations_int_; ++itr) {
        BernoulliLogsumCols(m_sumjl);
        m_Rjl_.setZero(nbVar_,Mparam_.nbcolclust_);
        for (int j = 0; j < nbVar_; ++j) {
          m_sumjl.row(j).maxCoeff(&maxIndex);
          m_Rjl_(j,maxIndex)=1;
        }
        v_Rl_ = m_Rjl_.colwise().sum();
        if((v_Rl_.array()<.00001).any()){
          empty_cluster = true;
          throw "Column clustering failed while runnin model.";
        }
        else
        {empty_cluster = false;}

        //M-step
        if(!Mparam_.fixedproportions_) {
          v_logPiek_=(v_Tk_.array()/nbSample_).log();
        }
        m_Ykl_old2_ = m_Ykl_;
        m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
        m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
        for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
          for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
            if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
            m_Akl_(k,l) = 1;
            else
              m_Akl_(k,l) = 0;
          }
        }
        //m_Akl_ = m_Ykl_.array()>=(v_Tk_*v_Rl_.transpose()/2).array();
        Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;
        if(((m_Ykl_-m_Ykl_old2_).array()/m_Ykl_.array()).abs().sum()<Mparam_.epsilon_int_)
        {
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

float BinaryLBModelequalepsilon::EstimateLikelihood()
{
	Likelihood_ = (dimprod_*(Epsilon_*std::log(Epsilon_/(1-Epsilon_)) + std::log(1-Epsilon_))
	          +v_Tk_.dot(v_logPiek_) + v_Rl_.dot(v_logRhol_)
	          -(m_Tik_.array()*(RealMin + m_Tik_.array()).log()).sum()
	          -(m_Rjl_.array()*(RealMin + m_Rjl_.array()).log()).sum())/dimprod_;

	return Likelihood_;
}
//Computer Fuzzy clustering criteria and set the terminate variable accordingly
void BinaryLBModelequalepsilon::likelihoodStopCriteria()
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

//Compute change in Alpha and set the terminate variable accordingly
void BinaryLBModelequalepsilon::ParameterStopCriteria()
{
  float relativechange = ((m_Ykl_-m_Ykl_old1_).array()/m_Ykl_old1_.array()).abs().sum();
  if(relativechange<Mparam_.epsilon_)
    StopAlgo = true;
  else
    StopAlgo = false;
}

MatrixBinary BinaryLBModelequalepsilon::GetArrangedDataClusters()
{
  Eigen::ArrayXi v_Zi = (GetRowClassificationVector()).array();
  Eigen::ArrayXi v_Wj = (GetColumnClassificationVector()).array();
  MatrixBinary m_clusterij = MatrixBinary::Zero(nbSample_,nbVar_);

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
  MatrixBinary temp = m_clusterij;

  for (int i = 0; i < nbSample_; ++i) {
    m_clusterij.row( rowincrement(v_Zi(i)) + nbindrows(v_Zi(i))) = temp.row(i);
    rowincrement(v_Zi(i))+=1;
  }

  return m_clusterij;

}
#ifndef RPACKAGE
void BinaryLBModelequalepsilon::DisplayCluster()
{
  CImg<unsigned char>  cluster(nbVar_,nbSample_,1,1,0);
  CImg<unsigned char>  data(nbVar_,nbSample_,1,1,0);

  MatrixBinary m_clusterij = GetArrangedDataClusters();

  // Assign value to images
  for (int i = 0; i < nbSample_; ++i) {
    for (int j = 0; j < nbVar_; ++j) {
      if(m_clusterij(i,j) == 0)
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

bool BinaryLBModelequalepsilon::InitCEMRows()
{
  // Initialization of various parameters
  int cols=std::min(100,int(nbVar_));
  m_Uil_ = MatrixReal::Zero(nbSample_,cols);
  SelectRandomColsfromdata(m_Uil_,cols);
  v_Ui_ = m_Uil_.rowwise().sum();
  m_Vjk_ =MatrixReal::Zero(nbVar_,Mparam_.nbrowclust_);
  v_Rl_ = MatrixReal::Ones(cols,1);
  v_Tk_ = MatrixReal::Zero(Mparam_.nbrowclust_,1);
  m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
  m_Rjl_ = MatrixReal::Zero(nbVar_,Mparam_.nbcolclust_);
  m_Akl_.setZero(Mparam_.nbrowclust_,cols);
  GenerateRandomMean(m_Akl_);
  W1_ = RealMax;
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal _m_sumik(nbSample_,Mparam_.nbrowclust_);
  VectorReal::Index _minIndex;
  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      _m_sumik = MatrixReal::Ones(nbSample_,cols)*m_Akl_.cast<float>().transpose()-2*m_Uil_*m_Akl_.cast<float>().transpose()+v_Ui_*MatrixReal::Ones(1,Mparam_.nbrowclust_);
      m_Tik_ = MatrixReal::Zero(nbSample_,Mparam_.nbrowclust_);
      for (int i = 0; i < nbSample_; ++i) {
        _m_sumik.row(i).minCoeff(&_minIndex);
        m_Tik_(i,_minIndex)=1;
      }

      v_Tk_ = m_Tik_.colwise().sum();
      if((v_Tk_.array()<.00001).any()){
        empty_cluster = true;
        throw "Row clustering failed inside Initialization.";
      }
      else
      {empty_cluster = false;}

      // M-step
      W1_old_ = W1_;
      W1_ = (m_Tik_.array()*_m_sumik.array()).sum();
      m_Ukl_ = m_Tik_.transpose()*m_Uil_;
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        for (int l = 0; l < cols; ++l) {
          if(m_Ukl_(k,l)>=v_Tk_(k)/2)
          m_Akl_(k,l) = 1;
          else
            m_Akl_(k,l) = 0;
        }
      }
      //m_Akl_ = m_Ukl_.array()>=(v_Tk_*MatrixReal(1,cols)/2).array();
      if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_) {
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

bool BinaryLBModelequalepsilon::InitCEMCols()
{
  //Determine row partition using CEM algorithm with equal proportions
  MatrixReal m_sumjl(nbVar_, Mparam_.nbcolclust_);
  MatrixReal m_Tk_Rl(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  MatrixReal m_Ulk(Mparam_.nbcolclust_,Mparam_.nbrowclust_);
  VectorReal::Index _minIndex;
  m_Vjk_=m_Dataij_.cast<float>().transpose()*m_Tik_;
  SelectRandomRows(m_Ulk);
  m_Ukl_ = m_Ulk.transpose();
  m_Akl_.resize(Mparam_.nbrowclust_,Mparam_.nbcolclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
      if(m_Ukl_(k,l)>=v_Tk_(k)/2)
      m_Akl_(k,l) = 1;
      else
        m_Akl_(k,l) = 0;
    }
  }

  try {
    for (int itr = 0; itr < Mparam_.nbinititerations_; ++itr) {
      // CE-step
      m_sumjl = MatrixReal::Ones(nbVar_,1)*v_Tk_.transpose()*m_Akl_.cast<float>() - 2*m_Vjk_*m_Akl_.cast<float>() +m_Xjl_.cast<float>();
      m_Rjl_.setZero(nbVar_,Mparam_.nbcolclust_);
      for (int j = 0; j < nbVar_; ++j) {
        m_sumjl.row(j).minCoeff(&_minIndex);
        m_Rjl_(j,_minIndex)=1;
      }
      v_Rl_ = m_Rjl_.colwise().sum();
      //Check for empty cluster
      if((v_Rl_.array()<.00001).any()){
        empty_cluster = true;
        throw "Column clustering failed inside Initialization.";
      }
      else
      {empty_cluster = false;}

      // M-step
      W1_old_ = W1_;
      W1_ = (m_Rjl_.array()*m_sumjl.array()).sum();
      m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
      m_Tk_Rl = v_Tk_*v_Rl_.transpose()/2.0;
      for (int k = 0; k < Mparam_.nbrowclust_; ++k) {
        for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
          if(m_Ykl_(k,l)>=m_Tk_Rl(k,l))
          m_Akl_(k,l) = 1;
          else
            m_Akl_(k,l) = 0;
        }
      }
      //m_Akl_ = m_Ykl_.array()>=(v_Tk_*v_Rl_.transpose()/2).array();

      if (std::abs((W1_-W1_old_)/W1_)<Mparam_.initepsilon_) {
        Epsilon_ = W1_/dimprod_;
        break;
      }
    }
    Epsilon_ = W1_/dimprod_;
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

void BinaryLBModelequalepsilon::SelectRandomColsfromdata(MatrixReal& _m_il,int cols)
{
  if(cols==Mparam_.nbcoldata_)
    _m_il=m_Dataij_.cast<float>();
  else{
    //random shuffle Algorithm
    VectorInteger _v_temp = RandSample(nbVar_,cols);

    for (int l = 0; l < cols; ++l){
      _m_il.col(l)=m_Dataij_.cast<float>().col(_v_temp(l));
      //_m_il.col(l)=m_Dataij_.cast<float>().col(l);
    }
  }
}

void BinaryLBModelequalepsilon::SelectRandomRows(MatrixReal& m_lk)
{
  VectorInteger v_temp = RandSample(nbVar_,Mparam_.nbcolclust_);
  for (int l = 0; l < Mparam_.nbcolclust_; ++l) {
    m_lk.row(l) = m_Vjk_.row(v_temp(l));
    //m_lk.row(l) = m_Vjk_.row(l);
  }
}

void BinaryLBModelequalepsilon::GenerateRandomMean(MatrixBinary& m_kl)
{
  VectorInteger _v_temp = RandSample(nbSample_,Mparam_.nbrowclust_);
  for (int k = 0; k < Mparam_.nbrowclust_; ++k){
    m_kl.row(k) = m_Uil_.cast<bool>().row(_v_temp(k));
    //m_kl.row(k) = m_Uil_.cast<bool>().row(k);
  }
}
