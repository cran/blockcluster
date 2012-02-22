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

/** @file BinaryLBModelequalepsilon.h
 *  @brief Declares concrete BinaryLBModelequalepsilon model class derived from ICoClustModel.
 **/

#ifndef BINARYLBMODELEQUALEPSILON_H_
#define BINARYLBMODELEQUALEPSILON_H_

/** @brief Concrete model class for Binary data.
 * This model assumes equal dispersion among various  co-clusters.
 *
 */
#include <limits.h>
#include <iostream>
#include "../typedefs/typedef.h"
#include "../Models/ICoClustModel.h"
#ifndef RPACKAGE
#include "../../../CImg/CImg.h"
#endif
class BinaryLBModelequalepsilon : public ICoClustModel
{
  public:
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param Mparam A constant reference to various ModelParameters.
     * */
    BinaryLBModelequalepsilon(MatrixBinary const& m_Dataij,ModelParameters const& Mparam);
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param rowlabels various labels for rows (-1  for unknown row label)
     * @param collabels various labels for columns (-1 for unknown column label)
     * @param Mparam A constant reference to various ModelParameters.
     * */
    BinaryLBModelequalepsilon(MatrixBinary const& m_Dataij,VectorInteger const & rowlabels,
                              VectorInteger const & collabels,ModelParameters const& Mparam);

    /** cloning */
    virtual BinaryLBModelequalepsilon* Clone(){return new BinaryLBModelequalepsilon(*this);}
    virtual void LogSumRows(MatrixReal & _m_sum);
    virtual void LogSumCols(MatrixReal & _m_sum);
    virtual void MStepFull();
    virtual bool EMRows();
    virtual bool CEMRows();
    virtual bool EMCols();
    virtual bool CEMCols();
    virtual bool SEMRows();
    virtual bool SEMCols();
    virtual void likelihoodStopCriteria();
    virtual float EstimateLikelihood();
    virtual void ParameterStopCriteria();
    virtual bool CEMInit();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    //virtual void UpdateAllUsingConditionalProbabilities();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    /**Return class mean BinaryLBModelequalepsilon::m_akl_ for all the blocks (co-clusters)*/
    const MatrixBinary& Getmean() const;
    /**Return Class despersion BinaryLBModelequalepsilon::m_epsilonkl_ for all the blocks (co-clusters) */
    float Getdispersion() const;
    /**Destructor*/
    inline virtual ~BinaryLBModelequalepsilon(){};

    const MatrixBinary& GetArrangedDataClusters();
#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void DisplayCluster();
#endif

  protected:
    //Variables involved in Bernoulli model
    MatrixBinary const& m_Dataij_;
    MatrixBinary m_ClusterDataij_;
    MatrixReal m_Xjl_,m_Xik_,m_Tk_Rl_;
    float dimprod_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_,m_Ukl_;
    VectorReal v_Ui_;
    MatrixReal m_Ykl_,m_Ykl_old2_,m_Ykl_old1_;
    MatrixBinary m_Akl_,m_Aklstart_,m_Aklmax_;
    float Epsilon_,Epsilonstart_,Epsilonmax_;
    float Likelihood_old;
    float W1_,W1_old_;

    //M-steps
    void MStepRows();
    void MStepCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool InitCEMRows();
    bool InitCEMCols();
    void InitBernoulliLogsumRows(MatrixReal & m_sum);
    void InitBernoulliLogsumCols(MatrixReal & m_sum);
    void SelectRandomColsfromdata(MatrixReal& m,int col);
    void SelectRandomRows(MatrixReal& m_lk);
    void GenerateRandomMean(MatrixBinary& m);
};


inline const MatrixBinary& BinaryLBModelequalepsilon::Getmean() const
{
  return m_Akl_;
}

inline float BinaryLBModelequalepsilon::Getdispersion() const
{
  return Epsilon_;
}

inline void BinaryLBModelequalepsilon::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Uil_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l) {
      if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
      m_Akl_(k,l) = 1;
      else
        m_Akl_(k,l) = 0;
    }
  }
  Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;

}

inline void BinaryLBModelequalepsilon::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }
  m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k) {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l) {
      if(m_Ykl_(k,l)>=m_Tk_Rl_(k,l))
      m_Akl_(k,l) = 1;
      else
        m_Akl_(k,l) = 0;
    }
  }
  Epsilon_= (m_Ykl_.array()-(v_Tk_*v_Rl_.transpose()).array()*m_Akl_.cast<float>().array()).abs().sum()/dimprod_;

}



#endif /* BINARYLBMODELEQUALEPSILON_H_ */
