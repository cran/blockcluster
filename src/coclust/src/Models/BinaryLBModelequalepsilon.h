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
#include "../../CImg/CImg.h"
#endif
class BinaryLBModelequalepsilon : public ICoClustModel
{
  public:
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b bayesian hyperparameters
     **/
    BinaryLBModelequalepsilon( MatrixBinary const& m_Dataij
                             , ModelParameters const& Mparam
                             , int a=1, int b=1);
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param rowlabels various labels for rows (-1  for unknown row label)
     * @param collabels various labels for columns (-1 for unknown column label)
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b bayesian hyperparameters
     **/
    BinaryLBModelequalepsilon( MatrixBinary const&  m_Dataij
                             , VectorInteger const & rowlabels
                             , VectorInteger const & collabels
                             , ModelParameters const& Mparam
                             , int a=1, int b=1);

    /** cloning */
    virtual BinaryLBModelequalepsilon* clone(){return new BinaryLBModelequalepsilon(*this);}
    virtual void logSumRows(MatrixReal & _m_sum);
    virtual void logSumCols(MatrixReal & _m_sum);
    virtual void mStepFull();
    virtual bool emRows();
    virtual bool cemRows();
    virtual bool emCols();
    virtual bool cemCols();
    virtual bool semRows();
    virtual bool semCols();
    virtual bool GibbsRows();
    virtual bool GibbsCols();
    virtual STK::Real estimateLikelihood();
    virtual void parameterStopCriteria();
    virtual bool cemInitStep();
    virtual bool emInitStep();
    virtual void finalizeOutput();
    virtual void consoleOut();
    virtual STK::Real iclCriteriaValue();
    //virtual void UpdateAllUsingConditionalProbabilities();
    virtual void modifyThetaStart();
    virtual void copyThetaStart();
    virtual void copyThetaMax();
    virtual void modifyThetaMax();
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;
    /**Return class mean BinaryLBModelequalepsilon::m_akl_ for all the blocks (co-clusters)*/
    MatrixBinary const& mean() const;
    /**Return Class despersion BinaryLBModelequalepsilon::m_epsilonkl_ for all the blocks (co-clusters) */
    STK::Real dispersion() const;
    /**Destructor*/
    virtual ~BinaryLBModelequalepsilon(){};

    MatrixBinary const& arrangedDataClusters();
#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void displayCluster();
#endif

  protected:
    //Variables involved in Bernoulli model
    int a_,b_;//hyper-parameters
    MatrixBinary const&  m_Dataij_;
    MatrixBinary m_ClusterDataij_;
    MatrixReal m_Xjl_,m_Xik_,m_Tk_Rl_;
    STK::Real dimprod_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_,m_Ukl_;
    VectorReal v_Ui_;
    MatrixReal m_Ykl_,m_Ykl_old2_,m_Ykl_old1_;
    MatrixBinary m_Akl_,m_Aklstart_,m_Aklmax_;
    STK::Real Epsilon_,Epsilonstart_,Epsilonmax_;
    STK::Real W1_,W1_old_;

    //M-steps
    void mStepRows();
    void mStepCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool initCEMRows();
    bool initCEMCols();
    bool initEMRows();
    bool initEMCols();
    void initBernoulliLogSumRows(MatrixReal & m_sum);
    void initBernoulliLogSumCols(MatrixReal & m_sum);
    void selectRandomColsFromData(MatrixReal& m,int col);
    void SelectRandomRows(MatrixReal& m_lk);
    void generateRandomMean(MatrixBinary& m);
};


inline MatrixBinary const& BinaryLBModelequalepsilon::mean() const
{ return m_Akl_;}

inline STK::Real BinaryLBModelequalepsilon::dispersion() const
{
  return Epsilon_;
}

inline void BinaryLBModelequalepsilon::mStepRows()
{
  if(!Mparam_.fixedproportions_) { v_logPiek_=(v_Tk_/nbSample_).log();}

  m_Ykl_   = m_Tik_.transpose()*m_Uil_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      m_Akl_(k,l) = (m_Ykl_(k,l)>=m_Tk_Rl_(k,l)) ? 1 :0;
    }
  }
  Epsilon_= (m_Ykl_- (v_Tk_*v_Rl_.transpose()).prod(m_Akl_.cast<STK::Real>()) ).abs().sum()/dimprod_;

}

inline void BinaryLBModelequalepsilon::mStepCols()
{
  if(!Mparam_.fixedproportions_) { v_logRhol_=(v_Rl_/nbVar_).log();}

  m_Ykl_   = m_Vjk_.transpose()*m_Rjl_;
  m_Tk_Rl_ = v_Tk_*v_Rl_.transpose()/2.0;
  for ( int k = 0; k < Mparam_.nbrowclust_; ++k)
  {
    for ( int l = 0; l < Mparam_.nbcolclust_; ++l)
    {
      m_Akl_(k,l) = (m_Ykl_(k,l)>=m_Tk_Rl_(k,l)) ? 1 : 0;
    }
  }
  Epsilon_= (m_Ykl_- (v_Tk_*v_Rl_.transpose()).prod(m_Akl_.cast<STK::Real>()) ).abs().sum()/dimprod_;

}



#endif /* BINARYLBMODELEQUALEPSILON_H_ */
