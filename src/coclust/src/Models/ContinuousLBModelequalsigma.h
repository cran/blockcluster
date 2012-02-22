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


/** @file ContinuousLBModelequalsigma.h
 *  @brief Declares concrete model class ContinuousLBModelequalsigma derived from ICoClustModel.
 **/

#ifndef CONTINUOUSLBMODELEQUALSIGMA_H_
#define CONTINUOUSLBMODELEQUALSIGMA_H_

/** @brief Concrete model class for continuous data.
 * This class does presume equal variance among various  co-clusters.
 */
#include "ICoClustModel.h"

class ContinuousLBModelequalsigma: public ICoClustModel
{
  public:
    ContinuousLBModelequalsigma( MatrixReal const& m_Dataij,ModelParameters const& Mparam);
    ContinuousLBModelequalsigma(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,
                       VectorInteger const & collabels,ModelParameters const& Mparam);
    /** cloning */
    virtual ContinuousLBModelequalsigma* Clone(){return new ContinuousLBModelequalsigma(*this);}
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
    virtual void ParameterStopCriteria();
    virtual float EstimateLikelihood();
    virtual bool CEMInit();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    //virtual void UpdateAllUsingConditionalProbabilities();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    const MatrixReal& GetArrangedDataClusters();
    virtual ~ContinuousLBModelequalsigma(){};

    /**Return means ContinuousLBModelequalsigma::m_Mukl_*/
    const MatrixReal& GetMean() const;
    /**Return Sigma ContinuousLBModelequalsigma::m_Sigma2kl_*/
    float GetSigma() const;
  protected:
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    float dimprod_;
    MatrixReal m_Dataij2_;
    MatrixReal m_Mukl_,m_Mukl2_;
    float Sigma2_,Sigma2start_,Sigma2max_;
    MatrixReal m_Muklold1_,m_Muklold2_,m_Muklstart_,m_Muklmax_;
    MatrixReal m_Vjk1_,m_Vjk2_;
    MatrixReal m_Uil1_,m_Uil2_;
    float Likelihood_old;

    //M-steps
    void MStepRows();
    void MStepCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool InitCEMCols();
    void SelectRandomRowsfromdata(MatrixReal &);
    void GenerateRandomMean(const MatrixReal & , MatrixReal &);
};

inline const MatrixReal& ContinuousLBModelequalsigma::GetMean() const
{
  return m_Mukl_;
}

inline float ContinuousLBModelequalsigma::GetSigma() const
{
  return Sigma2_;
}

inline void ContinuousLBModelequalsigma::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Mukl_ = (m_Tik_.transpose()*m_Uil1_).array()/(v_Tk_*(v_Rl_.transpose())).array();
  m_Mukl2_ = m_Mukl_.array().square();
  Sigma2_ = ((m_Tik_.transpose()*m_Uil2_).array().sum()-v_Tk_.transpose()*m_Mukl2_*v_Rl_)/dimprod_;
}

inline void ContinuousLBModelequalsigma::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  m_Mukl_ = (m_Vjk1_.transpose()*m_Rjl_).array()/(v_Tk_*(v_Rl_.transpose())).array();
  m_Mukl2_ = m_Mukl_.array().square();
  Sigma2_ = ((m_Vjk2_.transpose()*m_Rjl_).array().sum()-v_Tk_.transpose()*m_Mukl2_*v_Rl_)/dimprod_;
}
#endif /* CONTINUOUSLBMODELEQUALSIGMA_H_ */
