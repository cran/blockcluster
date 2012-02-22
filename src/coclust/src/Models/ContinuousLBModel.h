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


/** @file ContinuousLBModel.h
 *  @brief Declares concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#ifndef CONTINUOUSLBMODEL_H_
#define CONTINUOUSLBMODEL_H_

/** @brief Concrete model class for continuous data.
 * This class does not presume equal variance among various  co-clusters.
 */
#include "ICoClustModel.h"

class ContinuousLBModel: public ICoClustModel
{
  public:
    ContinuousLBModel( MatrixReal const& m_Dataij,ModelParameters const& Mparam);
    ContinuousLBModel(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,
                       VectorInteger const & collabels,ModelParameters const& Mparam);
    virtual ContinuousLBModel* Clone(){return new ContinuousLBModel(*this);}
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
    virtual ~ContinuousLBModel(){};

    /**Return various co-clusters mean ContinuousLBModel::m_Mukl_*/
    const MatrixReal& GetMean() const;
    /**Return various co-clusters Sigma ContinuousLBModel::m_Sigma2kl_*/
    const MatrixReal& GetSigma() const;
  protected:
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    float dimprod_;
    MatrixReal m_Dataij2_;
    MatrixReal m_Mukl_,m_Mukl2_,m_Sigma2kl_,m_Muklstart_,m_Muklmax_,m_Sigma2klstart_,m_Sigma2klmax_;
    MatrixReal m_Muklold1_,m_Muklold2_;
    MatrixReal m_Vjk1_,m_Vjk2_;
    MatrixReal m_Uil1_,m_Uil2_;
    float Likelihood_old;

    //M-steps
    void MStepRows();
    void MStepCols();
    //Internal steps during initialization
    bool InitCEMCols();
    void SelectRandomRowsfromdata(MatrixReal &);
    void GenerateRandomMean(const MatrixReal & , MatrixReal &);
};

inline const MatrixReal& ContinuousLBModel::GetMean() const
{
  return m_Mukl_;
}

inline const MatrixReal& ContinuousLBModel::GetSigma() const
{
  return m_Sigma2kl_;
}

inline void ContinuousLBModel::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  MatrixReal m_trkl = (v_Tk_*v_Rl_.transpose());
  m_Mukl_ = (m_Tik_.transpose()*m_Uil1_).array()/m_trkl.array();
  m_Mukl2_ = m_Mukl_.array().square();
  m_Sigma2kl_ = (m_Tik_.transpose()*m_Uil2_).array()/m_trkl.array() - m_Mukl2_.array();
}

inline void ContinuousLBModel::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  MatrixReal m_trkl = (v_Tk_*v_Rl_.transpose());
  m_Mukl_ = (m_Vjk1_.transpose()*m_Rjl_).array()/m_trkl.array();
  m_Mukl2_ = m_Mukl_.array().square();
  m_Sigma2kl_ = (m_Vjk2_.transpose()*m_Rjl_).array()/m_trkl.array() - m_Mukl2_.array();
}

#endif /* CONTINUOUSLBMODEL_H_ */
