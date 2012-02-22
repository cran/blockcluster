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

/** @file ContingencyLBModel.h
 *  @brief Declares concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#ifndef CONTINGENCYLBMODEL_H_
#define CONTINGENCYLBMODEL_H_

/** @brief Concrete model class for contingency data-sets.
 * This model assumes that the row and column effects are unknown.
 *
 */
#include "ICoClustModel.h"

class ContingencyLBModel: public ICoClustModel
{
  public:
    ContingencyLBModel( MatrixInteger const& m_Dataij,ModelParameters const& Mparam);
    ContingencyLBModel(MatrixInteger const& m_Dataij,VectorInteger const & rowlabels,
                       VectorInteger const & collabels,ModelParameters const& Mparam);
    virtual ~ContingencyLBModel(){};
    virtual ContingencyLBModel* Clone(){return new ContingencyLBModel(*this);}
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
    const MatrixInteger& GetArrangedDataClusters();
    /**Return Poisson Parameters ContingencyLBModel::m_Gammakl_*/
    const MatrixReal& GetGamma() const;
    /**Return Column proportions ContingencyLBModel::v_Piek_*/
  protected:
    //Variables involved in Poisson model
    MatrixInteger const& m_Dataij_;
    MatrixInteger m_ClusterDataij_;
    float DataSum_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Gammakl_, m_Gammaklold_, m_Gammakl1_, m_Gammakl1old_,m_Gammaklstart_,m_Gammaklmax_;
    VectorReal v_nul_;
    VectorReal v_muk_;
    VectorReal v_Ui_,v_Vj_;
    VectorReal v_Yl_,v_Yk_;
    MatrixReal m_Ykl_;
    float Likelihood_old;

    //M-steps
    void MStepRows();
    void MStepCols();

    // Functions to be used for internal computation in model initialization
    bool InitCEMRows();
    bool InitCEMCols();
    void SelectRandomColsfromdata(MatrixReal& m,int col);
    void GenerateRandomPoissonParameterRows(MatrixReal& m,int col);
    void GenerateRandomPoissonParameterCols(MatrixReal& m);
};

inline const MatrixReal& ContingencyLBModel::GetGamma() const
{
  return m_Gammakl_;
}

inline void ContingencyLBModel::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Ykl_ = m_Tik_.transpose()*m_Uil_;
  m_Gammakl_ = m_Ykl_.array()/(m_Ykl_.rowwise().sum()*v_Yl_.transpose()).array();
}

inline void ContingencyLBModel::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  m_Ykl_ = m_Vjk_.transpose()*m_Rjl_;
  m_Gammakl_ = m_Ykl_.array()/(v_Yk_*m_Ykl_.colwise().sum()).array();
}


#endif /* CONTINGENCYLBMODEL_H_ */
