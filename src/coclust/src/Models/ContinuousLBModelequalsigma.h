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
    virtual bool Estep();
    virtual bool CEstep();
    virtual void Mstep();
    virtual void likelihoodStopCriteria();
    virtual void ParameterStopCriteria();
    float EstimateLikelihood();
    virtual bool CEMInit();
    virtual bool FuzzyCEMInit();
    virtual bool RandomInit();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    MatrixReal GetArrangedDataClusters();
    virtual ~ContinuousLBModelequalsigma(){};

    /**Return means ContinuousLBModelequalsigma::m_Mukl_*/
    MatrixReal GetMean();
    /**Return Sigma ContinuousLBModelequalsigma::m_Sigma2kl_*/
    float GetSigma();
  protected:
    MatrixReal const& m_Dataij_;
    float dimprod_;
    MatrixReal m_Dataij2_;
    MatrixReal m_Mukl_,m_Mukl2_;
    float Sigma2_,Sigma2start_,Sigma2max_;
    MatrixReal m_Muklold1_,m_Muklold2_,m_Muklstart_,m_Muklmax_;
    MatrixReal m_Vjk1_,m_Vjk2_;
    MatrixReal m_Uil1_,m_Uil2_;
    float Likelihood_old;

  private:
    // Functions used to operate on data in intermediate steps when running the model
    bool EMRows();
    bool EMCols();
    bool CEMRows();
    bool CEMCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool InitCEMRows();
    void SelectRandomRowsfromdata(MatrixReal &);
    void GenerateRandomMean(const MatrixReal & , MatrixReal &);
};

inline MatrixReal ContinuousLBModelequalsigma::GetMean()
{
  return m_Mukl_;
}

inline float ContinuousLBModelequalsigma::GetSigma()
{
  return Sigma2_;
}

inline void ContinuousLBModelequalsigma::Modify_theta_start()
{
	m_Muklstart_ = m_Mukl_;
	Sigma2start_ = Sigma2_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

inline void ContinuousLBModelequalsigma::Copy_theta_start()
{
	m_Mukl_ = m_Muklstart_;
	Sigma2_ = Sigma2start_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
}

inline void ContinuousLBModelequalsigma::Copy_theta_max()
{
	m_Mukl_ = m_Muklmax_;
	Sigma2_ = Sigma2max_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

inline void ContinuousLBModelequalsigma::Modify_theta_max()
{
	m_Muklmax_ = m_Mukl_;
	Sigma2max_ = Sigma2_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}
#endif /* CONTINUOUSLBMODELEQUALSIGMA_H_ */
