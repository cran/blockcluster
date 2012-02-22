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
 * created on: Jan 30, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ContingencyLBModel_mu_i_nu_j.h
 *  @brief Declares concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#ifndef CONTINGENCYLBMODEL_MU_I_NU_J_H_
#define CONTINGENCYLBMODEL_MU_I_NU_J_H_

/** @brief Concrete model class for contingency data-sets.
 * This model assumes that the row and column effects are known.
 *
 */
#include "ICoClustModel.h"

class ContingencyLBModel_mu_i_nu_j: public ICoClustModel
{
  public:
    /**
     * Constructor
     * @param m_Dataij Constant reference to contingency data.
     * @param v_Mui Constant reference to Known row effects.
     * @param v_Nuj Constant reference to Known column effects.
     * @return
     */
    ContingencyLBModel_mu_i_nu_j(MatrixInteger const& m_Dataij,VectorReal const& v_Mui,VectorReal const& v_Nuj,ModelParameters const& Mparam);
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
    MatrixInteger GetArrangedDataClusters();
    /**Return Poisson Parameters ContingencyLBModel_mu_i_nu_j::m_Gammakl_*/
    MatrixReal GetGamma();
    /**Destructor*/
    inline virtual ~ContingencyLBModel_mu_i_nu_j(){};

  protected:
    //Variables involved in Bernoulli model
    MatrixInteger const& m_Dataij_;
    VectorReal const& v_Mui_;
    VectorReal const& v_Nuj_;
    float DataSum_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Gammakl_, m_Gammaklold_, m_Gammakl1_, m_Gammakl1old_,m_Gammaklstart_,m_Gammaklmax_;
    VectorReal v_nul_;
    VectorReal v_muk_;
    MatrixReal m_Ykl_;

    float Likelihood_old;

    // Functions to be used for internal computations
    void PoissonLogsumRows(MatrixReal & m);
    void PoissonLogsumCols(MatrixReal & m);
    bool EMRows();
    bool CEMRows();
    bool EMCols();
    bool CEMCols();


  private:
    //Utility functions
    VectorInteger PartRnd(int n,VectorReal proba);
    VectorReal Cumsum(VectorReal proba);
    MatrixReal Unifrnd(float a, float b, int row, int col);
};

inline MatrixReal ContingencyLBModel_mu_i_nu_j::GetGamma()
{
  return m_Gammakl_;
}

inline void ContingencyLBModel_mu_i_nu_j::Modify_theta_start()
{
	m_Gammaklstart_ = m_Gammakl_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

inline void ContingencyLBModel_mu_i_nu_j::Copy_theta_start()
{
	m_Gammakl_ = m_Gammaklstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
  m_Gammakl1_ = m_Gammakl_;
}

inline void ContingencyLBModel_mu_i_nu_j::Copy_theta_max()
{
	m_Gammakl_ = m_Gammaklmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

inline void ContingencyLBModel_mu_i_nu_j::Modify_theta_max()
{
	m_Gammaklmax_ = m_Gammakl_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}
#endif /* CONTINGENCYLBMODEL_MU_I_NU_J_H_ */
