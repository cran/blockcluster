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
    ContingencyLBModel_mu_i_nu_j( MatrixReal const& m_Dataij, VectorReal const& v_Mui,
                                  VectorReal const& v_Nuj, ModelParameters const& Mparam);
    ContingencyLBModel_mu_i_nu_j( MatrixReal const& m_Dataij, VectorInt const & rowlabels
                                , VectorInt const & collabels
                                , VectorReal const& v_Mui
                                , VectorReal const& v_Nuj
                                , ModelParameters const& Mparam);
    /**Destructor*/
    inline virtual ~ContingencyLBModel_mu_i_nu_j(){};

    /** cloning */
    virtual ContingencyLBModel_mu_i_nu_j* clone(){return new ContingencyLBModel_mu_i_nu_j(*this);}
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
    virtual void parameterStopCriteria();
    virtual STK::Real estimateLikelihood();
    virtual bool randomInit();
    virtual bool cemInitStep();
    virtual bool emInitStep();
    virtual void finalizeOutput();
    virtual void consoleOut();
    virtual void modifyThetaStart();
    virtual void copyThetaStart();
    virtual void copyThetaMax();
    virtual void modifyThetaMax();
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;
    MatrixReal const& arrangedDataClusters();
    /**Return Poisson Parameters ContingencyLBModel_mu_i_nu_j::m_Gammakl_*/
    MatrixReal const& gamma() const;

  protected:
    //Variables involved in Bernoulli model
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    VectorReal const& v_Mui_;
    VectorReal const& v_Nuj_;
    STK::Real DataSum_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Gammakl_, m_Gammaklold_, m_Gammakl1_, m_Gammakl1old_,m_Gammaklstart_,m_Gammaklmax_;
    VectorReal v_nul_;
    VectorReal v_muk_;
    VectorReal v_Ui_,v_Vj_;
    MatrixReal m_Ykl_;

    //M-steps
    void mStepRows();
    void mStepCols();
    // Functions to be used for internal computation in model initialization
    bool initCEMRows();
    bool initCEMCols();
    bool initEMRows();
    bool initEMCols();
    void selectRandomColsFromData(MatrixReal& m,int col);
    void randomPoissonParameterRows(MatrixReal& m,int col);
    void randomPoissonParameterCols(MatrixReal& m);
};

inline MatrixReal const& ContingencyLBModel_mu_i_nu_j::gamma() const
{ return m_Gammakl_;}

inline void ContingencyLBModel_mu_i_nu_j::mStepRows()
{
  if(!Mparam_.fixedproportions_) { v_logPiek_=(v_Tk_/nbSample_).log();}

  m_Ykl_     = m_Tik_.transpose()*m_Uil_;
  m_Gammakl_ = m_Ykl_/(m_Tik_.transpose()*v_Mui_*v_nul_.transpose());
  //m_Gammakl_ = m_Ykl_/(v_Tk_* v_Rl_.transpose());
}

inline void ContingencyLBModel_mu_i_nu_j::mStepCols()
{
  if(!Mparam_.fixedproportions_) { v_logRhol_=(v_Rl_/nbVar_).log();}

  m_Ykl_     = m_Vjk_.transpose()*m_Rjl_;
  m_Gammakl_ = m_Ykl_/(v_muk_*v_Nuj_.transpose()*m_Rjl_);
  //m_Gammakl_ = m_Ykl_/(v_Tk_* v_Rl_.transpose());
}

#endif /* CONTINGENCYLBMODEL_MU_I_NU_J_H_ */
