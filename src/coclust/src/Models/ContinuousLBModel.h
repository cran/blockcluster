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


/** @file ContinuousLBModel.h
 *  @brief Declares concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#ifndef CONTINUOUSLBMODEL_H_
#define CONTINUOUSLBMODEL_H_

#include "ICoClustModel.h"

/** @brief Concrete model class for continuous data.
 *  This class does not presume equal variance among various  co-clusters.
 */
class ContinuousLBModel: public ICoClustModel
{
  public:
    ContinuousLBModel( MatrixReal const& m_Dataij,ModelParameters const& Mparam);
    ContinuousLBModel(MatrixReal const& m_Dataij,VectorInteger const & rowlabels,
                       VectorInteger const & collabels,ModelParameters const& Mparam);
    virtual ~ContinuousLBModel(){};
    virtual ContinuousLBModel* clone(){return new ContinuousLBModel(*this);}
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
    virtual bool cemInitStep();
    virtual bool emInitStep();
    virtual void finalizeOutput();
    virtual void consoleOut();
    virtual void modifyThetaStart();
    virtual void copyThetaStart();
    virtual void copyThetaMax();
    virtual void modifyThetaMax();
    MatrixReal const& arrangedDataClusters();
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;

    /**Return various co-clusters mean ContinuousLBModel::m_Mukl_*/
    MatrixReal const& mean() const;
    /**Return various co-clusters Sigma ContinuousLBModel::m_Sigma2kl_*/
    MatrixReal const& sigma2() const;

  protected:
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    STK::Real dimprod_;
    MatrixReal m_Dataij2_;
    MatrixReal m_Mukl_, m_Sigma2kl_, m_Muklstart_, m_Muklmax_,m_Sigma2klstart_,m_Sigma2klmax_;
    MatrixReal m_Muklold1_,m_Muklold2_;
    MatrixReal m_Vjk1_,m_Vjk2_;
    MatrixReal m_Uil1_,m_Uil2_;

    //M-steps
    void mStepRows();
    void mStepCols();
    //Internal steps during initialization
    bool initCEMCols();
    bool initEMCols();
    void selectRandomRowsFromData(MatrixReal &);
    void generateRandomMean(const MatrixReal & , MatrixReal &);
};

inline MatrixReal const& ContinuousLBModel::mean() const
{ return m_Mukl_;}

inline MatrixReal const& ContinuousLBModel::sigma2() const
{ return m_Sigma2kl_;}

inline void ContinuousLBModel::mStepRows()
{
  if(!Mparam_.fixedproportions_) { v_logRhol_=(v_Rl_/nbVar_).log();}

  MatrixReal m_trkl = v_Tk_*v_Rl_.transpose();
  m_Mukl_           = (m_Tik_.transpose()*m_Uil1_)/m_trkl;
  m_Sigma2kl_       = (m_Tik_.transpose()*m_Uil2_)/m_trkl - m_Mukl_.square();
}

inline void ContinuousLBModel::mStepCols()
{
  if(!Mparam_.fixedproportions_) { v_logRhol_=(v_Rl_/nbVar_).log();}

  MatrixReal m_trkl = v_Tk_*v_Rl_.transpose();
  m_Mukl_           = (m_Vjk1_.transpose()*m_Rjl_)/m_trkl;
  m_Sigma2kl_       = ((m_Vjk2_.transpose()*m_Rjl_)/m_trkl) - m_Mukl_.square();
}

#endif /* CONTINUOUSLBMODEL_H_ */
