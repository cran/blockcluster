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


/** @file BinaryLBModel.h
 *  @brief Declares concrete model class BinaryLBModel derived from ICoClustModel.
 **/

#ifndef BinaryLBModel_H_
#define BinaryLBModel_H_

/** @brief Concrete model class for binary data.
 * This class assumes unequal dispersion among various co-clusters.
 *
 */
#include <limits.h>
#include <iostream>
#include "../typedefs/typedef.h"
#include "../Models/ICoClustModel.h"
#ifndef RPACKAGE
#include "../../../CImg/CImg.h"
#endif
class BinaryLBModel : public ICoClustModel
{
  public:
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param Mparam A constant reference to various ModelParameters.
     * */
    BinaryLBModel(MatrixBinary const& m_Dataij,ModelParameters const& Mparam);
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param rowlabels various labels for rows (-1  for unknown row label)
     * @param collabels various labels for columns (-1 for unknown column label)
     * @param Mparam A constant reference to various ModelParameters.
     * */
    BinaryLBModel(MatrixBinary const& m_Dataij,VectorInteger const & rowlabels,
                  VectorInteger const & collabels,ModelParameters const& Mparam);

    /** cloning */
    virtual BinaryLBModel* Clone(){return new BinaryLBModel(*this);}
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
    virtual bool CEMInit();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    //virtual void UpdateAllUsingConditionalProbabilities();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    float EstimateLikelihood();
    const MatrixBinary& GetArrangedDataClusters();
    /**Return class mean BinaryLBModel::m_akl_ for all the blocks (co-clusters)*/
    const MatrixBinary& Getmean() const;
    /**Return Class despersion BinaryLBModel::m_epsilonkl_ for all the blocks (co-clusters) */
    const MatrixReal& Getdispersion() const;
    /**Destructor*/
    inline virtual ~BinaryLBModel(){};


#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void DisplayCluster();
#endif

  protected:
    //Variables involved in Bernoulli model
    MatrixBinary const& m_Dataij_;
    MatrixBinary m_ClusterDataij_;
    float dimprod_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Alphakl_, m_Alphaklold_, m_Alphakl1_, m_Alphakl1old_,m_Alphaklmax_,m_Alphaklstart_;
    MatrixBinary m_akl_;
    MatrixReal m_epsilonkl_,m_epsilonklstart_,m_epsilonklmax_;
    float Likelihood_old;


    //M-steps
    void MStepRows();
    void MStepCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool InitCEMRows();
    bool InitCEMCols();
    void InitBernoulliLogsumRows(MatrixReal & m_sum);
    void InitBernoulliLogsumCols(MatrixReal & m_sum);
    void SelectRandomColsfromdata(MatrixReal& m,int col);
    void GenerateRandomBernoulliParameterRows(MatrixReal& m,int col);
    void GenerateRandomBernoulliParameterCols(MatrixReal& m);
};


inline const MatrixBinary& BinaryLBModel::Getmean() const
{
  return m_akl_;
}

inline const MatrixReal& BinaryLBModel::Getdispersion() const
{
  return m_epsilonkl_;
}

inline void BinaryLBModel::MStepRows()
{
  if(!Mparam_.fixedproportions_) {
    v_logPiek_=(v_Tk_.array()/nbSample_).log();
  }

  m_Alphakl_ = ((m_Tik_.transpose())*m_Uil_).array()/(v_Tk_*(v_Rl_.transpose())).array();
}

inline void BinaryLBModel::MStepCols()
{
  if(!Mparam_.fixedproportions_) {
    v_logRhol_=(v_Rl_.array()/nbVar_).log();
  }

  m_Alphakl_ = (m_Vjk_.transpose()*m_Rjl_).array()/(v_Tk_*v_Rl_.transpose()).array();
}

#endif /* BinaryLBModel_H_ */
