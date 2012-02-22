/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2011  Parmeet Singh Bhatia

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
 * Project:  coclust
 * created on: Nov 25, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

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
    /**Constructor
     * @param m_Dataij a constant reference on the data set
     * */
    BinaryLBModel(MatrixBinary const& m_Dataij,ModelParameters const& Mparam);

    virtual bool Estep();
    virtual bool CEstep();
    virtual void Mstep();
    virtual void likelihoodStopCriteria();
    virtual void ParameterStopCriteria();
    virtual bool CEMInit();
    virtual bool FuzzyCEMInit();
    virtual bool RandomInit();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    float EstimateLikelihood();
    MatrixBinary GetArrangedDataClusters();
    /**Return class mean BinaryLBModel::m_akl_ for all the blocks (co-clusters)*/
    MatrixBinary& Getmean();
    /**Return Class despersion BinaryLBModel::m_epsilonkl_ for all the blocks (co-clusters) */
    MatrixReal& Getdispersion();
    /**Destructor*/
    inline virtual ~BinaryLBModel(){};


#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void DisplayCluster();
#endif

  protected:
    //Variables involved in Bernoulli model
    MatrixBinary const& m_Dataij_;
    float dimprod_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Alphakl_, m_Alphaklold_, m_Alphakl1_, m_Alphakl1old_,m_Alphaklmax_,m_Alphaklstart_;
    MatrixBinary m_akl_;
    MatrixReal m_epsilonkl_,m_epsilonklstart_,m_epsilonklmax_;
    float Likelihood_old;

  private:
    // Functions used to operate on data in intermediate steps when running the model
    void BernoulliLogsumRows(MatrixReal & _m_sum);
    void BernoulliLogsumCols(MatrixReal & _m_sum);
    bool EMRows();
    bool CEMRows();
    bool EMCols();
    bool CEMCols();
    void Compute_FuzzyStoppingCriteria();
    void Compute_ChangeinAlpha();
    bool TerminateEEstep();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool InitCEMRows();
    bool InitCEMCols();
    void InitBernoulliLogsumRows(MatrixReal & _m_sum);
    void InitBernoulliLogsumCols(MatrixReal & _m_sum);
    void SelectRandomColsfromdata(MatrixReal& m,int col);
    void GenerateRandomBernoulliParameterRows(MatrixReal& m,int col);
    void GenerateRandomBernoulliParameterCols(MatrixReal& m);
};


inline MatrixBinary& BinaryLBModel::Getmean()
{
  return m_akl_;
}

inline MatrixReal& BinaryLBModel::Getdispersion()
{
  return m_epsilonkl_;
}

inline void BinaryLBModel::Modify_theta_start()
{
	m_Alphaklstart_ = m_Alphakl_;
	//m_aklstart_ = m_akl_;
	m_epsilonklstart_ = m_epsilonkl_;
	v_logPiekstart_ = v_logPiek_;
	v_logRholstart_ = v_logRhol_;

	m_Rjlstart_ = m_Rjl_;
}

inline void BinaryLBModel::Copy_theta_start()
{
	m_Alphakl_ = m_Alphaklstart_;
	//m_akl_ = m_aklstart_;
	m_epsilonkl_ = m_epsilonklstart_;
	v_logPiek_ = v_logPiekstart_;
	v_logRhol_ = v_logRholstart_;

	m_Rjl_ = m_Rjlstart_;

	//initialization
	v_Rl_ = m_Rjl_.colwise().sum();
	m_Alphakl1_ = m_Alphakl_;
}

inline void BinaryLBModel::Copy_theta_max()
{
	m_Alphakl_ = m_Alphaklmax_;
	//m_akl_ = m_aklmax_;
	m_epsilonkl_ = m_epsilonklmax_;
	v_logPiek_ = v_logPiekmax_;
	v_logRhol_ = v_logRholmax_;

	m_Tik_ = m_Tikmax_;
	m_Rjl_ = m_Rjlmax_;
	Likelihood_ = Lmax_;
}
inline void BinaryLBModel::Modify_theta_max()
{
	m_Alphaklmax_ = m_Alphakl_;
	//m_aklmax_ = m_akl_;
	m_epsilonklmax_ = m_epsilonkl_;
	v_logPiekmax_ = v_logPiek_;
	v_logRholmax_ = v_logRhol_;

	m_Rjlmax_ = m_Rjl_;
	m_Tikmax_ = m_Tik_;
	Lmax_ = Likelihood_;
}

#endif /* BinaryLBModel_H_ */
