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
#include "../../../CImg/CImg.h"
#endif
class BinaryLBModelequalepsilon : public ICoClustModel
{
  public:
    /**Constructor
     * @param m_Dataij a constant reference on the data set.
     * */
    BinaryLBModelequalepsilon(MatrixBinary const& m_Dataij,ModelParameters const& Mparam);

    virtual bool Estep();
    virtual bool CEstep();
    virtual void Mstep();
    virtual void likelihoodStopCriteria();
    virtual float EstimateLikelihood();
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
    /**Return class mean BinaryLBModelequalepsilon::m_akl_ for all the blocks (co-clusters)*/
    MatrixBinary Getmean();
    /**Return Class despersion BinaryLBModelequalepsilon::m_epsilonkl_ for all the blocks (co-clusters) */
    float Getdispersion();
    /**Destructor*/
    inline virtual ~BinaryLBModelequalepsilon(){};

    MatrixBinary GetArrangedDataClusters();
#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void DisplayCluster();
#endif

  protected:
    //Variables involved in Bernoulli model
    MatrixBinary const& m_Dataij_;
    MatrixReal m_Xjl_,m_Xik_,m_Tk_Rl_;
    float dimprod_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_,m_Ukl_;
    VectorReal v_Ui_;
    MatrixReal m_Ykl_,m_Ykl_old2_,m_Ykl_old1_;
    MatrixBinary m_Akl_,m_Aklstart_,m_Aklmax_;
    float Epsilon_,Epsilonstart_,Epsilonmax_;
    float Likelihood_old;
    float W1_,W1_old_;

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
    void SelectRandomColsfromdata(MatrixReal& m,int col);
    void SelectRandomRows(MatrixReal& m_lk);
    void GenerateRandomMean(MatrixBinary& m);
};


inline MatrixBinary BinaryLBModelequalepsilon::Getmean()
{
  return m_Akl_;
}

inline float BinaryLBModelequalepsilon::Getdispersion()
{
  return Epsilon_;
}

inline void BinaryLBModelequalepsilon::Modify_theta_start()
{
	m_Aklstart_ = m_Akl_;
	Epsilonstart_ = Epsilon_;
  v_logPiekstart_ = v_logPiek_;
  v_logRholstart_ = v_logRhol_;

  m_Rjlstart_ = m_Rjl_;
}

inline void BinaryLBModelequalepsilon::Copy_theta_start()
{
	m_Akl_ = m_Aklstart_;
	Epsilon_ = Epsilonstart_;
  v_logPiek_ = v_logPiekstart_;
  v_logRhol_ = v_logRholstart_;

  m_Rjl_ = m_Rjlstart_;

  //initialization
  v_Rl_ = m_Rjl_.colwise().sum();
  m_Ykl_old1_ = m_Ykl_;
}

inline void BinaryLBModelequalepsilon::Copy_theta_max()
{
	m_Akl_ = m_Aklmax_;
	Epsilon_ = Epsilonmax_;
  v_logPiek_ = v_logPiekmax_;
  v_logRhol_ = v_logRholmax_;

  m_Tik_ = m_Tikmax_;
  m_Rjl_ = m_Rjlmax_;
  Likelihood_ = Lmax_;
}

inline void BinaryLBModelequalepsilon::Modify_theta_max()
{
	m_Aklmax_ = m_Akl_;
	Epsilonmax_ = Epsilon_;
  v_logPiekmax_ = v_logPiek_;
  v_logRholmax_ = v_logRhol_;

  m_Rjlmax_ = m_Rjl_;
  m_Tikmax_ = m_Tik_;
  Lmax_ = Likelihood_;
}

#endif /* BINARYLBMODELEQUALEPSILON_H_ */
