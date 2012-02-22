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

/** @file ICoClustModel.h
 *  @brief This file declares the ICoClustModel abstract model class.
 **/

#ifndef ICOCLUSTMODEL_H_
#define ICOCLUSTMODEL_H_

#include <time.h>
#include "../typedefs/typedef.h"
#include "../InputParameters/InputParameters.h"
#include "../../../stkpp/include/Model.h"
//only need to declare initialization abstract class
class IInit;
/** @brief ICoClustModel is an abstract class which provides interface for various models.
 *  It provides interface for the most common
 *  functions in EM and CEM algorithms.
 */
class ICoClustModel : public STK::IModel
{
  protected:
    /** constructor*/
    ICoClustModel(ModelParameters const& Mparam)
    {
      Mparam_ = Mparam;
    }

  public:
    /** Destructor*/
    virtual ~ICoClustModel(){};

    /**Interface for Expectation step.*/
    virtual bool Estep() = 0;
    /**
     * Interface for Expectation and classification step. The model implementing CEM algorithm, must implement this
     * function. The Expectation and classification step is performed in the same function.
     */
    virtual bool CEstep() = 0;
    /**Interface for Maximization step.*/
    virtual void Mstep() = 0;
    /**
     * Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
     * ICoClustModel::StopAlgo parameter to either true or false depending on whether the change
     * in Likelihood is less than Mparam_.epsilon_ or not respectively.
     */
    virtual void likelihoodStopCriteria() = 0;
    /**
     * Interface for calculating Stopping condition using percentage Change in Parameter. This function will set the
     * ICoClustModel::StopAlgo parameter to either true or false depending on whether the change in Likelihood is less than
     * Mparam_.epsilon_ or not respectively.
     */
    virtual void ParameterStopCriteria() = 0;
    /**
     * Interface function for CEM initialization . It should initialize model parameters
     * using CEM algorithm.
     */
    virtual bool CEMInit() = 0;
    /**
     * Interface function for Random Initialization.It should initialize model parameters
     * using Random initialization.
     */
    virtual bool RandomInit() = 0;
    /**
     * Interface function for CEM initialization . It should initialize model parameters
     * using Fuzzy CEM algorithm.
     */
    virtual bool FuzzyCEMInit() = 0;
    /** This function will provide the current status of ICoClustModel::StopAlgo parameter.*/
    bool stopAlgo();
    /**
     * Interface for finalizing the output. This function will allow the model to finalize all the output parameters
     * after the stopping condition is reached.
     */
    virtual void FinalizeOutput() = 0;
    /**
     * Interface for output of model parameters on console.
     */
    virtual void ConsoleOut() = 0;
    /**
     * This function will return the row classification vector.
     * @return Row classification vector
     */
    VectorInteger& GetRowClassificationVector();
    /**
     * This function will return the column classification vector.
     * @return COlumn classification vector
     */
    VectorInteger& GetColumnClassificationVector();

    /**
     * This function will return the row proportions of mixtures models.
     * @return Row proportions of mixing models.
     */
    VectorReal & GetRowProportions();
    /**
     * This function will return the column proportions of mixtures models.
     * @return column proportions of mixing models.
     */
    VectorReal & GetColProportions();
    /**
     * This function will return the posterior probabilities for each row.
     * @return Row posterior probabilities (Rows represent original rows and columns represent class probabilities)
     */
    MatrixReal & GetRowPosteriorprob();
    /**
     * This function will return the posterior probabilities for each column.
     * @return Column posterior probabilities ((Rows represent original columns and columns represent class probabilities))
     */
    MatrixReal & GetColPosteriorprob();
    /**Set Epsilon*/
    void SetEpsilon(int phase);
    /**Get Rmessage*/
    std::string Rmessage(){return R_errormsg;}
    /**Set Rmessage*/
    void SetRmessage(std::string s){R_errormsg = s;}
    /**Get status of empty cluster*/
    bool isEmptyCluster(){return empty_cluster;}
    /**set the value for Empty Cluster*/
    void SetEmptyCluster(bool val){empty_cluster = val;}
    virtual void Modify_theta_start() = 0;
    virtual void Copy_theta_start() = 0;
    virtual void Copy_theta_max() = 0;
    virtual void Modify_theta_max() = 0;
    float GetLikelihood();
    virtual float EstimateLikelihood() = 0;
    /*
     * This will set reference for IInit::p_Init_ to the instantiated Initialization inside CoClustermain.cpp
     */
    void SetInitAlgo(IInit*& init);


  protected:
    std::string  R_errormsg;
    bool empty_cluster;
    ModelParameters Mparam_;
    float Likelihood_,Lmax_;
    VectorInteger v_Zi_;/**<Row classification vector*/
    VectorInteger v_Wj_;/**<Column classification vector*/
    MatrixReal m_Tik_,m_Rjl_,m_Rjlstart_,m_Rjlmax_,m_Tikmax_;
    VectorReal v_Tk_,v_Rl_;
    VectorReal v_Piek_,v_logPiek_,v_logPiekstart_,v_logPiekmax_,v_logRhol_,v_Rhol_,v_logRholstart_,v_logRholmax_;
    MatrixInteger m_Zik_,m_Wjl_;
    bool StopAlgo;/**<Boolean Variable used to identify whether the stopping condition is reached or not*/
    //IInit * p_Init_;

    // Utility Functions

    /**
     * Generate Random intergers
     * @param n Max integer value
     * @param k Number of random integers to be generated (should not exceed n)
     */
    VectorInteger RandSample(int n,int k);

  private:
    /** a CoClustModel should never be able to run
     * (it' s a kitchen without chef)
     * @return false in all case
     **/
    inline bool run() { return false;}
};

inline void ICoClustModel::SetEpsilon(int phase)
{
  if(phase== 0 )
    Mparam_.epsilon_ = Mparam_.eps_xem_;
  else {
    Mparam_.epsilon_ = Mparam_.eps_XEM_;
  }
}
inline bool ICoClustModel::stopAlgo()
{
  return StopAlgo;
}

inline VectorInteger & ICoClustModel::GetRowClassificationVector()
{
  // Calculate classification vector for rows
  VectorReal::Index maxIndex;
  for (int i = 0; i < nbSample_; ++i) {
	m_Tik_.row(i).maxCoeff(&maxIndex);
	v_Zi_(i) = maxIndex;
  }
  return v_Zi_;
}

inline VectorInteger & ICoClustModel::GetColumnClassificationVector()
{
  // Calculate classification vector for columns
  VectorReal::Index maxIndex;
  for (int j = 0; j < nbVar_; ++j) {
	m_Rjl_.row(j).maxCoeff(&maxIndex);
	v_Wj_(j) = maxIndex;
  }
  return v_Wj_;
}

inline VectorReal& ICoClustModel::GetRowProportions()
{
  return v_Piek_;
}

inline VectorReal& ICoClustModel::GetColProportions()
{
  return v_Rhol_;
}

inline MatrixReal& ICoClustModel::GetRowPosteriorprob()
{
	return m_Tik_;
}

inline MatrixReal& ICoClustModel::GetColPosteriorprob()
{
	return m_Rjl_;
}

/*inline void ICoClustModel::SetInitAlgo(IInit*& init)
{
  p_Init_ = init;
}*/

inline float ICoClustModel::GetLikelihood()
{
	return Likelihood_;
}

#endif /* ICOCLUSTMODEL_H_ */
