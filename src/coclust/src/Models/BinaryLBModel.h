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
#include "../../CImg/CImg.h"
#endif

class BinaryLBModel : public ICoClustModel
{
  public:
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b bayesian hyperparameters
     * */
    BinaryLBModel( MatrixBinary const&  m_Dataij
                 , ModelParameters const& Mparam
                 , int a=1,int b=1);
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param rowlabels various labels for rows (-1  for unknown row label)
     * @param collabels various labels for columns (-1 for unknown column label)
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b Bayesian hyperparameters
     * */
    BinaryLBModel( MatrixBinary const&  m_Dataij
                 , VectorInteger const& rowlabels
                 , VectorInteger const& collabels
                 , ModelParameters const& Mparam
                 , int a=1,int b=1);

    /** cloning */
    virtual BinaryLBModel* clone(){return new BinaryLBModel(*this);}
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
    virtual bool cemInitStep();
    virtual bool emInitStep();
    virtual void finalizeOutput();
    virtual STK::Real iclCriteriaValue();
    virtual void consoleOut();
    virtual void modifyThetaStart();
    virtual void copyThetaStart();
    virtual void copyThetaMax();
    virtual void modifyThetaMax();
    virtual STK::Real estimateLikelihood();
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;

    MatrixBinary const&  arrangedDataClusters();
    /**Return class mean BinaryLBModel::m_akl_ for all the blocks (co-clusters)*/
    MatrixBinary const&  mean() const;
    /**Return Class despersion BinaryLBModel::m_epsilonkl_ for all the blocks (co-clusters) */
    MatrixReal const& dispersion() const;
    /**Destructor*/
    inline virtual ~BinaryLBModel(){};

#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void displayCluster();
#endif

  protected:
    //Variables involved in Bernouilli model
    int a_,b_;//hyper-parameters
    MatrixBinary const&  m_Dataij_;
    MatrixBinary m_ClusterDataij_;
    MatrixReal m_Vjk_;
    MatrixReal m_Uil_;
    MatrixReal m_Alphakl_, m_Alphaklold_, m_Alphakl1_, m_Alphakl1old_,m_Alphaklmax_,m_Alphaklstart_;
    MatrixBinary m_akl_;
    MatrixReal m_epsilonkl_,m_epsilonklstart_,m_epsilonklmax_;

    void mStepRows();
    void mStepCols();

    void mGibbsStepRows();
    void mGibbsStepCols();

    // Functions used to operate on data in intermediate steps when running the Initialization
    bool initCEMRows();
    bool initCEMCols();
    bool initEMRows();
    bool initEMCols();
    void initBernoulliLogSumRows(MatrixReal & m_sum);
    void initBernoulliLogSumCols(MatrixReal & m_sum);
    void selectRandomColsFromData(MatrixReal& m,int col);
    void generateRandomBernoulliParameterRows(MatrixReal& m, int cols);
    void generateRandomBernoulliParameterCols(MatrixReal& m);
};


inline MatrixBinary const&  BinaryLBModel::mean() const
{ return m_akl_;}

inline MatrixReal const& BinaryLBModel::dispersion() const
{ return m_epsilonkl_;}

#endif /* BinaryLBModel_H_ */
