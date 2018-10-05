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


#ifndef ICOCLUSTMODEL_H_
#define ICOCLUSTMODEL_H_
/** @file ICoClustModel.h
 *  @brief This file declares the ICoClustModel abstract model class. All the concrete models
 *  classes are derived from this abstract class.
 **/

#include <time.h>
#include <vector>
#include "../typedefs/typedef.h"
#include "../InputParameters/InputParameters.h"

/** @brief This is an abstract class which provides interface for various models.
 *  It provides interfaces for the most common functions in EM, CEM  and SEM algorithms as well
 *  as provide implementation for common functionalities.
 */
class ICoClustModel
{
  private:
    /** Default constructor */
    ICoClustModel();

  public:
    /** Constructor
     *  @param Mparam ModelParameters
     **/
    ICoClustModel(ModelParameters const& Mparam);
    /** Constructor
     * @param Mparam ModelParameters
     * @param rowlabels Row clusters for each row (-1 for unknown cluster for each row)
     * @param collabels Column clusters for each column (-1 unknown cluster for each column)
     * */
    ICoClustModel( ModelParameters const& Mparam, VectorInt const & rowlabels,
                   VectorInt const & collabels);
    /** Destructor*/
    inline virtual ~ICoClustModel(){};

    //various interface functions
    /** Cloning interface*/
    virtual ICoClustModel* clone() = 0;
    /** FUll M-step interface*/
    virtual  void mStepFull() = 0;
    /**Interface for EM Algorithm for rows*/
    virtual bool emRows() = 0;
    /** Interface for calculating log sum for rows*/
    virtual void logSumRows(MatrixReal&) = 0;
    /** Interface for calculating log sum for columns*/
    virtual void logSumCols(MatrixReal&) = 0;
    /**Interface for EM Algorithm for Columns*/
    virtual bool emCols() = 0;
    /**Interface for CEM Algorithm for rows*/
    virtual bool cemRows() = 0;
    /**Interface for CEM Algorithm for Columns*/
    virtual bool cemCols() = 0;
    /**Interface for SEM Algorithm for rows*/
    virtual bool semRows() = 0;
    /**Interface for SEM Algorithm for columns*/
    virtual bool semCols() = 0;
    /**Interface for Gibbs Algorithm for rows*/
    virtual bool GibbsRows() = 0;
    /**Interface for Gibbs Algorithm for columns*/
    virtual bool GibbsCols() = 0;
    /**
     * Interface for calculating Stopping condition using percentage Change in Parameter values. This function will set the
     * ICoClustModel::stopAlgo_ parameter to either true or false depending on whether the change in Likelihood is less than
     * Mparam_.epsilon_ or not respectively.
     */
    virtual void parameterStopCriteria() = 0;
    /** Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
     * ICoClustModel::stopAlgo_ parameter to either true or false depending on whether the change
     * in Likelihood is less than Mparam_.epsilon_ or not respectively.
     */
    virtual void likelihoodStopCriteria();
    /** Interface function for CEM initialization . It will initialize model parameters
     * using CEM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool cemInitStep();
    /** Interface function for Random Initialization. It will initialize model parameters
     * using Random initialization.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool randomInit();
    /** Interface function for EM initialization . It will initialize model parameters
     * using EM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool emInitStep();
    /** Interface function for calculating ICL criteria value.
     * @return ICL criteria.
     */
    virtual STK::Real iclCriteriaValue();
    /** This function will provide the current status of ICoClustModel::stopAlgo_ parameter.*/
    bool stopAlgo();
    /** @return the number of free parmaters of the distribution of a block.*/
    virtual int nbFreeParameters() const = 0;
    /** Interface for finalizing the output. This function will allow the model to finalize
     * all the output parameters after the algorithm terminates.
     */
    virtual void finalizeOutput() = 0;
    /** Interface for output of model parameters on console. */
    virtual void consoleOut() = 0;

    //common computation getter/setter functions
    /** @return ModelParameters */
    ModelParameters const& modelParameters() const;
    /** @param Mparam ModelParameters */
    void setModelParameters(ModelParameters& Mparam);
    /** @param Epsilon Value to be set as Mparam_.epsilon_ */
    void setEpsilon(STK::Real Epsilon);
    /**Get Error message*/
    std::string errorMsg(){return Error_msg_;}
    /**Set Error message*/
    void setMsg(const std::string& s){ Error_msg_ = s;}
    /**Get status of empty cluster*/
    bool isEmptyCluster() const {return empty_cluster_;}
    /**set the value for Empty Cluster*/
    void setEmptyCluster(bool val){empty_cluster_ = val;}
    /** Modify the value of parameters during xem step if better parameters are found.*/
    virtual void modifyThetaStart() = 0;
    /** Copy the value of parameters after xem step to be used during XEM step.*/
    virtual void copyThetaStart() = 0;
    /** Copy the value of parameters after XEM step to be set as final extimated parameters.*/
    virtual void copyThetaMax() = 0;
    /** Modify the value of parameters during XEM step if better parameters are found.*/
    virtual void modifyThetaMax() = 0;
    /** Estimate likelihood value*/
    virtual STK::Real estimateLikelihood() = 0;
    /** E-step for rows*/
    bool eStepRows();
    /** E-step for columns*/
    bool eStepCols();
    /** CE-step for Rows*/
    bool ceStepRows();
    /** CE-step for columns*/
    bool ceStepCols();
    /** SE-step for Rows*/
    bool seStepRows();
    /** SE-step for columns*/
    bool seStepCols();
    /** Draw conditional Row classification vector*/
    void rowClassMatrixDraw();
    /** Draw conditional column classification vector*/
    void colClassMatrixDraw();
    /**Set the known and unknown row labels for semi-supervised coclustering*/
    void setRowLabels(VectorInt const&);
    /**Set the known  and unknown column labels for semi-supervised coclustering*/
    void setColLabels(VectorInt const&);
    /** Get the likelihood value*/
    STK::Real likelihood() const;
    /** This function will return the row classification vector.
     * @return Row classification vector
     */
    VectorInt const& rowClassificationVector() const;
    /** This function will return the column classification vector.
     * @return COlumn classification vector
     */
    VectorInt const& columnClassificationVector() const;
    /** This function will return the row proportions of mixtures models.
     * @return Row proportions of mixing models.
     */
    VectorReal const& rowProportions() const;
    /** This function will return the column proportions of mixtures models.
     * @return column proportions of mixing models.
     */
    const VectorReal & colProportions() const;
    /** This function will return the posterior probabilities for each row.
     * @return Row posterior probabilities (Rows represent original rows and columns represent class probabilities)
     */
    const MatrixReal & rowPosteriorProb() const;
    /** This function will return the posterior probabilities for each column.
     * @return Column posterior probabilities ((Rows represent original columns and columns represent class probabilities))
     */
    const MatrixReal & colPosteriorProb() const;

    template<class T>
    void arrangedDataCluster(T&,const T&);

  protected:
    std::string  Error_msg_;
    /// number of sample (number of rows) and of variables (number of columns)
    int nbSample_,nbVar_;
    //variables use in case of semi-supervised co-clustering
    STK::Array1D<int> UnknownLabelsRows_, UnknownLabelsCols_ ;
    STK::Array1D<std::pair<int,int> > knownLabelsRows_, knownLabelsCols_;
    VectorInt v_nbRowClusterMembers_,v_nbColClusterMembers_;

    ModelParameters Mparam_;
    STK::Real likelihood_,Lmax_;
    bool empty_cluster_;
    MatrixReal m_Tik_, m_Rjl_
             , m_Tikstart_, m_Rjlstart_
             , m_Tikmax_, m_Rjlmax_;
    // sum by column and by row of the posterior probabilities
    VectorReal v_Tk_,v_Rl_;
    // proportions and log proportions
    VectorReal v_Piek_, v_Rhol_, v_logPiek_, v_logRhol_
             , v_logPiekstart_, v_logRholstart_
             , v_logPiekmax_, v_logRholmax_;
    /** Row and column classification matrices respectively*/
    MatrixInt m_Zik_,m_Wjl_;
    /**Row and column classification vector*/
    VectorInt v_Zi_, v_Wj_;
    /** Boolean Variable used to identify whether the stopping condition
     *  is reached or not*/
    bool stopAlgo_;
    STK::Real dimprod_;

    // Utility Functions
    /** Generate k Random integers in the range 0 to n-1
     *  @param n Max integer value
     *  @param k Number of random integers to be generated (should not exceed n)
     */
    VectorInt randSample(int n,int k);
    /** compute the vector v_Tk_ and check if the size block is not too small
     *  @return true if the size block is under the threshold, false otherwise
     **/
    bool finalizeStepRows();
    /** compute the vector v_Rl_ and check if the size block is not too small
     *  @return true if the size block is under the threshold, false otherwise
     **/
    bool finalizeStepCols();
    /** Parameter finalization common for all models*/
    void commonFinalizeOutput();
    /** generate random m_Tik_
     *  @return the minimal number of individuals in class */
    int randomFuzzyTik();
    /** generate random m_Rjl_
     *  @return the minimal number of individuals in class */
    int randomFuzzyRjl();

    VectorInt partRnd(int n,VectorReal proba);
    VectorReal cumSum(VectorReal proba);
    PointReal unifRnd(STK::Real a, STK::Real b, int col);

  private:
    //make assignment operator private
    ICoClustModel& operator=(const ICoClustModel&);
};

inline ModelParameters const& ICoClustModel::modelParameters() const
{ return Mparam_;}

inline void ICoClustModel::setModelParameters(ModelParameters& Mparam)
{ Mparam_ = Mparam;}

inline void ICoClustModel::setEpsilon(STK::Real Epsilon)
{ Mparam_.epsilon_ = Epsilon;}

inline bool ICoClustModel::stopAlgo()
{ return stopAlgo_;}

inline  VectorInt const& ICoClustModel::rowClassificationVector() const
{ return v_Zi_;}

inline  VectorInt const& ICoClustModel::columnClassificationVector() const
{ return v_Wj_;}

inline VectorReal const& ICoClustModel::rowProportions() const
{ return v_Piek_;}

inline VectorReal const& ICoClustModel::colProportions() const
{ return v_Rhol_;}

inline MatrixReal const& ICoClustModel::rowPosteriorProb() const
{return m_Tik_;}

inline MatrixReal const& ICoClustModel::colPosteriorProb() const
{return m_Rjl_;}

inline STK::Real ICoClustModel::likelihood() const
{ return likelihood_;}

template<class T>
void ICoClustModel::arrangedDataCluster(T& m_ClusterDataij_,const T& m_Dataij_)
{
  STK::VectorXi v_Zi = v_Zi_;
  STK::VectorXi v_Wj = v_Wj_;
  m_ClusterDataij_.resize(nbSample_,nbVar_);
  m_ClusterDataij_ = 0;
  //Rearrange data into clusters

  VectorInt rowincrement(Mparam_.nbrowclust_, 0);
  VectorInt nbindrows(Mparam_.nbrowclust_+1, 0);

  for ( int k = 1; k < Mparam_.nbrowclust_; ++k)
  { nbindrows[k] = (v_Zi==(k-1)).count()+nbindrows[k-1];}

  VectorInt colincrement(Mparam_.nbcolclust_, 0);
  VectorInt nbindcols(Mparam_.nbcolclust_+1, 0);

  for ( int l = 1; l < Mparam_.nbcolclust_; ++l)
  { nbindcols[l]= (v_Wj==(l-1)).count()+nbindcols[l-1];}

  for ( int j = 0; j < nbVar_; ++j)
  {
    m_ClusterDataij_.col(colincrement[v_Wj[j]] + nbindcols[v_Wj[j]]) = m_Dataij_.col(j);
    colincrement[v_Wj[j]]+=1;
  }
  T temp = m_ClusterDataij_;

  for ( int i = 0; i < nbSample_; ++i)
  {
    m_ClusterDataij_.row( rowincrement[v_Zi[i]] + nbindrows[v_Zi[i]]) = temp.row(i);
    rowincrement[v_Zi[i]]+=1;
  }

}

#endif /* ICOCLUSTMODEL_H_ */
