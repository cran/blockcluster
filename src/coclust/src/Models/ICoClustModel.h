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

/** @brief ICoClustModel is an abstract class which provides interface for various models.
 *  It provides interfaces for the most common functions in EM, CEM  and SEM algorithms as well
 *  as provide implementation for common functionalities.
 */
class ICoClustModel
{
    public:
    /**Constructor
     * @param Mparam ModelParameters
     * */
    ICoClustModel(ModelParameters const& Mparam);
    /** Constructor
     * @param Mparam ModelParameters
     * @param rowlabels Row clusters for each row (-1 for unknown cluster for each row)
     * @param collabels Column clusters for each column (-1 unknown cluster for each column)
     * */
    ICoClustModel(ModelParameters const& Mparam,VectorInteger const & rowlabels,
                  VectorInteger const & collabels);
    /** Destructor*/
    virtual ~ICoClustModel(){};

    //various interface functions
    /** Cloning interface*/
    virtual ICoClustModel* Clone() = 0;
    /** FUll M-step interface*/
    virtual  void MStepFull() = 0;
    /**Interface for EM Algorithm for rows*/
    virtual bool EMRows() = 0;
    /** Interface for calculating log sum for rows*/
    virtual void LogSumRows(MatrixReal&) = 0;
    /** Interface for calculating log sum for columns*/
    virtual void LogSumCols(MatrixReal&) = 0;
    /**Interface for EM Algorithm for Columns*/
    virtual bool EMCols() = 0;
    /**Interface for CEM Algorithm for rows*/
    virtual bool CEMRows() = 0;
    /**Interface for CEM Algorithm for Columns*/
    virtual bool CEMCols() = 0;
    /**Interface for SEM Algorithm for rows*/
    virtual bool SEMRows() = 0;
    /**Interface for SEM Algorithm for columns*/
    virtual bool SEMCols() = 0;
    /**
     * Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
     * ICoClustModel::StopAlgo parameter to either true or false depending on whether the change
     * in Likelihood is less than Mparam_.epsilon_ or not respectively.
     */
    virtual void likelihoodStopCriteria() = 0;
    /**
     * Interface for calculating Stopping condition using percentage Change in Parameter values. This function will set the
     * ICoClustModel::StopAlgo parameter to either true or false depending on whether the change in Likelihood is less than
     * Mparam_.epsilon_ or not respectively.
     */
    virtual void ParameterStopCriteria() = 0;
    /**
     * Interface function for CEM initialization . It will initialize model parameters
     * using CEM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool CEMInit();
    /**
     * Interface function for Random Initialization. It will initialize model parameters
     * using Random initialization.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool RandomInit();
    /**
     * Interface function for CEM initialization . It will initialize model parameters
     * using Fuzzy CEM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool FuzzyCEMInit();
    /** This function will provide the current status of ICoClustModel::StopAlgo parameter.*/
    bool stopAlgo();
    /**
     * Interface for finalizing the output. This function will allow the model to finalize
     * all the output parameters after the algorithm terminates.
     */
    virtual void FinalizeOutput() = 0;
    /**
     * Interface for output of model parameters on console.
     */
    virtual void ConsoleOut() = 0;
    /**
     * @return ModelParameters
     */
    ModelParameters const& GetModelParameters() const;
    /**
     * @param Mparam ModelParameters
     */
    void SetModelParameters(ModelParameters& Mparam);
    /**
     *
     * @param Epsilon Value to be set as Mparam_.epsilon_
     */
    void SetEpsilon(float Epsilon);
    /**Get Error message*/
    std::string GetErrormsg(){return Error_msg_;}
    /**Set Error message*/
    void SetErrormsg(const std::string& s){Error_msg_ = s;}
    /**Get status of empty cluster*/
    bool isEmptyCluster() const {return empty_cluster;}
    /**set the value for Empty Cluster*/
    void SetEmptyCluster(bool val){empty_cluster = val;}
    /** Modify the value of parameters during xem step if better parameters are found.*/
    virtual void Modify_theta_start() = 0;
    /** Copy the value of parameters after xem step to be used during XEM step.*/
    virtual void Copy_theta_start() = 0;
    /** Copy the value of parameters after XEM step to be set as final extimated parameters.*/
    virtual void Copy_theta_max() = 0;
    /** Modify the value of parameters during XEM step if better parameters are found.*/
    virtual void Modify_theta_max() = 0;
    /** Estimate likelihood value*/
    virtual float EstimateLikelihood() = 0;


    //common computation functions
    /** E-step for rows*/
    bool ERows();
    /** E-step for columns*/
    bool ECols();
    /** CE-step for Rows*/
    bool CERows();
    /** CE-step for columns*/
    bool CECols();
    /** SE-step for Rows*/
    bool SERows();
    /** SE-step for columns*/
    bool SECols();
    /** Draw conditional Row classification vector*/
    void RowClassMatrixdraw();
    /** Draw conditional column classification vector*/
    void ColClassMatrixdraw();


    //common setter functions
    /**Set the known and unknown row labels for semi-supervised coclustering*/
    void SetRowLabels(const VectorInteger&);
    /**Set the known  and unknown column labels for semi-supervised coclustering*/
    void SetColLabels(const VectorInteger&);

    //common getter functions
    /** Get the likelihood value*/
    float GetLikelihood() const;
    /**
     * This function will return the row classification vector.
     * @return Row classification vector
     */
    const VectorInteger& GetRowClassificationVector() const;
    /**
     * This function will return the column classification vector.
     * @return COlumn classification vector
     */
    const VectorInteger& GetColumnClassificationVector() const;
    /**
     * This function will return the row proportions of mixtures models.
     * @return Row proportions of mixing models.
     */
    const VectorReal& GetRowProportions() const;
    /**
     * This function will return the column proportions of mixtures models.
     * @return column proportions of mixing models.
     */
    const VectorReal & GetColProportions() const;
    /**
     * This function will return the posterior probabilities for each row.
     * @return Row posterior probabilities (Rows represent original rows and columns represent class probabilities)
     */
    const MatrixReal & GetRowPosteriorprob() const;
    /**
     * This function will return the posterior probabilities for each column.
     * @return Column posterior probabilities ((Rows represent original columns and columns represent class probabilities))
     */
    const MatrixReal & GetColPosteriorprob() const;

  protected:
    std::string  Error_msg_;
    int nbSample_,nbVar_;

    //variables use in case of semi-supervised co-clustering
    std::vector<int> UnknownLabelsRows_, UnknownLabelsCols_ ;
    std::vector<std::pair<int,int> > knownLabelsRows_, knownLabelsCols_;
    VectorInteger v_nbRowClusterMembers_,v_nbColClusterMembers_;

    bool empty_cluster;
    ModelParameters Mparam_;
    float Likelihood_,Lmax_;
    VectorInteger v_Zi_;/**<Row classification vector*/
    VectorInteger v_Wj_;/**<Column classification vector*/
    MatrixReal m_Tik_,m_Rjl_,m_Rjlstart_,m_Tikstart_,m_Rjlmax_,m_Tikmax_;
    VectorReal v_Tk_,v_Rl_;
    VectorReal v_Piek_,v_logPiek_,v_logPiekstart_,v_logPiekmax_,v_logRhol_,v_Rhol_,v_logRholstart_,v_logRholmax_;
    MatrixInteger m_Zik_,m_Wjl_;/**Row and column classification matrices respectively*/
    bool StopAlgo;/**<Boolean Variable used to identify whether the stopping condition is reached or not*/

    // Utility Functions
    /**
     * Generate k Random integers in the range 0 to n-1
     * @param n Max integer value
     * @param k Number of random integers to be generated (should not exceed n)
     */
    VectorInteger RandSample(int n,int k);
    /** Parameter finalization common for all models*/
    void CommonFinalizeOutput();

  private:
    //make assignment operator private
    ICoClustModel& operator=(const ICoClustModel&);

};

inline ModelParameters const& ICoClustModel::GetModelParameters() const
{
  return Mparam_;
}

inline void ICoClustModel::SetModelParameters(ModelParameters& Mparam)
{
  Mparam_ = Mparam;
}

inline void ICoClustModel::SetEpsilon(float Epsilon)
{
  Mparam_.epsilon_ = Epsilon;
}

inline bool ICoClustModel::stopAlgo()
{
  return StopAlgo;
}

inline const VectorInteger & ICoClustModel::GetRowClassificationVector() const
{
  return v_Zi_;
}

inline const VectorInteger & ICoClustModel::GetColumnClassificationVector() const
{
  return v_Wj_;
}

inline const VectorReal& ICoClustModel::GetRowProportions() const
{
  return v_Piek_;
}

inline const VectorReal& ICoClustModel::GetColProportions() const
{
  return v_Rhol_;
}

inline const MatrixReal& ICoClustModel::GetRowPosteriorprob() const
{
	return m_Tik_;
}

inline const MatrixReal& ICoClustModel::GetColPosteriorprob() const
{
	return m_Rjl_;
}

inline float ICoClustModel::GetLikelihood() const
{
	return Likelihood_;
}

inline void ICoClustModel::SetRowLabels(const VectorInteger & rowlabels)
{
  int cluster, length = rowlabels.size();
  for ( int i = 0; i < length; ++i) {
    cluster = rowlabels(i);
    if (cluster<0) {
      UnknownLabelsRows_.push_back(i);
    } else {
      knownLabelsRows_.push_back(std::pair<int,int>(i,cluster));
      v_nbRowClusterMembers_(cluster)+=v_nbRowClusterMembers_(cluster);
      v_Zi_(i) = 1;
      m_Zik_(i,cluster) = 1;
    }
  }
}

inline void ICoClustModel::SetColLabels(const VectorInteger & collabels)
{
  int cluster , length = collabels.size();
  for ( int j = 0; j < length; ++j) {
    cluster = collabels(j);
    if (cluster<0) {
      UnknownLabelsCols_.push_back(j);
    } else {
      knownLabelsCols_.push_back(std::pair<int,int>(j,cluster));
      v_nbColClusterMembers_(cluster)+=v_nbColClusterMembers_(cluster);
      v_Wj_(j) = 1;
      m_Wjl_(j,cluster) = 1;
    }
  }
}

inline void ICoClustModel::CommonFinalizeOutput()
{
  // Calculate row and column proportions
  if(!Mparam_.fixedproportions_){
    v_Piek_ = v_logPiek_.array().exp();
    v_Rhol_ = v_logRhol_.array().exp();
  }
  else
  {
    v_Piek_ = (1.0/Mparam_.nbrowclust_)*MatrixReal::Ones(Mparam_.nbrowclust_,1);
    v_Rhol_ = (1.0/Mparam_.nbcolclust_)*MatrixReal::Ones(Mparam_.nbcolclust_,1);
  }

  VectorReal::Index maxIndex;
  // Calculate classification vector for rows
  for ( int i = 0; i < nbSample_; ++i) {
  m_Tik_.row(i).maxCoeff(&maxIndex);
  v_Zi_(i) = maxIndex;
  }

  // Calculate classification vector for columns
  for ( int j = 0; j < nbVar_; ++j) {
  m_Rjl_.row(j).maxCoeff(&maxIndex);
  v_Wj_(j) = maxIndex;
  }

}

#endif /* ICOCLUSTMODEL_H_ */
