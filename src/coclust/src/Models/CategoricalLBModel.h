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

#ifndef CATEGORICALLBMODEL_H_
#define CATEGORICALLBMODEL_H_
/**@file CategoricalLBModel.h
 * @brief 
 */
#include "ICoClustModel.h"

class CategoricalLBModel:public ICoClustModel
{
  public:
    CategoricalLBModel(MatrixInteger const& m_Dataij,ModelParameters const& Mparam,int a=1,int b=1);
    CategoricalLBModel(MatrixInteger const& m_Dataij,VectorInteger const & rowlabels,
                       VectorInteger const & collabels,ModelParameters const& Mparam,int a=1,int b=1);
    virtual CategoricalLBModel* Clone(){return new CategoricalLBModel(*this);}
    virtual void LogSumRows(MatrixReal & m_sum);
    virtual void LogSumCols(MatrixReal & m_sum);
    virtual void MStepFull(){};
    virtual bool EMRows();
    virtual bool CEMRows();
    virtual bool EMCols();
    virtual bool CEMCols();
    virtual bool SEMRows();
    virtual bool SEMCols();
    virtual void likelihoodStopCriteria();
    virtual void ParameterStopCriteria();
    virtual bool RandomInit();
    virtual float ICLCriteriaValue();
    virtual void FinalizeOutput();
    virtual void ConsoleOut();
    virtual void Modify_theta_start();
    virtual void Copy_theta_start();
    virtual void Copy_theta_max();
    virtual void Modify_theta_max();
    virtual float EstimateLikelihood();
    const MatrixInteger& GetArrangedDataClusters();
    inline const std::vector<MatrixReal>& Getmean(){return m3_Alphahkl_;}
    virtual ~CategoricalLBModel();

  protected:
    int a_,b_;//hyper-parameters
    const MatrixInteger& m_Dataij_;
    MatrixInteger m_ClusterDataij_;
    int r_; //number of categories
    std::vector<MatrixReal> m3_Alphahkl_,m3_Alphahklold_,m3_Alphahkl1_,
    m3_Alphahkl1old_,m3_Alphahklstart_,m3_Alphahklmax_,m3_logAlhphahkl_;
    std::vector<MatrixBinary> m3_Yhij_,m3_Yijh_,m3_Yjih_;//different ways to store data
    float Likelihood_old;

    virtual void MStepRows();
    virtual void MStepCols();


};



#endif /* CATEGORICALLBMODEL_H_ */
