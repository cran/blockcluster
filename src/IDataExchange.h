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
 * Project:  RCocluster
 * created on: Feb 22, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ICoClust.h
 *  @brief 
 **/

#ifndef IDataExchange_H_
#define IDataExchange_H_
#include <iostream>
#include "conversion.h"
#include "coclust/src/typedefs/typedef.h"
#include "coclust/src/Models/ICoClustModel.h"
#include "coclust/src/InputParameters/InputParameters.h"
#include <Rcpp.h>

/** @brief
 *
 */
class IDataExchange
{
  public:
    IDataExchange(){InitializeParamEnum();}
    void InitializeParamEnum();
    virtual void DataInput(Rcpp::S4 & obj) = 0;
    virtual void Output(Rcpp::S4& obj,ICoClustModel*, bool) = 0;
    void SetInput(Rcpp::S4 & obj);
    Strategy& GetStrategy(){return strategy_;}
    StrategyParameters& GetStrategyParameters(){return Stratparam_;}
    ModelParameters& GetModelParameters(){return Mparam_;}
    const VectorInteger& GetRowLabels() const {return v_rowlabels_;}
    const VectorInteger& GetColLabels() const {return v_collabels_;}
    virtual ~IDataExchange();
  protected:
    Strategy strategy_;
    StrategyParameters Stratparam_;
    ModelParameters Mparam_;
    VectorInteger v_rowlabels_,v_collabels_;
    std::map<std::string,Algorithm> S_Algorithm;
    std::map<std::string,StopCriteria> S_StopCriteria;
    std::map<std::string,DataType> S_DataType;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model> S_Model;
};
#endif /*IDataExchange_H_*/
