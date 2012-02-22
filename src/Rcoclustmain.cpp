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
 * created on: Jan 10, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file Rcoclustermain.cpp
 *  @brief 
 **/

#include <Rcpp.h>
#include <time.h>
#include <map>
#include "ICoClust.h"
#include "ContinuousCoClust.h"
#include "BinaryCoClust.h"
#include "ContingencyCoClust.h"
#include "coclust/src/enumerations/enumerations.h"
using namespace std;
RcppExport SEXP CoClustmain(SEXP robj)
{
  BEGIN_RCPP
  Rcpp::S4 CoClustobj(robj);
  Rcpp::S4 strategy(CoClustobj.slot("strategy"));
  //Initialize enums
  std::map<std::string,DataType> S_DataType;
  //Datatype
  S_DataType["binary"] = Binary;
  S_DataType["contingency"] = Contingency;
  S_DataType["continuous"] = Continuous;

  DataType datatype = S_DataType[Rcpp::as<std::string>(CoClustobj.slot("datatype"))];
  ICoClust * p_Coclust_;

  switch (datatype) {
    case Binary:
      p_Coclust_ = new BinaryCoClust();
      break;
    case Contingency:
      p_Coclust_ = new ContingencyCoClust();
      break;
    case Continuous:
      p_Coclust_ = new ContinuousCoClust();
      break;
    default:
      break;
  }
  p_Coclust_->InitializeEnum();
  p_Coclust_->SetStrategy(CoClustobj,strategy);
  p_Coclust_->ModelInput(CoClustobj,strategy);
  p_Coclust_->SetModel();
  double begin =clock();
  bool success = p_Coclust_->Run();
  p_Coclust_->Output(CoClustobj,success);
  CoClustobj.slot("time") = double(clock()-begin)/CLOCKS_PER_SEC;

  delete p_Coclust_;
  return CoClustobj;
  END_RCPP
}
