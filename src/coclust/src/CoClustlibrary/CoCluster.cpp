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
 * created on: Nov 28, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file CoCluster.cpp
 *  @brief Implements CoCluster class.
 **/

#include "CoCluster.h"
#include "../../../stkpp/include/STKernel.h"

bool CoCluster::run()
{
#ifdef COVERBOSE
  std::cout<<"Starting Co-Clustering.."<<std::endl;
#endif
  // Assigning model to stopping criteria
  if (p_StopCriteria_) p_StopCriteria_->SetModel(p_Model_);
  else
    throw STK::runtime_error(_T("In CoCluster::run(), StopCriteria is not instantiated\n"));

  //Assigning model to Initialization method
  if(p_Init_) p_Init_->SetModel(p_Model_);
  else
    throw STK::runtime_error(_T("In CoCluster::run(), initialization is not instantiated\n"));
  // setting model to Algo
  if (p_Algo_) {
    p_Algo_->SetModel(p_Model_);
    p_Algo_->SetStopCriteria(p_StopCriteria_);
    p_Algo_->SetInitialization(p_Init_);
  }
  else
    throw STK::runtime_error(_T("In CoCluster::run(), Algorithm is not instantiated\n"));
  if(p_Algo_->run())
    return true;
  else
    return false;
}


