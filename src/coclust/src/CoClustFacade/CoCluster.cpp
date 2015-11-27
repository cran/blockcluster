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


/** @file CoCluster.cpp
 *  @brief Implements CoCluster class.
 **/

#include "CoCluster.h"

bool CoCluster::run()
{
#ifdef COVERBOSE
  std::cout<<"Starting Co-Clustering.."<<std::endl;
#endif

  if (!p_Strategy_||!p_Algo_||!p_Model_||!p_Init_){
#ifdef COVERBOSE
  std::cout<<"Invalid pointer(s).."<<std::endl;
#endif
    return false;
  }

  p_Algo_->SetModel(p_Model_);
  p_Init_->SetModel(p_Model_);
  p_Strategy_->SetInit(p_Init_);
  p_Strategy_->SetAlgo(p_Algo_);
  p_Strategy_->SetModel(p_Model_);
try {
  if(p_Strategy_->run())
  {
    p_Model_->SetErrormsg("Co-Clustering successfully terminated!");
    return true;
  }
} catch (std::exception & e) {
  p_Model_->SetErrormsg(e.what());
#ifdef COVERBOSE
  std::cout<<e.what();
#endif
}
  return false;
}


