/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Author : Parmeet Bhatia
 Contact : bhatia.parmeet@gmail.com , serge.iovleff@stkpp.org

 */


/** @file XStrategyAlgo.cpp
 *  @brief Implements XStrategyAlgo class.
 **/

#include "XStrategyAlgo.h"
#include <iostream>
#include <fstream>

void XStrategyAlgo::writeVector(int it, MatrixReal classes, std::string type)
{
  std::ostringstream stream;
  stream << "data/" << it << "-" << type << ".txt";
  std::string str = stream.str();
  std::cout << str << std::endl;
  std::ofstream o_file(str.c_str());
  std::cout << classes.sizeCols() << std::endl;
  std::cout << classes.sizeRows() << std::endl;
  for (int i = 0; i < classes.sizeRows(); ++i)
  {
    int maxInd = 0;
    int maxVal = classes(i, 0);
    for (int j = 0; j < classes.sizeCols(); ++j)
    {
      if (classes(i, j) > maxVal) maxInd = j;
    }
    o_file << maxInd << std::endl;
  }
}

bool XStrategyAlgo::run()
{
#ifdef COVERBOSE
  std::cout<<"Entering XStrategy::run()"<<"\n";
#endif
  STK::Real Lmax = -RealMax, L1 = -RealMax, Lcurrent;
  int ntry_empty=0;
  bool non_empty = false;
  bool isinitialized = false;
  for ( int itry = 0; itry < Stratparam_.nbtry_; ++itry)
  {
    non_empty = false;
    L1 = -RealMax;
    //set espilon to eps_xem
    p_Model_->setEpsilon(p_Model_->modelParameters().eps_xem_);
    for ( int ixem = 0; ixem < Stratparam_.nbxem_; ++ixem)
    {
      ntry_empty = 0;
      p_Model_->setEmptyCluster(true);
      while(ntry_empty<100)
      {
        if (p_Init_->run())
        {
          isinitialized = true;
          for ( int itr = 0; itr < Stratparam_.nbiter_xem_; ++itr)
          {
            if(p_Algo_->run()) // run initial algorithm was successful
            {
              (p_Model_->*Stratparam_.Stop_Criteria)();
              if (p_Model_->stopAlgo()) { break;} // algo convergence
            }
            else { break;} // algo divergence
          }
        }
        if (!p_Model_->isEmptyCluster()) break;
        ntry_empty++;
      }
      // all initialization failed
      if(!isinitialized) return false;
      // compute current likelihood
      Lcurrent = p_Model_->estimateLikelihood();
#ifdef COVERBOSE
      std::cout<<"In XStrategy::run(). ixem =" << ixem <<"\n";
      std::cout<<"Initialization terminated. Lcurrent =" << Lcurrent <<"\n";
      if (p_Model_->isEmptyCluster())
      { std::cout<<"Warning: Empty cluster after initialization.\n";}
#endif
      //
      if (Lcurrent>=L1)
      {
        non_empty = true;
#ifdef COVERBOSE
        std::cout<<"In XStrategy::run(). ixem =" << ixem <<"\n";
        p_Model_->consoleOut();
        std::cout<<"L1=" << L1 << "\n";
#endif
        L1 = Lcurrent;
        p_Model_->modifyThetaStart();
      }
    } // ixem
    if(non_empty) // there is a theta start available
    {
      //set epsilon to esp_XEM
      p_Model_->setEpsilon(p_Model_->modelParameters().eps_XEM_);
      p_Model_->copyThetaStart();
      for ( int itr = 0; itr < Stratparam_.nbiter_XEM_; ++itr)
      {
        if(p_Algo_->run())
        {
          (p_Model_->*Stratparam_.Stop_Criteria)();
          if (p_Model_->stopAlgo()) { break;}
        }
        else { break;}
      }
    }
    else { continue;}
    Lcurrent = p_Model_->estimateLikelihood();
    if(!p_Model_->isEmptyCluster()&&Lcurrent>Lmax)
    {
      Lmax = Lcurrent;
      p_Model_->modifyThetaMax();
    }
  } // iTry
  // look at the value of lmax
  if(Lmax > -RealMax)
  {
    p_Model_->copyThetaMax();
    p_Model_->finalizeOutput();
#ifdef COVERBOSE
    std::cout << "XStrategy::run() over.\n";
    std::cout << "Lmax=" << Lmax << "\n";
    p_Model_->consoleOut();
#endif
    return true;
  }
  else
  {
#ifdef COVERBOSE
    std::cout<<"XStrategy::run() over. All  try fails\n";
#endif
    return false;
  }
}
