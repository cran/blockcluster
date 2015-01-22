/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2013  Parmeet Singh Bhatia

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

 Contact : bhatia.parmeet@gmail.com , parmeet.bhatia@inria.fr
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
  std::cout << classes.cols() << std::endl;
  std::cout << classes.rows() << std::endl;
  for (int i = 0; i < classes.rows(); ++i)
  {
    int maxInd = 0;
    int maxVal = classes(i, 0);
    for (int j = 0; j < classes.cols(); ++j)
    {
      if (classes(i, j) > maxVal) maxInd = j;
    }
    o_file << maxInd << std::endl;
  }
}

bool XStrategyAlgo::run()
{
#ifdef COVERBOSE
  std::cout<<"Running XStrategy.."<<"\n";
#endif
  float Lmax = -RealMax,L1 = -RealMax,Lcurrent;
  int ntry_empty=0;
  bool non_empty = false;
  bool isinitialized = false;
  for ( int itry = 0; itry < Stratparam_.nbtry_; ++itry) {
    non_empty = false;
    L1 = -RealMax;
    //set espilon to eps_xem
    p_Model_->SetEpsilon(p_Model_->GetModelParameters().eps_xem_);
    for ( int ixem = 0; ixem < Stratparam_.nbxem_; ++ixem) {
      ntry_empty = -1;
      p_Model_->SetEmptyCluster(true);
      while(p_Model_->isEmptyCluster()&&ntry_empty<100)
      {
        if (p_Init_->run())
        {
          isinitialized = true;
            for ( int itr = 0; itr < Stratparam_.nbiter_xem_; ++itr)
            {

              if(p_Algo_->run())
              {
                (p_Model_->*Stratparam_.Stop_Criteria)();
                if(p_Model_->stopAlgo())
                {
                   break;
                }
              }
              else
              {
                break;
              }
            }
            //p_Model_->FinalizeOutput();
        }
        ntry_empty++;
      }

      if(!isinitialized) return false;
      Lcurrent = p_Model_->EstimateLikelihood();
      if(!p_Model_->isEmptyCluster()&&Lcurrent>=L1)
      {
        non_empty = true;
        L1 = Lcurrent;
        p_Model_->Modify_theta_start();
      }
    }
    //set epsilon to esp_XEM
    if(non_empty){
      p_Model_->SetEpsilon(p_Model_->GetModelParameters().eps_XEM_);
      p_Model_->Copy_theta_start();
        for ( int itr = 0; itr < Stratparam_.nbiter_XEM_; ++itr)
        {

          if(p_Algo_->run())
          {
            (p_Model_->*Stratparam_.Stop_Criteria)();
            if(p_Model_->stopAlgo())
            {
              break;
            }
          }
          else
          {
            break;
          }
        }
        //p_Model_->FinalizeOutput();
    }
    Lcurrent = p_Model_->EstimateLikelihood();
    if(!p_Model_->isEmptyCluster()&&Lcurrent>Lmax)
    {
      Lmax = Lcurrent;
      p_Model_->Modify_theta_max();
    }
  }
  if(!p_Model_->isEmptyCluster()){
    p_Model_->Copy_theta_max();
    p_Model_->FinalizeOutput();
#ifdef COVERBOSE
    std::cout<<"Algorithm over."<<"\n";
    p_Model_->ConsoleOut();
    std::cout<<"\nLmax:"<<Lmax<<"\n";
#endif
    return true;
  }
  else
  {
    return false;
  }
}
