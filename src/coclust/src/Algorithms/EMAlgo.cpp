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
 * created on: Nov 25, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file EMAlgo.cpp
 *  @brief Implements EMAlgo class derived from IAlgo.
 **/

#include "EMAlgo.h"
#include "../Models/ICoClustModel.h"
/**
 * The run method do the following in sequence:
 * -# It initialized the model parameters using ICoClustModel::Initialize interface.
 * -# It iterates between ICoClustModel::Estep and ICoClustModel::Mstep until either maximum number of iterations Aparam_.nbiterations_ is reached
 * or ICoClustModel::StopAlgo condition becomes true.
 * -# It finalizes the model output using ICoClustModel::Finalize interface.
 */
bool EMAlgo::run()
{
    if (p_Init_->run()) {
      p_Model_->SetEpsilon(1);
#ifdef COVERBOSE
    std::cout<<"Running Algorithm.."<<std::endl;
#endif
    bool flag = false , clusteringfail = false;
    for (int itr = 0; itr < Aparam_.nbiter_XEM_; ++itr)
    {

      if(!p_Model_->Estep())
      {
#ifdef COVERBOSE
        std::cout<<"Global iterations:"<<itr<<"\n";
#endif
        clusteringfail = true;
        break;
      }
      p_Model_->Mstep();
      p_StopCriteria_->run();
      if(p_Model_->stopAlgo())
      { flag = true;
        break;
      }
    }

    if (!clusteringfail) {
#ifdef COVERBOSE
    std::cout<<"Finalizing output.."<<std::endl;
#endif
    p_Model_->FinalizeOutput();
#ifdef COVERBOSE
    p_Model_->ConsoleOut();
#endif
#ifdef COVERBOSE
    if (!flag) std::cout<<"Algorithm over(Maximum iterations reached)."<<std::endl;
    else  std::cout<<"Algorithm stop successfully by epsilon criteria!"<<std::endl;
#endif
#ifdef RPACKAGE
    if (!flag) p_Model_->SetRmessage("Algorithm over(Maximum iterations reached");
    else  p_Model_->SetRmessage("Algorithm stop successfully by epsilon criteria!");
#endif
    return true;
    }
    else
      return false;
    }
#ifdef COVERBOSE
    std::cout<<"Initialization Failed.."<<std::endl;
#endif
    return false;
}
