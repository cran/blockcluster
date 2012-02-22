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

/**
 * @file XEMStrategyAlgo.cpp
 * @brief Implements XEM strategy algorithm class XEMStrategyAlgo derived from IAlgo.
 */

#include "XEMStrategyAlgo.h"

bool XEMStrategyAlgo::run()
{
#ifdef COVERBOSE
	std::cout<<"Running XEM Strategy Algorithm.."<<"\n";
#endif
	float Lmax = -RealMax,L1 = -RealMax,Lcurrent;
	int ntry_empty=0;
	bool non_empty = false;
	for (int itry = 0; itry < Aparam_.nbtry_; ++itry) {
		non_empty = false;
		L1 = -RealMax;
		//set espilon to eps_xem
		p_Model_->SetEpsilon(0);
		for (int ixem = 0; ixem < Aparam_.nbxem_; ++ixem) {
			ntry_empty = -1;
			p_Model_->SetEmptyCluster(true);
			while(p_Model_->isEmptyCluster()&&ntry_empty<100)
			{
				if (p_Init_->run())
				{
				    for (int itr = 0; itr < Aparam_.nbiter_xem_; ++itr)
				    {

				      if(p_Model_->Estep())
				      {
					      p_Model_->Mstep();
					      p_StopCriteria_->run();
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
				    p_Model_->FinalizeOutput();
				}
				ntry_empty++;
			}
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
		  p_Model_->SetEpsilon(1);
			p_Model_->Copy_theta_start();
		    for (int itr = 0; itr < Aparam_.nbiter_XEM_; ++itr)
		    {

		      if(p_Model_->Estep())
		      {
			      p_Model_->Mstep();
			      p_StopCriteria_->run();
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
		    p_Model_->FinalizeOutput();
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
