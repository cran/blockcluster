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

/** @file CoCluster.h
 *  @brief Declares class CoCluster which acts as entry point to do co-clustering.
 */

#ifndef COCLUSTER_H_
#define COCLUSTER_H_

#include "../Algorithms/IAlgo.h"
class IInit;
class ICoClustModel;
class IStopCriteria;


/**@brief This class acts as entry point to perform co-clustering. This is the only class which is instantiated inside
 * the main function. This class provides functions for setting algorithm, model and stopping criteria. It also implements
 * display co-clusters function in case of binary data. Refer to file CoClustermain.cpp to see how it works.
 */
class CoCluster
{
  public:
    /**Constructor*/
    inline CoCluster(){};
    /** This function perform Co-Clustering.*/
    bool run();
    /**This function sets the algorithm to be run.*/
    void SetAlgo(IAlgo *& Algo);
    /**This function sets the model to be run.*/
    void SetModel(ICoClustModel *& Model);
    void SetInit(IInit *& init);
    /**This function sets the stopping criteria to be used*/
    void SetStopCriteria(IStopCriteria *& stop);
    ~CoCluster(){};

  private:
    /** Main pointer on the algorithm instantiated in the main function or by the
     *  user of this class.
     **/
    IAlgo * p_Algo_;
    ICoClustModel * p_Model_;
    IInit * p_Init_;
    IStopCriteria * p_StopCriteria_;
};

inline void CoCluster::SetAlgo(IAlgo *& Algo)
{
  p_Algo_ = Algo;
}

inline void CoCluster::SetStopCriteria(IStopCriteria *& stop)
{
  p_StopCriteria_ = stop;
}

inline void CoCluster::SetModel(ICoClustModel *& Model)
{
  p_Model_ = Model;
}

inline void CoCluster::SetInit(IInit *& init)
{
  p_Init_ = init;
}
#endif /* COCLUSTER_H_ */
