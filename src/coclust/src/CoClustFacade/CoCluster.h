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

/** @file CoCluster.h
 *  @brief Declares class CoCluster which acts as entry point to do co-clustering.
 */

#ifndef COCLUSTER_H_
#define COCLUSTER_H_

#include "../Strategy/IStrategy.h"
class IInit;
class ICoClustModel;
class IAlgo;


/**@brief This class acts as entry point to perform co-clustering. This class provides functions for setting
 * algorithm, model and Initialization. This class also provides a CoCluster::run method that sets various pointers
 * in various classes and finally call the IAlgo::run method to perform co-clustering.
 */
class CoCluster
{
  public:
    /**Default Constructor*/
    inline CoCluster(){};
    /** This function perform Co-Clustering.*/
    bool run();
    /**It sets the algorithm to be run.*/
    void SetStrategy(IStrategy * );
    /**It  sets the model to be run.*/
    void SetModel(ICoClustModel * );
    /** It sets the Initialization method */
    void SetInit(IInit * );
    /** It sets the algorithm to be run */
    void SetAlgo(IAlgo * );
    /**Default Destructor*/
    ~CoCluster(){};

  private:
    IStrategy * p_Strategy_; /**Pointer to algorithm*/
    ICoClustModel * p_Model_; /**Pointer to Model*/
    IInit * p_Init_; /**Pointer to Initialization*/
    IAlgo * p_Algo_; /**Pointer to Algorithm*/
};

inline void CoCluster::SetStrategy(IStrategy * strategy)
{
  p_Strategy_ = strategy;
}

inline void CoCluster::SetModel(ICoClustModel * Model)
{
  p_Model_ = Model;
}

inline void CoCluster::SetInit(IInit * Init)
{
  p_Init_ = Init;
}

inline void CoCluster::SetAlgo(IAlgo * algo)
{
  p_Algo_ = algo;
}
#endif /* COCLUSTER_H_ */
