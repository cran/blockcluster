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


/** @file IStrategy.h
 *  @brief Declares interface class IStrategy for algorithms.
 **/

#ifndef ISTRATEGY_H_
#define ISTRATEGY_H_

#include "../Models/ICoClustModel.h"
#include "../Initialization/IInit.h"
#include "../Algorithms/IAlgo.h"

/** @brief This is is interface(abstract) class for running various algorithms.
 */
class IStrategy
{
  public:
    /** Constructor*/
    inline IStrategy(StrategyParameters const& Aparam)
    {
      Stratparam_ = Aparam;
    }
    /** Virtual destructor*/
    inline virtual ~IStrategy(){};
    /** Interface for running the algorithm.
     *  @return @c false an erro occur, @c false otherwise
     **/
    virtual bool run() = 0;
    /** Function to set the model
     *  @param model the model to set
     **/
    void SetModel(ICoClustModel * model);
    void SetInit(IInit * init);
    void SetAlgo(IAlgo* algo);

  protected:
    ICoClustModel * p_Model_; /**<Pointer to Model to be run.*/
    IInit * p_Init_; /**Pointer to initialization to be run*/
    IAlgo * p_Algo_; /**Pointer to algorithm to be run*/
    StrategyParameters Stratparam_;
};


inline void IStrategy::SetModel(ICoClustModel * model)
{
  p_Model_ = model;
}

inline void IStrategy::SetInit(IInit * init)
{
  p_Init_ = init;
}

inline void IStrategy::SetAlgo(IAlgo *  algo)
{
  p_Algo_ = algo;
}

#endif /* ISTRATEGY_H_ */
