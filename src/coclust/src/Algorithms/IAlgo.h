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

/** @file IAlgo.h
 *  @brief Declares interface class IAlgo for algorithms.
 **/

#ifndef IALGO_H_
#define IALGO_H_

#include "../Models/ICoClustModel.h"
#include "../StoppingCriteria/IStopCriteria.h"
#include "../Initialization/IInit.h"

#include "../../../stkpp/include/Sdk.h"

/** @brief This is is interface(abstract) class for algorithms. This class provides
 * interface function for running
 * the algorithm, for setting the model and stopping criteria.
 *
 * It derive from the interface base class IRunnerBase. If some
 * error occur in the run method. It will be catch and the message stored
 * in msg_eror_.
 */
class IAlgo : public STK::IRunnerBase
{
  public:
    /** Constructor*/
    inline IAlgo(AlgoParameters const& Aparam)
    {
      Aparam_ = Aparam;
    }
    /** Virtual destructor*/
    inline virtual ~IAlgo(){};
    /** Interface for running the algorithm.
     *  @return @c false an erro occur, @c false otherwise
     **/
    virtual bool run() = 0;
    /** Function to set the model
     *  @param model the model to set
     **/
    void SetModel(ICoClustModel *& model);
    /** Function to set stopping criteria
     * @param stop the stopping criteria to set
     **/
    void SetStopCriteria(IStopCriteria *& stop);
    /*
     * Function to set Initialization
     */
    void SetInitialization(IInit *& init);

  protected:
    ICoClustModel * p_Model_; /**<Pointer to Model to be run.*/
    IStopCriteria * p_StopCriteria_;/**<Pointer to stopping criteria.*/
    IInit * p_Init_;/**<Pointer to initialization.*/
    AlgoParameters Aparam_;
};


inline void IAlgo::SetModel(ICoClustModel *& model)
{
  p_Model_ = model;
}

inline void IAlgo::SetStopCriteria(IStopCriteria *& stop)
{
  p_StopCriteria_ = stop;
}

inline void IAlgo::SetInitialization(IInit *& init)
{
  p_Init_ = init;
}
#endif /* IALGO_H_ */
