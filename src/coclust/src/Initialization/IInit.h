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
 * Project:  cocluster
 * created on: Dec 19, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file IInit.h
 *  @brief Declares abstract class IInit for Initializations.
 **/

#ifndef IINIT_H_
#define IINIT_H_

/** @brief This is an interface class for initialization algorithms. It have only one abstract function IInit::run which
 * when called, initialize the various model parameters by calling the corresponding Initialization function
 * declared in IModel. Every initialization algorithm must derive from this abstract class
 * and provide implementation for IInit::run method.
 */
#include "../Models/ICoClustModel.h"
class IInit
{
  protected:
    /** constructor */
    IInit(){};
    ICoClustModel * p_Model_;
  public:
    /** Interface for running initialization */
    virtual bool run() = 0;
    /*
     * Refer IInit::p_Model_ to the instantiated model class inside CoClustermain.cpp
     */
    void SetModel(ICoClustModel *& model);
    /** Virtual Destructor*/
    virtual ~IInit(){};
};

inline void IInit::SetModel(ICoClustModel *& model)
{
  p_Model_ = model;
}
#endif /* IINIT_H_ */
