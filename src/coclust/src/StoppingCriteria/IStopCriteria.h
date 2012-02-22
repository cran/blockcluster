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
 * created on: Dec 27, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file IStopCriteria.h
 *  @brief This file declares abstract class IStopCriteria which provides interface for various stopping criteria.
 **/

#ifndef ISTOPCRITERIA_H_
#define ISTOPCRITERIA_H_

#include "../Models/ICoClustModel.h"
/** @brief IStopCriteria is an abstract class that provides interface for different stopping criteria. Every stoppping criteria
 * must be inherited from this base class. It has only one virtual function IStopCriteria::run which is interface function
 *  to run the stopping criteria.
 */
class IStopCriteria
{
  public:
    /**Constructor*/
    IStopCriteria(){};
    /** This is interface function to run the stopping criteria*/
    virtual void run() = 0;
    /** This function set the pointer IStopCriteria::p_Model_ to given model.*/
    void SetModel(ICoClustModel *& model);
    /**Destructor*/
    virtual ~IStopCriteria(){};

  protected:
    ICoClustModel * p_Model_;

};

inline void IStopCriteria::SetModel(ICoClustModel *& model)
{
  p_Model_ = model;
}
#endif /* ISTOPCRITERIA_H_ */
