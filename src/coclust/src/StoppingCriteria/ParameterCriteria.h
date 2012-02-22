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

/** @file ParameterCriteria.h
 *  @brief This file declares and implements class ParameterCriteria.
*/
#ifndef PARAMETERCRITERIA_H_
#define PARAMETERCRITERIA_H_

#include "IStopCriteria.h"
/** @brief ParameterCriteria is an concrete class that calls ICoClustModel::ParameterStopCriteria() function in its
 * ParameterCriteria::run method.
 *
 */
class ParameterCriteria : public IStopCriteria
{
  public:
    ParameterCriteria(){};
    virtual void run();
    virtual ~ParameterCriteria(){};
};

inline void ParameterCriteria::run()
{
  p_Model_->ParameterStopCriteria();
}

#endif /* PARAMETERCRITERIA_H_ */
