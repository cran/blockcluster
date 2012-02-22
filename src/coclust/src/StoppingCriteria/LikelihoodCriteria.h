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

/** @file LikelihoodCriteria.h
 *  @brief This file declares and implements class LikelihoodCriteria.
 **/

#ifndef LIKELIHOODCRITERIA_H_
#define LIKELIHOODCRITERIA_H_

#include "IStopCriteria.h"
/** @brief LikelihoodCriteria is an concrete class that calls ICoClustModel::likelihoodStopCriteria() function in its
 * LikelihoodCriteria::run method.
 */
class LikelihoodCriteria : public IStopCriteria
{
  public:
    LikelihoodCriteria(){};
    virtual void run();
    virtual ~LikelihoodCriteria(){};
};

inline void LikelihoodCriteria::run()
{
  p_Model_->likelihoodStopCriteria();
}
#endif /* LIKELIHOODCRITERIA_H_ */
