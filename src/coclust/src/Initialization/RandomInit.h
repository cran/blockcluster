/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

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
 * created on: Feb 3, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file RandomInit.h
 *  @brief Declares Random initialization class RandomInit derived from IInit.
 **/

#ifndef RANDOMINIT_H_
#define RANDOMINIT_H_

/** @brief This class provides functionalities for Random initialization. It is  derived from IInit
 * abstract class.
 *
 */
#include "IInit.h"

class RandomInit: public IInit
{
  public:
    inline RandomInit(){};
    virtual bool run();
    inline virtual ~RandomInit(){};
};

inline bool RandomInit::run()
{
  return p_Model_->RandomInit();
}
#endif /* RANDOMINIT_H_ */
