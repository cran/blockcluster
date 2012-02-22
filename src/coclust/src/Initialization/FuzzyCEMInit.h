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

/** @file FuzzyCEMInit.h
 *  @brief Declares Fuzzy CEM initialization class FuzzyCEMInit derived from IInit.
 **/

#ifndef FUZZYCEMINIT_H_
#define FUZZYCEMINIT_H_

/** @brief This class provides functionalities for Fuzzy CEM initialization. It is  derived from IInit
 * abstract class.
 *
 */
#include "IInit.h"

class FuzzyCEMInit: public IInit
{
  public:
    inline FuzzyCEMInit(){};
    virtual bool run();
    inline virtual ~FuzzyCEMInit(){};
};

inline bool FuzzyCEMInit::run()
{
  return p_Model_->FuzzyCEMInit();
}
#endif /* FUZZYCEMINIT_H_ */
