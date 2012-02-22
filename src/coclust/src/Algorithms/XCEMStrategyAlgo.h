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

/**
 * @file XCEMStrategyAlgo.h
 * @brief Declares XCEMStrategyAlgo algorithm class.
 */
#ifndef XCEMSTRATEGYALGO_H_
#define XCEMSTRATEGYALGO_H_

/** @brief This is a concrete class derived from IAlgo and declares XCEM strategy.
 * It contains only one function that is required to run XCEM strategy algorithm.
 */

#include "IAlgo.h"

class XCEMStrategyAlgo: public IAlgo
{
  public:
    XCEMStrategyAlgo(AlgoParameters const& Aparam)
    :IAlgo(Aparam)
    {};
    virtual bool run();
    virtual ~XCEMStrategyAlgo(){};
};

#endif /* XCEMSTRATEGYALGO_H_ */
