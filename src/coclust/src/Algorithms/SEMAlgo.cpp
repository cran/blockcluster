/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2013  Parmeet Singh Bhatia

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : bhatia.parmeet@gmail.com , parmeet.bhatia@inria.fr
 */

/** @file SEMAlgo.cpp
 *  @brief In this file .
 **/

#include "SEMAlgo.h"

SEMAlgo::SEMAlgo()
{
  // TODO Auto-generated constructor stub

}

SEMAlgo::~SEMAlgo()
{
  // TODO Auto-generated destructor stub
}

bool SEMAlgo::run()
{
  if (p_Model_->SEMRows()) {
    if (p_Model_->SEMCols()) {
      return true;
    }
  }
  return false;
}
