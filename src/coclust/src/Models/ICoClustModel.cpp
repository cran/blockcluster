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
 * created on: Feb 22, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file ICoClustModel.cpp
 *  @brief This file only initializes the static members of ICoClustModel.
 **/

#include "ICoClustModel.h"


VectorInteger ICoClustModel::RandSample(int n,int k)
{
  //random shuffle Algorithm
  srand(float(clock()));
  int random,temp;
  VectorInteger v_temp(n),v_randint(k);
  for (int j = 0; j < n; ++j) {
    v_temp(j)=j;
  }

  for (int l = 0; l < k; ++l){
    random=std::rand()%(n-l);
    v_randint(l) = v_temp(random);
    //swap elements
    temp = v_temp(n-l-1);
    v_temp(n-l-1)=v_temp(random);
    v_temp(random)=temp;
  }
  return v_randint;
}

