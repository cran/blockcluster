/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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

    Contact : Serge.Iovleff@stkpp.org
*/

/*
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: implement the Interface base class IModel.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IModel.h
 *  @brief In this file we implement the Interface base class IModel.
 **/

#include "../include/STK_IModel.h"

namespace STK
{
//Constructor
IModel::IModel( Integer const& nbSample, Integer const& nbVar)
              : nbSample_(nbSample)
              , nbVar_(nbVar)
              , logLikelihood_(-Arithmetic<Real>::infinity())
              , nbFreeParameter_(0)
{}

//Destructor
IModel::~IModel() {}

/* set the default value of the parameters of the model */
void IModel::setDefault( Integer const& nbSample, Integer const& nbVar)
{
  nbSample_ = nbSample;
  nbVar_ = nbVar;
  logLikelihood_ = -Arithmetic<Real>::infinity();
  nbFreeParameter_ = 0;
}

} // namespace STK
