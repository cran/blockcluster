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
 * created on: Dec 22, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file InputParameters.cpp
 *  @brief This file implement methods of InputParameters class.
 **/

#include "InputParameters.h"
#ifdef STK_DMANAGER
InputParameters::InputParameters(STK::Integer const& level=1) : STK::IPage("InputOptions", level, true)
{
  InitializeEnum();
   options_.reserve(20);
   options_.push_back(STK::Option("DataType", STK::Option::string_, false));

   options_.push_back(STK::Option("nbinititerations", STK::Option::integer_, true));
   options_.push_back(STK::Option("initepsilon", STK::Option::real_, true));

   options_.push_back(STK::Option("nbtry", STK::Option::integer_, true));
   options_.push_back(STK::Option("nbxem", STK::Option::integer_, true));
   options_.push_back(STK::Option("epsilon_xemstart", STK::Option::real_, true));
   options_.push_back(STK::Option("epsilon_xem", STK::Option::real_, true));
   options_.push_back(STK::Option("nbiterations_xemstart", STK::Option::integer_, true));
   options_.push_back(STK::Option("nbiterations_xem", STK::Option::integer_, true));

   options_.push_back(STK::Option("nbiterations_int", STK::Option::integer_, true));
   options_.push_back(STK::Option("epsilon_int", STK::Option::real_, true));

   options_.push_back(STK::Option("nbrowclust", STK::Option::integer_, false));
   options_.push_back(STK::Option("nbcolclust", STK::Option::integer_, false));
   //options_.push_back(STK::Option("fixedproportions", STK::Option::integer_,true));
   options_.push_back(STK::Option("Algorithm", STK::Option::string_, true));
   options_.push_back(STK::Option("ModelName", STK::Option::string_, true));
   options_.push_back(STK::Option("StopCriteria", STK::Option::string_, true));
   //options_.push_back(STK::Option("sparse", STK::Option::integer_, true));
   options_.push_back(STK::Option("Initialization", STK::Option::string_, true));
   options_.push_back(STK::Option("DataFileName", STK::Option::string_, false));
   options_.push_back(STK::Option("OptionalFileNames", STK::Option::lstring_, true));
}

void InputParameters::ReadFromOptionFile(std::string optionfilename)
{

  STK::ReadWritePages rw(optionfilename);
  InputParameters ip_page;
  rw.addPage(ip_page);
  rw.read(optionfilename);
  strategy_.DataType_ = S_DataType[rw.p_page("InputOptions")->option("DataType").get(STK::String())];

  Mparam_.nbinititerations_ = rw.p_page("InputOptions")->option("nbinititerations").get(STK::Integer());
  Mparam_.initepsilon_ = rw.p_page("InputOptions")->option("initepsilon").get(STK::Real());

  Aparam_.nbtry_ = rw.p_page("InputOptions")->option("nbtry").get(STK::Integer());
  Aparam_.nbxem_ = rw.p_page("InputOptions")->option("nbxem").get(STK::Integer());
  Mparam_.eps_xem_ = rw.p_page("InputOptions")->option("epsilon_xemstart").get(STK::Real());
  Mparam_.eps_XEM_ = rw.p_page("InputOptions")->option("epsilon_xem").get(STK::Real());
  Aparam_.nbiter_xem_ = rw.p_page("InputOptions")->option("nbiterations_xemstart").get(STK::Integer());
  Aparam_.nbiter_XEM_ = rw.p_page("InputOptions")->option("nbiterations_xem").get(STK::Integer());

  Mparam_.epsilon_int_ = rw.p_page("InputOptions")->option("epsilon_int").get(STK::Real());
  Mparam_.nbiterations_int_ = rw.p_page("InputOptions")->option("nbiterations_int").get(STK::Integer());

  strategy_.Algo_ = S_Algorithm[rw.p_page("InputOptions")->option("Algorithm").get(STK::String())];
  strategy_.Model_ = S_Model[rw.p_page("InputOptions")->option("ModelName").get(STK::String())];
  strategy_.stopcriteria_ = S_StopCriteria[rw.p_page("InputOptions")->option("StopCriteria").get(STK::String())];
  Mparam_.nbrowclust_ = rw.p_page("InputOptions")->option("nbrowclust").get(STK::Integer());
  Mparam_.nbcolclust_ = rw.p_page("InputOptions")->option("nbcolclust").get(STK::Integer());
  //Mparam_.fixedproportions_ = rw.p_page("InputOptions")->option("fixedproportions").get(STK::Integer());
  strategy_.Init_ = S_Init[rw.p_page("InputOptions")->option("Initialization").get(STK::String())];
  datafilename_ = rw.p_page("InputOptions")->option("DataFileName").get(STK::String());
  optionalfilenames_ = rw.p_page("InputOptions")->option("OptionalFileNames").get(std::list<std::string>());

  //set fixproportions
  switch (strategy_.Model_) {
    case pik_rhol_epsilonkl:
      Mparam_.fixedproportions_ = false;
      break;
    case pik_rhol_epsilon:
      Mparam_.fixedproportions_ = false;
      break;
    case pi_rho_epsilon:
      Mparam_.fixedproportions_ = true;
      break;
    case pi_rho_epsilonkl:
      Mparam_.fixedproportions_ = true;
      break;
    case pik_rhol_sigma2:
      Mparam_.fixedproportions_ = false;
      break;
    case pik_rhol_sigma2kl:
      Mparam_.fixedproportions_ = false;
      break;
    case pi_rho_sigma2:
      Mparam_.fixedproportions_ = true;
      break;
    case pi_rho_sigma2kl:
      Mparam_.fixedproportions_ = true;
      break;
    case pik_rhol_unknown:
      Mparam_.fixedproportions_ = false;
      break;
    case pik_rhol_known:
      Mparam_.fixedproportions_ = false;
      break;
    case pi_rho_unknown:
      Mparam_.fixedproportions_ = true;
      break;
    case pi_rho_known:
      Mparam_.fixedproportions_ = true;
      break;
    default:
      Mparam_.fixedproportions_ = false;
      break;
  }

  //rw.write(std::cout);
}
#endif
