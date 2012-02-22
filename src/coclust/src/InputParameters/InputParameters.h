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

/** @file InputParameters.h
 *  @brief This file declares InputParameters class.
 **/

#ifndef INPUTPARAMETERS_H_
#define INPUTPARAMETERS_H_

/** @brief InputParameters contains all the parameters common between all the types of models. All the members in this class
 * are declared static and hence they are global. This class publically inherits IPage class defined in STKpp
 * project http://www.stkpp.org/html/classSTK_1_1IPage.html which allow easy access to reading and writing option values
 * from text file provided by the user.
 */
#include <string>
#include <map>
#include <list>
#include "../typedefs/typedef.h"
#ifdef STK_DMANAGER
#include "../../../stkpp/include/DManager.h"
#endif
#include "../enumerations/enumerations.h"

struct AlgoParameters
{
    //iterations for various stages of XEMStrategy and between EM algo for Rows and Columns
    int nbtry_;
    int nbxem_;
    int nbiter_xem_;
    int nbiter_XEM_;
};

struct ModelParameters
{
    //iterations and epsilon inside model
    float epsilon_int_ ;
    int nbiterations_int_;
    float eps_xem_;
    float eps_XEM_;

    //iterations and epsilon inside initialization
    int nbinititerations_;
    float initepsilon_;

    //espilon set either to eps_xem_ or esp_XEM_ depending on where we are in XEMstrategy
    float epsilon_;

    //other input options, self explanatory
    int nbrowclust_;
    int nbcolclust_;
    int nbrowdata_;
    int nbcoldata_;

    //proportion
    bool fixedproportions_;
};

struct Strategy
{
    Algorithm Algo_;
    StopCriteria stopcriteria_;
    Initialization Init_;
    Model Model_;
    DataType DataType_;
};
#ifdef STK_DMANAGER
class InputParameters : public STK::IPage
{
  public:
    /**Constructor*/
    InputParameters(STK::Integer const& level);
    /** This function reads the value of various parameters from text file.
     * @param optionfilename Path to the input option file.
    */
    void ReadFromOptionFile(std::string optionfilename);
    /**This function initializes enumerations*/
    void InitializeEnum();
    /**Get Model parameters*/
    const ModelParameters& GetModelparameters();
    /**Get Algorithm parameters*/
    const AlgoParameters& GetAlgoparameters();
    /**Get Strategy*/
    const Strategy& GetStrategy();
    /**Get filename*/
    std::string GetDatafilename(){return datafilename_;}
    /**Get optionalfilenames*/
    std::list<std::string>& GetOptionalfilenames(){return optionalfilenames_;}
    void SetModelparameters(ModelParameters const&);
    /**Set Algorithm parameters*/
    void SetAlgoparameters(AlgoParameters const&);
    /**Set strategy*/
    void SetStrategy(Strategy const&);
    /**Destructor*/
    ~InputParameters(){};

  protected:
    ModelParameters Mparam_;
    AlgoParameters Aparam_;
    Strategy strategy_;
    std::map<std::string,Algorithm> S_Algorithm;
    std::map<std::string,StopCriteria> S_StopCriteria;
    std::map<std::string,DataType> S_DataType;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model> S_Model;
    std::string datafilename_;
    std::list<std::string> optionalfilenames_;

};
#else
class InputParameters
{
  public:
    /**Constructor*/
    InputParameters(){};
    /**This function initializes enumerations*/
    void InitializeEnum();
    /**Get Model parameters*/
    const ModelParameters& GetModelparameters();
    /**Get Algorithm parameters*/
    const AlgoParameters& GetAlgoparameters();
    /**Get Strategy*/
    const Strategy& GetStrategy();
    /**Get filename*/
    std::string GetDatafilename(){return datafilename_;}
    /**Get optionalfilenames*/
    std::list<std::string>& GetOptionalfilenames(){return optionalfilenames_;}
    /**Set Model parameters*/
    void SetModelparameters(ModelParameters const&);
    /**Set Algorithm parameters*/
    void SetAlgoparameters(AlgoParameters const&);
    /**Set strategy*/
    void SetStrategy(Strategy const&);
    /**Destructor*/
    ~InputParameters(){};

  protected:
    ModelParameters Mparam_;
    AlgoParameters Aparam_;
    Strategy strategy_;
    std::map<std::string,Algorithm> S_Algorithm;
    std::map<std::string,StopCriteria> S_StopCriteria;
    std::map<std::string,DataType> S_DataType;
    std::map<std::string,Initialization> S_Init;
    std::map<std::string,Model> S_Model;
    std::string datafilename_;
    std::list<std::string> optionalfilenames_;

};
#endif

inline void InputParameters::InitializeEnum()
{
  //Datatype
  S_DataType["Binary"] = Binary;
  S_DataType["Contingency"] = Contingency;
  S_DataType["Continuous"] = Continuous;

  //Algorithm
  S_Algorithm["XEMStrategy"] = XEMStrategy;
  S_Algorithm["BEM2"] = BEM2;
  S_Algorithm["XCEMStrategy"] = XCEMStrategy;
  S_Algorithm["BCEM"] = BCEM;

  //StopCriteria
  S_StopCriteria["Parameter"] = Parameter;
  S_StopCriteria["Likelihood"] = Likelihood;

  //Initialization
  S_Init["CEMInit"] = e_CEMInit;
  S_Init["FuzzyCEMInit"] = e_FuzzyCEMInit;
  S_Init["RandomInit"] = e_RandomInit;

  //Models
  S_Model["pi_rho_epsilon"] = pi_rho_epsilon;
  S_Model["pik_rhol_epsilon"] = pik_rhol_epsilon;
  S_Model["pi_rho_epsilonkl"] = pi_rho_epsilonkl;
  S_Model["pik_rhol_epsilonkl"] = pik_rhol_epsilonkl;
  S_Model["pi_rho_unknown"] = pi_rho_unknown;
  S_Model["pik_rhol_unknown"] = pik_rhol_unknown;
  S_Model["pi_rho_known"] = pi_rho_known;
  S_Model["pik_rhol_known"] = pik_rhol_known;
  S_Model["pi_rho_sigma2"] = pi_rho_sigma2;
  S_Model["pik_rhol_sigma2"] = pik_rhol_sigma2;
  S_Model["pi_rho_sigma2kl"] = pi_rho_sigma2kl;
  S_Model["pik_rhol_sigma2kl"] = pik_rhol_sigma2kl;
}

inline const ModelParameters& InputParameters::GetModelparameters()
{
  return Mparam_;
}

inline const AlgoParameters& InputParameters::GetAlgoparameters()
{
  return Aparam_;
}

inline const Strategy& InputParameters::GetStrategy()
{
  return strategy_;
}

inline void InputParameters::SetAlgoparameters(AlgoParameters const& Aparam)
{
  Aparam_ = Aparam;
}

inline void InputParameters::SetModelparameters(ModelParameters const& Mparam)
{
  Mparam_ = Mparam;
}

inline void InputParameters::SetStrategy(Strategy const& strat)
{
  strategy_ = strat;
}
#endif /* INPUTPARAMETERS_H_ */
