
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
 * created on: Dec 14, 2011
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file typedef.h
 *  @brief This file define all the typedefs used in cocluster project.
 **/


#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <limits>

#ifdef EIGENCONTAINERS
#include "../../../Eigen/Dense"
#else
#include "../../../stkpp/include/Arrays.h"
#endif

class InputParameters;

/*
 * definition for numeric limits
 */
#define RealMax std::numeric_limits<float>::max()
#define RealMin std::numeric_limits<float>::min()

/*
 * Typedefs for Matrix and vector containers
 */

//Matrix containers
typedef Eigen::MatrixXf MatrixReal;
typedef Eigen::MatrixXi MatrixInteger;
typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixBinary;
//Vector Containers
typedef Eigen::VectorXf VectorReal;
typedef Eigen::VectorXi VectorInteger;
typedef Eigen::Matrix<bool,Eigen::Dynamic,1> VectorBinary;
//2D array containers
typedef Eigen::ArrayXXf Array2DReal;

/*
 * Alias for some classes
 */

typedef InputParameters IP;

#endif /* TYPEDEF_H_ */
