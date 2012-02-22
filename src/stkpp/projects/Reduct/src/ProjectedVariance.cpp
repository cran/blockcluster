///*
// * ProjectedVariance.cpp
// *
// *  Created on: 26 juil. 2011
// *      Author: aude
// */
//
//#include "../include/STK_ProjectedVariance.h"
//
//#include "../../Algebra/include/STK_Svd.h"
//#include "../../Algebra/include/STK_LinAlgebra1D.h"
//#include "../../Algebra/include/STK_LinAlgebra2D.h"
//
//#include "../../STKernel/include/STK_Misc.h"
//#include "../../STKernel/include/STK_ITContainer1D.h"
//
//#include <cmath>
//
//#include <iostream>
//#include <vector>
//
//
//
//
//namespace STK
//{
///*
// * Constructor.
// * @param data the input data set
// */
//ProjectedVariance::ProjectedVariance(Matrix const* p_data)
//                                    : ILinearReduct(p_data)
//{
//
//}
//
///*
// * Destructor
// */
//ProjectedVariance::~ProjectedVariance()
//{}
//
//
//void ProjectedVariance::run( Integer const& dim)
//{
//  dim_ = dim;
//
//  svdData();
//
//  computeAxis();
//
//  maximizeIndex();
//
//  dataReduced();
//}
//
//
//void ProjectedVariance::run( Vector const* weights, Integer const& dim)
//{
//  dim_ = dim;
//
//  p_weights_ = weights;
//
//  svdWData();
//
//  computeAxis();
//
//  wmaximizeIndex();
//
//  dataReduced();
//}
//
//
//void ProjectedVariance::svdData()
//{
//  svdData_ = Svd(*p_data_);
//}
//
//void ProjectedVariance::svdWData()
//{
//  svdData_ = Svd(*p_data_);
//  Matrix const* p_w_data;
//  for(int i=0;i<svdData_.nrowU_;i++)
//  {
//    for(int j=0;j<svdData_.ncolU_;j++)
//    {
//      (*p_w_data)(i,j) = (*p_data_)(i,j) * sqrt((*p_weights_)[i]);
//    }
//  }
//  svdData_ = Svd(*p_w_data);
//}
//
//
//void ProjectedVariance::computeAxis()
//{
//  for(int i=0;i<svdData_.ncolV_;i++)
//  {
//    for(int j=0;j<dim_;j++)
//    {
//      axis_(i,j) = svdData_.V_(i,j);
//    }
//  }
//}
//
//
//
//void ProjectedVariance::dataReduced()
//{
//  /** We determine the dim_ first vectors of the final matrix of data (after the dimension reduction) **/
//  for(int i=0;i<dim_;i++)
//  {
//    /** We look down each row of the matrix p_data_ **/
//    //for(int j=0;j<p_data_->sizeHo();j++)
//    for(int j=0;j<svdData_.nrowU_;j++)
//    {
//      (*p_reduced_)(j,i) = 0;
//      /** We look down each column of the matrix p_data and in parallel each column of the matrix V **/
//      //for(int k=0;k<p_data_->sizeVe();k++)
//      for(int k=0;k<svdData_.ncolU_;k++)
//      {
//          (*p_reduced_)(j,i) = (*p_reduced_)(j,i) + (*p_data_)(j,k) * axis_(i,k);
//      }
//    }
//  }
//}
//
//
//
//
//void ProjectedVariance::maximizeIndex()
//{
//  for(int i=0;i<dim_;i++)
//  {
//    index_values_[i] = svdData_.D_(i,i) * svdData_.D_(i,i);
//  }
//}
//
//
//
//void ProjectedVariance::wmaximizeIndex()
//{
//  for(int i=0;i<dim_;i++)
//  {
//    index_values_[i] = svdData_.D_(i,i) * svdData_.D_(i,i);
//  }
//}
//
//
//} // namespace STK
