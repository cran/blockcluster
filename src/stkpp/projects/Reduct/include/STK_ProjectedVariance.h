/*
 * STK_ProjectedVariance.h
 *
 *  Created on: 20 juil. 2011
 *      Author: aude
 */

#ifndef STK_PROJECTEDVARIANCE_H_
#define STK_PROJECTEDVARIANCE_H_

#include "STK_ILinearReduct.h"

#include "../../Algebra/include/STK_Svd.h"




namespace STK
{

/** @ingroup Reduct
 *  @brief A ProjectedVariance is an implementation of the abstract
 *  @c ILinearReduct class.
 *
 *  ProjectedVariance (PCA) is the best method to reduce the dimension of data.
 *  This method computes the principal components.
 *  The number of principal components is the same of the number of variables.
 *  Then, the best components are selected thanks to the index of inertia.
**/


class ProjectedVariance : public ILinearReduct
{
  public:


    ProjectedVariance(Matrix const* p_data);
    /**
     * Destructor
     */
    ~ProjectedVariance();


    /** Compute the Index.
     *  @param dim the dimension of the reduced data set
     */
    virtual void run( Integer const& dim );
    /**
     * Compute the weighted index.
     * @param weights the weights to used
     * @param dim number of Axis to compute
     */
    virtual void run( Vector const* weights, Integer const& dim);

    /** Do svd on the matrix p_data_**/
    void svdData();

    /** Do svd on the matrix p_data_ with the weights of the samples **/
    void svdWData();


  protected:


    /** Find the axis by maximizing the Index. */
    virtual void maximizeIndex();

    /** Find the axis by maximizing the weighed Index. */
    virtual void wmaximizeIndex();

    /** Compute the new coordinates of samples **/
    void dataReduced();

    /** Object which contains the matrix resulting from svd of p_data_ **/
    Svd svdData_;

  private:

    /** Compute the axis using the matrix V_ of the svd of data **/
    void computeAxis();


};

} // namespace STK

#endif /* STK_PROJECTEDVARIANCE_H_ */
