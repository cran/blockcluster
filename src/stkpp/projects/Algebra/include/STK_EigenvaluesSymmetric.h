/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  stkpp::Algebra
 * Purpose:  Define The EigenvaluesSymmetric Class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_EigenvaluesSymmetric.h
 *  @brief In this file we define the EigenvaluesSymmetric class (for a
 * symmetric matrix).
 **/
 
#ifndef STK_EIGENVALUESSYMMETRIC_H
#define STK_EIGENVALUESSYMMETRIC_H

#include "../../Sdk/include/STK_IRunnerPtr2D.h"

#include "../../Arrays/include/STK_MatrixSquare.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief The class EigenvaluesSymmetric compute the eigenvalue Decomposition
 *  of a symmetric Matrix.
 * 
 *  The decomposition of a symmetric matrix require
 *  - Input:  A symmetric matrix A of size (n,n)
 *  - Output:
 *     -# P Matrix of size (n,n).
 *     -# D Vector of dimension n
 *     -# \f$ A = PDP' \f$
 *  The matrix A can be copied or overwritten by the class.
 *
 *  The 2-norm (operator norm) of the matrix is given. if the 2-norm is less
 *  than the arithmetic precision of the type @c Real, the rank is set to 0.
 *  Thus the user can be faced with a rank 0 matrix and with a norm and a
 *  determinant very small.
 **/
class EigenvaluesSymmetric : public IRunnerPtr2D<Real, Matrix>
{
  /** concrete type of the runner */
  typedef IRunnerPtr2D<Real, Matrix> IRunnerPtrMatrix;

  public:
    /** Default constructor.
     *  @param A The symmetric matrix to decompose.
     **/
    EigenvaluesSymmetric( MatrixSquare const* A);
    
    /** Copy constructor.
     *  @param S the EigenValue to copy
     **/
    EigenvaluesSymmetric(const EigenvaluesSymmetric& S);

    /** virtual destructor */
    virtual ~EigenvaluesSymmetric();

    /** @brief Diagonalization of P_
     *  @return @c true if no error occur, @c false otherwise
     * */
    virtual bool run();

    /** @brief weighted diagonalization of P_
     *  @param weights the weights of each rows
     *  @return @c true if no error occur, @c false otherwise
     * */
    virtual bool run(Vector const& weights);

    /** Operator = : overwrite the EigenvaluesSymmetric with S.
     *  @param S EigenvaluesSymmetric to copy
     *  @return a reference on this
     **/
    EigenvaluesSymmetric& operator=(const EigenvaluesSymmetric &S);

    /** norm of the matrix
     * @return the norm of the matrix
     **/
    inline Real norm()  const { return norm_;}
    /** rank of the matrix
     *  @return the rank of the matrix
     **/
    inline Integer rank()  const { return rank_;}
    /** determinant of the Matrix
     * @return the determinant of the Matrix
     **/
    inline Real det()  const { return det_;}
    /** get rotation matrix
     *  @return the rotation matrix
     **/
    inline MatrixSquare const& rotation() const{ return P_;}
    /** get eigenvalues
     *  @return the eigenvalues
     **/
    inline Vector const& eigenvalues() const { return D_;}
    /** get the first index of the rows/columns
     *  @return the index of the first row/column
     **/
    inline Integer const& first() const{ return first_;}
    /** get last index of the rows/columns
     *  @return the index of the last row/column
     **/
    inline Integer const& last() const { return last_;}

    /** Compute the generalized inverse of the symmetric matrix.
     * The result is allocated dynamically and is not liberated by this
     * class.
     * @return A pointer on the generalized inverse.
     */
    MatrixSquare* ginv();

    /** Compute the generalized inverse of the symmetric matrix and put
     *  the result in res.
     *  @param res the generalized inverse of the Matrix.
     */
    void ginv(MatrixSquare& res);

  protected:
    /** P_ Square matrix or the eigenvectors. P_ is initialized using the
     *  matrix passe to the constructor.
     *  The initialization can be done by reference, in this case A will be
     *  overwritten in the diagonalization process.
     */
    MatrixSquare P_;
    /// Array of the eigenvalues
    Vector D_;
    /// first row/col of P_
    Integer first_;
    /// last row/col of P_
    Integer last_;
    /// norm of the matrix
    Real norm_;
    /// rank of the matrix
    Integer rank_;
    /** determinant of the Matrix */
    Real det_;
    /// Everything OK during computation ?
    bool error_;

  private:
    /** Temporary vector. Values of the sub-diagonal. */
    Vector F_;
    /** @brief compute the tri-diagonalization of P_ */
    void tridiagonalize();
    /** @brief compute the Householder matrix and P */
    void compHouse();
    /** computing the diagonalization of D_ and F_ */
    void diagonalize();
    /** compute rank, norm and determinant. */
    void compEstimates();
};

} // Namespace STK

#endif //STK_EIGENVALUESSYMMETRIC_H
