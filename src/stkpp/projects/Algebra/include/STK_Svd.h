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
 * Project:  Algebra
 * Purpose:  Define The Svd Class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Svd.h
 *  @brief In this file we define the Svd Class.
 **/
 
#ifndef STK_SVD_H
#define STK_SVD_H

#include "../../Arrays/include/STK_Matrix.h"
#include "../../Arrays/include/STK_MatrixSquare.h"

namespace STK
{

/** @ingroup Algebra
 *  @brief The class Svd compute the Singular Value Decomposition
 *  of a Matrix with the Golub-Reinsch Algorithm.
 * 
 *  The method take as:
 *  - input: A matrix A(nrow,ncol)
 *  - output:
 *    -# U Matrix (nrow,ncolU).
 *    -# D Vector (ncol)
 *    -# V Matrix (ncol,ncol).
 *  and perform the decomposition: 
 *  - A = UDV' (transpose V).
 *  U can have more cols than A,
 *  and it is possible to ompute some (all) vectors of Ker(A).
 **/
class Svd
{
  protected:
    /* containers */
    Matrix       U_;    ///< U_ matrix
    MatrixSquare V_;    ///< V_ square matrix
    Point        D_;    ///< Array of the singular values

    /* dimensions */
    Integer ncolV_;   ///< Number of cols (and rows) of V
    Integer ncolD_;   ///< Number of cols of D_
    Integer ncolU_;   ///< Number of cols of U
    Integer nrowU_;   ///< Number of rows of U_

    /* flags */
    bool withU_;        ///< Compute U_ ?
    bool withV_;        ///< Compute V_ ?
    bool ref_;          ///< Is this structure just a pointer on U_ ?

    /* results */
    Real    norm_;       ///< norm of the matrix (largest singular value)
    Integer rank_;       ///< rank of the matrix
    bool    error_;      ///< Everything OK during computation ?

  public :   
    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if true, we save the left housolder transforms
     *  in U_.
     *  @param withV if true, we save the right housolder transforms
     *  in V_.
     **/
    Svd( Matrix const&    A       = Matrix()
       , bool      ref     = false
       , bool      withU   = true
       , bool      withV   = true
       );

    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    Svd( const Svd &S);

    /** destructor.
     **/
    virtual ~Svd();
 
    /** Operator = : overwrite the Svd with S.
     *  @param S the Svd to copy
     **/
    Svd& operator=(const Svd &S);

    /** clear U_.                                                     */
    void clearU();

    /** clear V_.                                                     */
    void clearV();

    /** clear U_, V_ and D_.                                          */
    void clear();

    /** Compute the svd of the Matrix A and copy the data
     *  see the corresponding constructor Take care that if U_ was previously
     *  a reference, it cannot be modified.
     *  @param A is the matrix to decompose.
     *  @param withU if true, we save the left housolder transforms
     *  in U_.
     *  @param withV if true, we save the right housolder transforms
     *  in V_.
     **/    
    void newSvd( Matrix const&    A       = Matrix()
               , bool      withU   = true
               , bool      withV   = true
               );

    /** Computing the bidiagonalisation of M.
     *  The diagonal and the subdiagonal are stored in D and F
     *  @param M the matrix to bidiagonalize, the matrix is overwritten
     *  with the left and right Householder vectors.
     *  The method return the estimate of the inf norm of M.
     *  @param D the element of the diagonal
     *  @param F the element of the surdiagnal
     **/
    static Real bidiag(const Matrix& M, Point& D, Vector& F);

    /** right eliminate the element on the subdiagonal of the row nrow
     *  @param D the diagonal of the matrix
     *  @param F the surdiagonal of the matrix
     *  @param nrow the number of the row were we want to rightEliminate
     *  @param V a right orthogonal Square Matrix
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static void rightEliminate( Point& D
                              , Vector& F
                              , Integer const& nrow
                              , MatrixSquare& V
                              , bool withV = true
                              , Real const& tol = Arithmetic<Real>::epsilon()
                              );

    /** left eliminate the element on the subdiagonal of the row nrow
     *  @param D the diagonal of the matrix
     *  @param F the surdiagonal of the matrix
     *  @param nrow the number of the row were we want to rightEliminate
     *  @param U a left orthogonal Matrix
     *  @param withU true if we want to update U
     *  @param tol the tolerance to use
     **/
    static void leftEliminate( Point& D
                             , Vector& F
                             , Integer const& nrow
                             , Matrix& U
                             , bool withU = true
                             , Real const& tol = Arithmetic<Real>::epsilon()
                             );

    /** Computing the diagonalisation of a bidiagnal matrix
     *  @param D the diagoanl of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param U a left orthogonal Matrix
     *  @param withU true if we want to update U
     *  @param V a right orthogonal Square Matrix
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static bool diag( Point& D
                    , Vector& F
                    , Matrix& U
                    , MatrixSquare& V
                    , bool withU = true
                    , bool withV = true
                    , Real const& tol = Arithmetic<Real>::epsilon()
                    );

    /// Number of rows of U_
    inline Integer nrowU() const { return U_.sizeVe();}
    /// Number of cols of U_
    inline Integer ncolU() const { return U_.sizeHo();}
    /// Number of rows of D_
    inline Integer ncolD() const { return ncolD_;}
    /// Number of rows of V_
    inline Integer nrowV() const { return V_.sizeVe();}
    /// Number of cols of V_
    inline Integer ncolV() const { return V_.sizeHo();}
    /// Norm of the matrix
    inline Real normSup()  const { return norm_;}
    /// rank of the matrix
    inline Integer rank()  const { return rank_;}
    /// if an error occur during svd()
    inline bool error()    const { return error_;}
    /// get U (const)
    inline Matrix const&       getU() const { return U_;}
    /// get V (const)
    inline MatrixSquare const& getV() const { return V_;}
    /// get D (const)
    inline const Point&        getD() const { return D_;}
    
    /// get U
    inline Matrix&       getU() { return U_;}
    /// get V
    inline MatrixSquare& getV() { return V_;}
    /// get D
    inline Point&        getD() { return D_;}
    
  private:
    /// Values of the Surdiagonal
    Vector F_;
    /// Initialize the containers
    void init();
    /// Svd main steps
    void compSvd();
    /// Compute U (if withU_ is true)
    void compU();
    /// Compute V (if withV_ is true)
    void compV();
};

} // namespace STK

#endif
// STK_SVD_H
