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
 * Project:  stkpp::STatistiK::Law
 * Purpose:  Interface base class for all univariate probabilities laws.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Law_ITUnivariate.h
 *  @brief In this file we define the interface base class ITUnivariate for all
 *  probabilities laws.
 **/

#ifndef STK_LAW_ITUNIVARIATE_H
#define STK_LAW_ITUNIVARIATE_H

// RandBase header
#include "STK_Law_ILawBase.h"
#include "../../Sdk/include/STK_ITContainer2D.h"

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Interface base class for all the univariate probabilities laws.
 *
 *  A general probability law (discrete or real) possess a
 *  probability density law. It can be simulated.
 *
 * Interface base class for the unidimensional probabilities laws. Every derived
 * class have to furnish :
 * - a random generator method
 * @code TYPE rand() const @endcode
 * - a quantile method (inverse cdf)
 * @code TYPE icdf(Real const& p) const @endcode .
 *
 * The derived objects should also furnish the same static functions.
 * Instantiation of a derived object is interesting when one want to
 * simulate independent identical distributed random variates : the
 * creation of the object initialize all parameter-dependent variables.
 **/
template <class TYPE>
class ITUnivariate : public ILawBase
{
  protected:
    /** Constructor.
     *  @param name the name of the law
     **/
    ITUnivariate(String const& name) : ILawBase(name)
    { ;}

  public:
    /** Virtual destructor. **/
    virtual ~ITUnivariate() { ;}

    /** Pseudo-random Real law generator for a one dimensional
     *  container of @c TYPE.
     *  @param A the container to store the random numbers
     **/
    template< class Container1D>
    void rand1D( ITContainer1D< TYPE, Container1D>& A) const
    {
      // get dimensions
      const Integer first = A.first(), last = A.last();
      // generate and set random variables
      for (Integer i=first; i<=last; i++) A[i] = rand();
    }

    /** Pseudo-random Real law generator for a two dimensional
     *  container of Real.
     *  @param A the container to store the random numbers
     **/
    template < class TContainer2D>
    void rand2D( ITContainer2D< TYPE, TContainer2D>& A) const
    {
      // get dimensions
      const Integer firstRow = A.firstRow(), lastRow = A.lastRow();
      const Integer firstCol = A.firstCol(), lastCol = A.lastCol();
      // generate and set random variables
      for (Integer j=firstCol; j<=lastCol; j++)
        for (Integer i=firstRow; i<=lastRow; i++)
          A(i, j) = rand();
    }

  private:
    /** @brief compute the cumulative distribution function
     *  Give the probability that a random variate is less or equal
     *  to t.
     *  @param t the value to compute the cdf.
     *  @return the value of the cdf
     **/
    virtual Real cdf(Real const& t) const =0;

    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x the value to compute the pdf.
     *  @return the value of the pdf
     **/
    virtual Real pdf(TYPE const& x) const =0;
    
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x the value to compute the lpdf.
     *  @return the value of the log-pdf
     **/
    virtual Real lpdf(TYPE const& x) const = 0;
    /**
     *  Generate a @c TYPE random variate .
     **/
    virtual TYPE rand() const =0;

    /** @brief inverse cumulative distribution function
     *  Compute the Real quantile t such that the probability of a random
     *  variate less to t is less or equal to p.
     *  @param p value of the probability giving the quantile
     **/
    virtual Real icdf(Real const& p) const=0;
};

} // namespace Law

} //  namespace STK

#endif /*STK_LAW_ITUNIVARIATE_H*/
