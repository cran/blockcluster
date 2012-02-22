/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2010  Serge Iovleff

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
 * Project:  stkpp::Regress
 * created on: 25 juin 2010
 * Purpose:  Compute the coefficient of a B-Spline curves.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_BSplineCoefficients.h
 *  @brief In this file we define the BSplineCoefficients class.
 **/

#ifndef STK_BSPLINECOEFFICIENTS_H
#define STK_BSPLINECOEFFICIENTS_H

#include "../../Arrays/include/STK_Matrix.h"

namespace STK
{


/** @brief Compute the regression splines coefficients.
 * The BSplineCoefficients class compute the coefficients of a B-Spline curve
 * using the de Boor's algorithm. The knots can be uniform (the default) or
 * periodic.
 *
 * The input data set is a vector of size n and the output matrix of the
 * coefficients @c Coefficients() is a matrix of size (n, nbControlPoints).
 */
class BSplineCoefficients
{
  public:
    /** Method to use for positioning the knots. */
    enum KnotsPosition
    {
      uniform_  ///< uniform knots
    , periodic_ ///< periodic knots
    , density_  ///< knots using density of the data
    , unknown_  ///< unknown method
    };

    /** convert a String to a TypeGraph.
     *  @param type the type of graph in a string
     *  @return the TypeGraph represented by the String @c type. If the string
     *  does not match any known name, the @c unknown_ type is returned.
     **/
    static KnotsPosition StringToKnotsPosition( String const& type);

    /** convert a KnotsPosition to a String.
     *  @param type the KnotsPosition we want to convert to a string
     *  @return the string associated to this KnotsPosition
     **/
    static String KnotsPositionToString( KnotsPosition const& type);

    /**
     *  Constructor : initialize the data members. The number of knots is given
     *  by the formula nbKnots = nbControlPoints + degree +1.
     *  @param p_data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-Spline curves
     *  @param position method to use for positioning the knots
     **/
    BSplineCoefficients( Vector const* p_data
                       , Integer const& nbControlPoints
                       , Integer const& degree = 3
                       , KnotsPosition const& position = uniform_
                       );

    /**
     *  Constructor : initialize the data members. The number of knots is given
     *  by the formula nbKnots = nbControlPoints + degree +1.
     *  @param data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-Spline curves
     *  @param position method to use for positioning the knots
     **/
    BSplineCoefficients( Vector const& data
                       , Integer const& nbControlPoints
                       , Integer const& degree = 3
                       , KnotsPosition const& position = uniform_
                       );

    /** Destructor. */
    virtual ~BSplineCoefficients();

    /** run the computations. */
    void run();

    /** Compute the coefficients of the B-Spline curve for the given values.
     *  @param p_data the input data values
     *  @param nbControlPoints number of control points
     *  @param degree degree of the B-Spline curves
     *  @param position method to use for positioning the knots
     **/
    void setData( Vector const* p_data
                , Integer const& nbControlPoints
                , Integer const& degree = 3
                , KnotsPosition const& position = uniform_
                );

    /** give the degree of the B-Spline curves. */
    inline Integer const& degree() const { return degree_;}
    /** give the number of knots of the B-Spline curves. */
    inline Integer const& nbKnots() const { return nbKnots_;}
    /** give the number of control points of the B-Spline curves. */
    inline Integer const& nbControlPoints() const { return nbControlPoints_;}
    /** give the knots of the B-Spline curves. */
    inline Vector const& knots() const { return knots_;}
    /** The computed coefficients of the B-Spline curves.
     *  This is a matrix of size (p_data_->range(), 0:lastControlPoints).
     **/
    inline Matrix const& coefficients() const { return coefficients_;}

  protected:
    /** the input data set */
    Vector const* p_data_;
    /** number of knots of the B-Spline curves.*/
    Integer nbKnots_;
    /** Index of the last knots of the B-Spline curves. This is nbKnots_ - 1.*/
    Integer lastKnot_;
    /** number of control points of the B-Spline curves.*/
    Integer nbControlPoints_;
    /** Index of the last control point of the B-Spline curves.
     *  This is nbControlPoints_ - 1.
     **/
    Integer lastControlPoint_;
    /** degree of the B-splines curves. */
    Integer degree_;
    /** Method used in order to position the knots. */
    KnotsPosition position_;
    /** Vector of the knots */
    Vector knots_;
    /** Matrix of the coefficients */
    Matrix coefficients_;

  private:
    /** Minimal value of the knots */
    Real minValue_;
    /** Maximal value of the knots */
    Real maxValue_;
    /** compute the position of the knots of the B-Spline curves.*/
    void computeKnots();
    /** compute the position of the uniform knots.*/
    void computeUniformKnots();
    /** compute the position of the periodic knots.*/
    void computePeriodicKnots();
    /** compute the position of the density knots.*/
    void computeDensityKnots();
    /** Compute the coefficients of the B-Spline curves.*/
    void computeCoefficients();
    /** Compute a row of the coefficients
     * @param irow index of the row
     **/
    void computeCoefficientsRow(Integer const& irow);
};

}

#endif /* STK_BSPLINECOEFFICIENTS_H */
