/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  Analysis
 * Purpose:  Implementation of functions around the beta function ratio
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Funct_betaRatio.cpp
 *  @brief In this file we implement functions around the Beta ratio
 *  function.
 **/

#include "../include/STK_Funct_betaRatio.h"
#include "../include/STK_Funct_gammaRatio.h"
#include "../include/STK_ISerie.h"
#include "../include/STK_Algo.h"
#include "../include/STK_Funct_gamma.h"
#include "../include/STK_Funct_util.h"
#include "../include/STK_Funct_raw.h"

#include "../../Arrays/include/STK_Vector.h"

namespace STK
{

namespace Funct
{

/** @ingroup Analysis
 *  @brief Compute the continued fraction of the beta function.
 *
 *  Compute the continued fraction:
 * \f[
 *    \frac{a_1}{b_1+}\frac{a_2}{b_2+}\frac{a_3}{b_3+}\ldots
 * \f]
 * with \f$ a_1 = 1 \f$, \f$ a_{n+1} c_n d_n /((a+2n-2)(a+2n))\f$
 * and \f$ b_{n+1} = (1-(a+b)x) + c_n + d_{n+1}\f$, where
 * \f[
 *   c_n = \frac{n(b-n)x}{a+2n-1} \qquad
 *   d_n = \frac{(a+n-1)(a+b+n-1)x}{a+2n-1}.
 * \f]
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 *  @param iterMax maximal number of iteration
 *  @return the value of the beta ratio function using the continued fraction
 *  method
 */
static Real betaRatio_cf( Real const& a
                        , Real const& b
                        , Real const& x
                        , bool lower_tail = true
                        , Integer const& iterMax = 1000
                        )
{
  // Constants
  Real sum = a+b, y = (0.5 -x) + 0.5, res;
  /* Compute the Function:
   *  B(a,b,x,y) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} x^{a} (1-x)^{b}
   */
  Real bt = ( Const::_1_SQRT2PI_*sqrt((a*b)/sum)
        * exp( ( gammaLnStirlingError(sum)
                -gammaLnStirlingError(a)
                -gammaLnStirlingError(b)
                )
              -(dev0(a, sum*x) + dev0(b, y*sum))
              )
        );
  if (bt == 0.) return lower_tail ? 0. : 1.;

  // initialize b_n
  Real bn = (1+a*y-b*x)/(a+1); // b_1
  // initialize numerator
  Real Nold = 0, Ncur = 1.; // = a_1 = 1
  // initialize denominator
  Real Dold = 1, Dcur = bn; // = b_1

  // normalize if necessary
  if (abs(Dcur) > 4294967296.0)
  {
    Ncur /= Dcur;
    Dold /= Dcur;
    Dcur = 1.;
  }

  // auxiliary constant variable
  const Real cte= (0.5-sum*x)+0.5;
  // dn (can be reused)
  Real dn = (a*sum*x)/(a+1);
  // result
  Real cf = 1/Dcur;
  // auxiliary variables
  Real ap2n = a+2, apn = a+1;
  // iterations
  for (Integer n=1; /*n<=iterMax*/; n++, apn++, ap2n += 2)
  {
    // compute auxiliaries variables
    Real nx = n*x, cn = n*(b*x-nx)/(ap2n-1);
    // compute a_n
    Real an = (cn/ap2n)*(dn/(ap2n-2));
    // check trivial case)
    if (an == 0.) break;
    // update d_n
    dn = apn*(sum*x+nx)/(ap2n+1);
    // compute bn
    bn = cte + cn + dn;
    // compute numerator and denominator
    Real anew = bn * Ncur + an * Nold;
    Nold = Ncur;
    Ncur = anew;
    anew = bn * Dcur + an * Dold;
    Dold = Dcur;
    Dcur = anew;
    // normalize if necessary
    if (abs(Dcur) > 4294967296.0)
    {
      Ncur /= Dcur;
      Nold /= Dcur;
      Dold /= Dcur;
      Dcur = 1.;
    }
    // normalize if necessary with 2^32
    if (abs(Ncur) > 4294967296.0)
    {
      Ncur /= 4294967296.0;
      Nold /= 4294967296.0;
      Dold /= 4294967296.0;
      Dcur /= 4294967296.0;
    }
    // test D_n not too small
    if (abs(Dcur) != 0)
    {
      Real cfold = cf;
      cf = Ncur/Dcur;
      // check cv
      if (abs(cf - cfold) < abs(cf)*Arithmetic<Real>::epsilon())
      {
        res = lower_tail ? bt*cf/a : 1-bt*cf/a;
        break;
      }
    }
  }
#ifdef STK_DEBUG
  // result
  res = lower_tail ? bt*cf/a : 1-bt*cf/a;
  if (!Arithmetic<Real>::isFinite(res))
  {
    std::cout << "in cf : bt= " << bt
              << " cf = " << cf
              << " a  = " << a
              << "lt =  " << lower_tail
              << "\n";
  }
  return res;
#else
  return lower_tail ? bt*cf/a : 1-bt*cf/a;
#endif

}

/* @ingroup Analysis
 *  @brief compute the coefs of the beta Ratio function in
 *  its asymptotic expansion.
 *
 *  Given the n-1 first coefs, compute and return the n-th
 *  coefficient. The value is stored back to the Vector. The
 *  formula for computing the coefficient is
 *  \[
 *     a_{n}
 *      = \frac{b-1}{(2n+1)!}+\frac{1}{n}
 *        \sum_{k=1}^{n-1} \left[
 *                          \frac{kb-n}{(2k+1)!} p_{n-k}
 *                         \right]
 *  \f]
 *  @param A vector of dimension 1:(n-1), n>=2, of the coefs
 *  @param b second parameter of the beta function
 */
//static Real coefs_ae(Real const& b, Vector &A)
//{
//  // get the number of exising coefs
//  Integer n = A.size();
//  // initialize the sum
//  Real sum = 0.0;
//  // update sum
//  for (Integer k=1; k<n; k++)
//  {
//    sum += (k*b/n-1)*A[n-k]/(factorial(2*k+1));
//  }
//  // compute result
//  sum += (b-1)/factorial(2*n+1);
//  // save it in A
//  A.push_back(sum);
//  // and return it
//  return sum;
//}

#define d1(z) (0.5)
#define d2(z) (z/8. - 1./20.)
#define d3(z) (z*(z/48.-1./40.)+1./105.)
#define d4(z) (z*(z*(z/384.-1./160.)+101./16800.)-3./1400.)
#define d5(z) (z*(z*(z*(z/3840.-1./960.)+61./33600.)-13./8400.) \
              +1./1925.)
#define d6(z) (z*(z*(z*(z*(z/46080.-1./7680.)+143./403200.) \
              -59./112000.)+7999./19404000.)-691./5255250.)
#define d7(z) (z*(z*(z*(z*(z*(z/645120.-1./76800.)+41./806400.) \
              -11./96000.)+5941./38808000.)-2357./21021000.)+6./175175.)
#define d8(z) (z*(z*(z*(z*(z*(z*(z/10321920.-1./921600.)+37./6451200.) \
              -73./4032000.)+224137./6209280000.)-449747./10090080000.) \
              +52037./1681680000.)-10851./1191190000.)


/** @ingroup Analysis
 *  @brief Compute the incomplete beta function ratio I_x(a,b)
 *  using its asymptotic expansion.
 *
 * Compute the incomplete beta function ratio I_x(a,b)
 * using it' asymptotic expansion.
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param xm1 true if we want to evaluate the function at 1-x
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 *  @return the value of the beta ratio function using its asymptotic expansion
 **/
static inline Real betaRatio_ae( Real const& a
                               , Real const& b
                               , Real const& x
                               , bool xm1
                               , bool lower_tail
                               )
{
  // Compute b-1
  Real bm1 = b-1;
  // compute \nu = a+(b-1)/2
  Real nu = a+bm1/2.;
  // compute D = -\nu*log(x)
  Real D = (xm1 ? -log1p(-x) : -log(x)) * nu;
  // check trivial case
  if (D == 0.) return lower_tail ? 1. : 0.;
  if (Arithmetic<Real>::isInfinite(D)) return lower_tail ? 0. : 1.;
  // compute 12*\nu^2
  nu *= 12*nu;
  // variables for the series
  Real term = 1.;
  Real ak_term = 1.;
  Real a_k = 1.;
  // numerator and denominator
  Real den = term;
  Real num = 0.;
  // update term
  term *= bm1;
  // Generalized bernouilli coefs
  // coef 1 (d_1 = 1/2)
  term *= (b)*(b+1)/nu;
  Real coef = term/2; // d_1 = 1/2
  den += coef;
  a_k += (ak_term *= D/(b+1));
  num += coef * a_k;
  // coef 2
  term *= (b+2)*(b+3)/nu;
  coef = d2(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+2));
  a_k += (ak_term *= D/(b+3));
  num += coef * a_k;
  // coef 3
  term *= (b+4)*(b+5)/nu;
  coef = d3(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+4));
  a_k += (ak_term *= D/(b+5));
  num += coef * a_k;
  // coef 4
  term *= (b+6)*(b+7)/nu;
  coef = d4(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+6));
  a_k += (ak_term *= D/(b+7));
  num += coef * a_k;
  // coef 5
  term *= (b+8)*(b+9)/nu;
  coef = d5(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+8));
  a_k += (ak_term *= D/(b+9));
  num += coef * a_k;
  // coef 6
  term *= (b+10)*(b+11)/nu;
  coef = d6(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+10));
  a_k += (ak_term *= D/(b+11));
  num += coef * a_k;
  // coef 7
  term *= (b+12)*(b+13)/nu;
  coef = d7(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+12));
  a_k += (ak_term *= D/(b+13));
  num += coef * a_k;
  // coef 8
  term *= (b+14)*(b+15)/nu;
  coef = d8(bm1)*term;
  den += coef;
  a_k += (ak_term *= D/(b+14));
  a_k += (ak_term *= D/(b+15));
  num += coef * a_k;
  // result (P or Q ?)
  return lower_tail ? gammaRatioQ(b, D)
                    + (num / den) * poisson_pdf_raw(b, D)
                    : gammaRatioP(b, D)
                    - (num / den) * poisson_pdf_raw(b, D);
}

/** @ingroup Analysis
 *  @brief Compute the recurrence formula of the incomplete beta
 *  ratio function.
 *
 *  Compute the recurrence formula of the incomplete beta
 *  ratio function I_x(a,b).
 *
 *  @param s a+b
 *  @param a first parameter
 *  @param x value to evaluate the function
 *  @param n number of step
 **/
static inline Real serie_up( Real const& s
                           , Real const& a
                           , Real const& x
                           , Integer const& n
                           )
{
  // initilize d0 and sum
  Real sum = 1./a;
  Real di  = 1./a;
  // compute \sum_{i=0}^{n-1} \frac{d_i x^i}{a}
  for (Integer i=1; i<n; i++)
    sum += (di *= x*(s-1+i)/(a+i));
  // return result
  return sum;
}

/** @ingroup Analysis
 *  @brief Compute the incomplete beta function ratio I_x(a,b)
 *  using its recurrence formula and its asymptotic expansion.
 *
 * Compute the incomplete beta function ratio I_x(a,b)
 * using its recurrence formula and its asymptotic expansion.
 *
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param xm1 true if we want to evaluate the function at 1-x
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 **/
static inline Real betaRatio_up( Real const& a
                               , Real const& b
                               , Real const& x
                               , bool xm1
                               , bool lower_tail
                               )
{
  // number of iterations
  Integer n = Integer(a);
  // compute residual
  Real a0 = a-n;
  // zero case
  if (!a0) { --n; a0=1.;}
  // sum
  Real s= a0+b;
  // compute x and 1-x
  Real x0, y;
  if (xm1)
  { x0 = (0.5 -x) + 0.5; y =x;}
  else
  { x0 = x; y = (0.5 -x) + 0.5;}
  /* Compute the constant
   *  B(a,b,x,y) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} x^a (1-x)^b
   * using saddle-point expansion
   */
  Real bt = ( Const::_1_SQRT2PI_*sqrt((a0*b)/s)
            * exp( ( gammaLnStirlingError(s)
                    -gammaLnStirlingError(a0)
                    -gammaLnStirlingError(b)
                   )
                   -(dev0(a0, s*x0) + dev0(b, y*s))
                 )
            );
  // check trivial case
  if (!bt)
    return betaRatio_ae(b, a0, x, !xm1, !lower_tail);
  // compute serie
    bt *= serie_up(s, a0, x0, n);
  // general case
  return  lower_tail ? betaRatio_ae(b, a0, x, !xm1, !lower_tail) - bt
                     : betaRatio_ae(b, a0, x, !xm1, !lower_tail) + bt;
}

/** @ingroup Analysis
 *  @brief Compute the incomplete beta function ratio I_x(a,b)
 *  using its series representation.
 *
 * Compute the incomplete beta function ratio I_x(a,b)
 * using it's serie representation.
 * \f[
 * \frac{x^a}{B(a,b)}
 * \left(  \frac{1}{a}
 *       + \sum_{n=1}^{\infty} \frac{(1-b)(2-b)\ldots(n-b)}{n! (a+n)}
 * \right)
 * \f]
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 *  @return the value of the beta ratio function using its series representation
 **/
static Real betaRatio_sr( Real const& a
                        , Real const& b
                        , Real const& x
                        , bool lower_tail
                        )
{
  // Constant
  Real sum = a+b;
  // compute B(a,b,x,y) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} x^{a}
  Real bt = ( Const::_1_SQRT2PI_*sqrt((a*b)/sum)
            * exp( sum*x
                 + ( gammaLnStirlingError(sum)
                    -(gammaLnStirlingError(a)+gammaLnStirlingError(b))
                   )
                 -(dev0(a, sum*x) + dev0(b, sum))
                 )
            );
  // check factor
  if (!bt) return lower_tail ? 0. : 1.;
  // parameters
  Real term, n= 1., c = 1.;
  // sum = \sum_{n=1}^\infty (1-b/n) x^n /(a+n)
  sum = 0.0;
  do
  {
    sum += (term =  (c *= (1-b/n)*x) / (a+n));
    ++n;
  }
  while (abs(term) > abs(sum) * Arithmetic<Real>::epsilon());
  // return result
  return lower_tail ? bt*(1/a + sum)
                    : (1-bt/a)-(bt*sum);
}

/** @ingroup Analysis
 *  @brief compute the odd coefs of the beta Ratio function
 *  serie expansion.
 *
 *  Given the 2t first coefs, compute and return the (2t+1)-th
 *  coefficient. The value is stored to the back of the Vector.
 *
 *  @param A vector of dimension 1:n, n>=2, of the coefs
 *  @param std the binomial standard deviation \f$ \sqrt{pq} \f$
 *  @param qmp the difference q-p
 */
static Real coefs_odd_se( Real const& std, Real const& qmp, Vector &A)
{
  // get the number of exising coefs (should be even)
  Integer n = A.size();
  Integer t = n/2;
  // know case
  if (t==1)
  {
    // compute result
    Real sum = (qmp/3+std)*(qmp/3-std)/(4*std);
    // save it in A
    A.push_back(sum);
    // and return it
    return sum;
  }
  // initialize the sums
  Real sum, sum1 =0, sum2 =0;
  // update sum
  for (Integer k=1, k1=2; k<t; k++, k1++)
  {
    sum1 += A[k1]*A[n+1-k];
    sum2 += A[k] *A[n-k];
  }
  // compute the result
  sum = ( ((qmp*A[n] - A[t]*A[t])/2 - sum2)/(t+1)
          - sum1- A[t+1]*A[t+1]/2
        )/std;
  // save it in A
  A.push_back(sum);
  // and return it
  return sum;
}

/** @ingroup Analysis
 *  @brief compute the even coefs of the beta Ratio function
 *  serie expansion.
 *
 *  Given the 2t+1 first coefs, compute and return the (2t+2)-th
 *  coefficient. The value is stored to the back of the Vector.
 *
 *  @param A vector of dimension 1:n, n>=3, of the coefs
 *  @param std the binomial standard deviation \f$ \sqrt{pq} \f$
 *  @param qmp the difference q-p
 */
static Real coefs_even_se( Real const& std, Real const& qmp, Vector &A)
{
  // get the number of exising coefs (should be odd)
  Integer n = A.size();
  Integer t = n/2;
  // know case
  if (t==1)
  {
    // compute result
    Real sum = -qmp*(1+A[3]/std)*2/15;
    // save it in A
    A.push_back(sum);
    // and return it
    return sum;
  }
  // initialize the sums
  Real sum, sum1 =0, sum2 =0;
  // update sum
  for (Integer k=1, k1=2; k<=t; k++, k1++)
  {
    sum1 += A[k1]*A[n+1-k];
    sum2 += A[k] *A[n-k];
  }
  // compute the result
  sum = ( (qmp*A[n] - 2*sum2)/(n+2) - sum1)/std;
  // save it in A
  A.push_back(sum);
  // and return it
  return sum;
}

/** @ingroup Analysis
 *  @brief Compute the incomplete beta function ratio I_x(a,b)
 *  using its serie expansion.
 *
 * Compute the incomplete beta function ratio I_x(a,b)
 * using it's serie expansion.
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param xm1 @c true if x is in [q,p], @c false otherwise
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 *  @param iterMax Maximal number of iteration (default is 20)
 **/
static Real betaRatio_se( Real const& a
                        , Real const& b
                        , Real const& x
                        , bool xm1
                        , bool lower_tail
                        , Integer const& iterMax = 20
                        )
{
  // parameters
  Real s  = a+b;
  Real p  = a/s, q  = b/s;
  Real sx = s*x;
  Real sy = s-sx;
  //Real D  = a*log(x/p)+b*log(y/q);
  Real D  = (xm1 ? dev0(a, sy)+dev0(b, sx) : dev0(a, sx)+dev0(b, sy));

  Real z2 = 2 * D;
  Real z  = (((x<p)&&!xm1)||((x>q)&&xm1)) ? -sqrt(z2) : sqrt(z2);

  //  Compute normal cdf
  Real pnorm = lower_tail ?  normal_cdf_raw( z) : normal_cdf_raw(-z);
#ifdef STK_DEBUG
  if (!Arithmetic<Real>::isFinite(pnorm))
  {
    stk_cout << _T(" a= ") << a
             << _T(" b= ") << b
             << _T(" x= ") << x
             << _T(" p= ") << p
             << _T(" q= ") << q
             << _T(" D= ") << D
             << _T("\n");
  }
#endif
  //  Compute normal pdf
  Real dnorm = Const::_1_SQRT2PI_ * exp(-D);
  // check large values of z
  if (!dnorm) return pnorm;
  // epsilon value for the serie
  Real eps = Arithmetic<Real>::epsilon()*pnorm;

  // auxiliary variables
  Real std = sqrt(a)*sqrt(b)/s, qmp = q - p;

  // serie coefs
  Vector A(2);
  A.reserve(iterMax); // reserve enough space
  A[1] = std;
  A[2] = qmp/3;
  // variables for the series
  Real odd_term  = 1;
  Real even_term = 2/sqrt(s);

  Real sum1_term = z, sum1 = z;
  Real sum2_term = 1, sum2 = 1;

  // variables for the numerator and the denominator
  Real num = qmp*even_term/3; // = A[2]*even_term;
  Real den = std;             // = A[1];

  // compute ratio
  Real res = num/den;

  // computation
  for (Integer k=1; k<iterMax; k++)
  {
    // odd_term : (2k+1)!!/(a+b)^k
    odd_term *= (2*k+1)/s;
    // even_term (2*k)!!/(a+b)^{k+1/2}
    even_term *= 2*(k+1)/s;

    // compute a_{2k+1}
    Real a_n = coefs_odd_se(std, qmp, A);
    // check result
    if (Arithmetic<Real>::isFinite(a_n))
    {
      // denominator \sum_0^k (2j+1)!!a_{2j+1}/(a+b)^j
      den += odd_term * a_n;
      // num+=\sum_1^k (2k+1)!!a_{2k+1}/(a+b)^k (\sum_1^k z^{2j-1}/(2j-1)!!)
      num += odd_term * a_n * sum1;

      // compute a_{2k+2}
      a_n = coefs_even_se(std, qmp, A);
      // check result
      if (Arithmetic<Real>::isFinite(a_n))
      {
        // update sum2_term = z^{2k}/(2k)!!
        // update sum2      = \sum_0^k z^{2j}/(2j)!!
        sum2 += (sum2_term *= z2/(2*k));
        // num+=\sum_0^k (2k+2)!!a_{2k+2}/(a+b)^k (\sum_1^k z^{2j}/(2j)!!)
        num += even_term * a_n * sum2;
      }
    }
    // test cv
    Real old_res = res;
    res = num/den;
    if (abs(old_res - res)*dnorm < eps)
    { break;}

    // update sum1_term = z^{2k+1}/(2k+1)!!
    // update sum1      = \sum_1^{k+1} z^{2j-1}/(2j-1)!!
    sum1 += (sum1_term *= z2/(2*k+1));
  }
#ifdef STK_DEBUG
  // result
  Real result = lower_tail ? pnorm - res * dnorm
                           : pnorm + res * dnorm;
  if (!Arithmetic<Real>::isFinite(result))
  {
    std::cout << "Error in betRatio_se\n"
              << "pnorm = " << pnorm
              << " res = " << res
              << " dnorm  = " << dnorm
              << "lt =  " << lower_tail
              << "\n";
  }
  return result;
#else
  // result
    return lower_tail ? pnorm - res * dnorm
                      : pnorm + res * dnorm;
#endif
}

/** @ingroup Analysis
 *  Compute the beta ratio function.
 *  \f[
 *     I_x(a,b) = \int_0^x u^{a-1} (1-u)^{b-1} du
 *  \f]
 *  for \f$ 0\leq x \leq 1\f$.
 *
 *  @param a first parameter, must be >0
 *  @param b second parameter, must be >0
 *  @param x value to evaluate the function
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 **/
Real betaRatio( Real const& a, Real const& b, Real const& x
              , bool lower_tail)
{
  // Check if a and x are available
  if (  Arithmetic<Real>::isNA(a)
     || Arithmetic<Real>::isNA(b)
     || Arithmetic<Real>::isNA(x)
     ) return Arithmetic<Real>::NA();
  // Negative parameter not allowed
  if ((a<=0)||(b<=0))
    throw domain_error("Funct::betaRatio(a,b,x) "
                            "Negative parameter a or b");
  // trivial case
  if (x<=0) return lower_tail ? 0. : 1.;
  if (x>=1) return lower_tail ? 1. : 0.;
  // compute 1-x
  Real y = (0.5 - x) + 0.5;
  // parameters
  Real p=a/(a+b);
  Real q=b/(a+b);
  // small a or b
  if (min(a,b)<=1)
  {
    // case b<=1 and a>=15
    if (a>=15)
      return betaRatio_ae(a, b, x, false, lower_tail);
    // case a<=1 and b>=15
    if (b>=15)
      return betaRatio_ae(b, a, x, true, !lower_tail);
    // ((a<15) and (b<15) and (min(a,b)<=1))
    return (x<=0.5) ? betaRatio_sr(a, b, x, lower_tail)
                    : betaRatio_sr(b, a, y, !lower_tail);
  }
  // (1<a<=15) and (1<b<=15)
  if (max(a,b)<=15)
  {
    // general case
    return (x < p) ? betaRatio_cf(a, b, x,  lower_tail)
                   : betaRatio_cf(b, a, y, !lower_tail);
  }
  // ((1<a<=15) and (15<b)) or ((15<a) and (1<b<=15))
  if (min(a,b)<=15)
  {
    return (a<=15) ? betaRatio_up(a, b, x, false, lower_tail)
                   : betaRatio_up(b, a, x, true, !lower_tail);
  }

//  Real x0, y0, a0, b0, p0, q0;
//  bool flag;
//  if (x>p)
//  { x0=y; y0=x; a0=b; b0=a; flag = !lower_tail; p0=q; q0=p;}
//  else
//  { x0=x; y0=y; a0=a; b0=b; p0=p; q0=q; flag = lower_tail;}
//  // center region for x
//  if ((x0>1.03*p0)&&(y0<0.97*q0))
//    return betaRatio_se(a0-1, b0-1, x0, false, flag);
//  return betaRatio_cf(a0, b0, x0, flag);
  // (a>15) and (b>15) and (a<b)
  if (a<b)
  {
    // center region for x in [p , q]
    if ((x>=p-0.01*(q-p))&&(y<=q+0.01*(q-p)))
      return betaRatio_se(a-1, b-1, x, false, lower_tail);
    // extreme values
    return (x<p) ? betaRatio_cf(a, b, x,  lower_tail)
                 : betaRatio_cf(b, a, y, !lower_tail);
  }
  // (a>15) and (b>15) and (b<a)
  // center region for x in [q, p]
  if ((x>=q-0.01*(p-q))&&(y<=p+0.01*(p-q)))
    return betaRatio_se(b-1, a-1, x, true, !lower_tail);
  // extreme values
  return (x<p) ? betaRatio_cf(a, b, x,  lower_tail)
               : betaRatio_cf(b, a, y, !lower_tail);
}

} // namespace Funct

} // namespace STK
