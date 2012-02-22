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
 * Project:  stkpp::STKernel::Base
 * Purpose:  Define the Proxy classes for wrapping any types.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_Proxy.h
 *  @brief In this file we define the proxy classes for wrapping
 *  the STKpp fundamental types.
 **/

#ifndef STK_PROXY_H
#define STK_PROXY_H

namespace STK 
{

/** @ingroup Base
 *  The IProxy classes is an interface class for wrapping the fundamental
 * types.
 **/
 struct IProxy
 { /** Pure virtual destructor.*/
   virtual ~IProxy() {}
 };
 
/** @ingroup Base
 *  The ConstProxy class allow to surdefine operators and methods for
 *  every kind of constant class and type without using @c dynamic_cast.
 *  
 *  This class allow:
 *  - To avoid the use of @c dynamic_cast
 *  - To surdefine the predefined operators (like << and >>) for the
 *  C++ fundamental types.
 *  - To handle in a transparent way the NA values.
 **/
 template<class TYPE>
 class ConstProxy : public IProxy
 {
   protected:
     TYPE const& x_; ///< A constant reference on the object wrapped

   public:
     /** Default constructor : create a reference of the value x.
      *  @param x the value to wrap
      **/
     inline ConstProxy(TYPE const& x) : x_(x)
     { ;}

     /** Virtual destructor.
      **/
     virtual ~ConstProxy() {}

     /** Constant conversion operator to TYPE.
      **/
     inline operator TYPE const &() const
     { return x_;}
};

/** @ingroup Base
 *  The Proxy classe allow to surdefine operators and methods for
 *  every kind of class without using dynamic_cast.
 *  
 *  This class allow:
 *  - To avoid the use of @c dynamic_cast
 *  - To surdefine the predefined operators (like << and >>) for the
 *  C++ fundamental types.
 *  - To handle in a transparent way the NA values.
 **/
 template<class TYPE>
 class Proxy : public IProxy
 {
   protected:
     TYPE& x_; ///< A reference on the object wrapped

   public:
     /** Default constructor : create a reference of the Real x.
      *  @param x the object to wrap
      */
     inline Proxy(TYPE &x) : x_(x)
     { ;}

     /** overwrite with a TYPE.
      *  @param x the value to wrapp
      */
     inline Proxy& operator=(TYPE const& x)
     {
       x_ = x;
       return *this;
     }

     /** Virtual destructor.
      **/
     virtual ~Proxy() {}

     /** overwrite with a Proxy<TYPE>.
      *  @param x a wrapped value to copy
      **/
     inline Proxy& operator=(const Proxy<TYPE>& x)
     {
       x_ = x;
       return *this;
     }

     /** overwrite with a ConstProxy<TYPE>.
      *  @param x a constant wrapped value to set
      **/
     inline Proxy& operator=(const ConstProxy<TYPE>& x)
     {
       x_ = x;
       return *this;
     }

    /** Conversion operator to TYPE.
     **/
     inline operator TYPE &()
     { return x_;}
     
     /** Constant Conversion operator to TYPE.
      */
     inline operator TYPE const &() const
     { return x_;}
};

} // namespace STK

#endif /*STK_PROXY_H*/
