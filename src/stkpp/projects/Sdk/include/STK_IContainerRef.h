/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004  Serge Iovleff

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
 * Project: stkpp::Sdk
 * Purpose:  Define the Reference Container interface class.
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_IContainerRef.h
 *  @brief This is an internal header file, included by other
 *  Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class, except if you want to create your own Container
 *  Class.
 **/

#ifndef STK_ICONTAINERREF_H
#define STK_ICONTAINERREF_H

namespace STK
{
/** @ingroup TContainer
 *  @brief Base class for all referencing containers.
 *
 *  The IContainerRef class is the base class for all containers that
 *  can be referenced. If a container is derived from this class, then
 *  this container can be a reference (a wrapper) on an other container.
 *  A container R is a reference of the container A, if it wrap the
 *  data contained in A. In this case, the boolean @c ref_ is true.
 **/
class IContainerRef
{
  private:
    /** Is it a "true" container or a wrapper ?
     *  ref_ should be @c false if this own its own data, @c true otherwise.
     **/
    bool ref_;

  protected:
    /** Default constructor We have to specify the member ref_.
     *  @param ref : false if this own its own data.
     **/
    inline IContainerRef( bool ref) : ref_(ref)
    { ;}

    /** Copy constructor
     * @param T : The container to copy.
     * @param ref : is this a wrapper of T ?
     **/
    inline IContainerRef( const IContainerRef& T, bool ref) : ref_(ref)
    { ;}

  public:
    /** Virtual destructor.                                                 */
    inline virtual ~IContainerRef() { ;}

    /** is this own its own data ?
     *  @return @c true if this is reference container, @c false otherwise
     **/
    inline bool isRef() const { return ref_;}

  protected:
    /** Modify the container : can become a reference or the owner of
     *  the data. To use with care if we want to avoid memory leak.
     *  @param ref : false if this own its own data.
     **/
    inline void setRef(bool ref)
    { ref_ = ref;}
};

} // namespace STK

#endif
// STK_CONTAINERREF_H
