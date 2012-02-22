/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2011  Serge Iovleff

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
 * Project:  stkpp::Sdk
 * Purpose:  main interface base class for running class on 2D Container.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_IRunnerPtr2D.h
 *  @brief In this file we define the Interface base class for all the
 *  running classes on data set in 2D container given by ptr.
 **/

#ifndef STK_IRUNNERPTR2D_H
#define STK_IRUNNERPTR2D_H

#include "STK_IRunner.h"
#include "STK_ITContainer2D.h"

namespace STK
{

/** @ingroup Sdk
 *  @brief Abstract base class for all classes having a
 *  @code bool run(); @endcode method on a 2D Container data set.
 *  The data set will be passed by pointer.
 **/
template < class TYPE, class TContainer2D>
class IRunnerPtr2D : virtual public IRunnerBase
{
  /** Type of the container containing the data*/
  typedef ITContainer2D<TYPE, TContainer2D> Container2D;

  protected:
    /** default constructor */
    IRunnerPtr2D( Container2D const* p_data)
                : p_data_(p_data)
    { }
    /** copy constructor */
    IRunnerPtr2D( IRunnerPtr2D const& runner)
                : IRunnerBase(runner)
                , p_data_(runner.p_data_)
    { }

  public:
    /** destructor*/
    virtual ~IRunnerPtr2D() { }

    /** give the data set
     * @return a constant pointer on the data set.
     **/
    inline TContainer2D const* p_data() const
    { return (p_data_) ? p_data_->asPtrLeaf() : 0;}

    /** Set the data set and update the container.
     *  @param p_data The data set to run
     */
    virtual void setData( Container2D const* p_data)
    {
      p_data_ = p_data;
      update();
    }

    /** run the weighted computations.
     *  @param weights the weights of the samples
     *  @return @c true if no error occur during the running process, @c false
     *  otherwise
     **/
    virtual bool run( typename Container2D::TContainerVe const& weights) =0;

  protected:
    /** A pointer on the data set. */
    Container2D const* p_data_;

    /** update the runner.
     * This virtual method will be called when the state of the runner will
     * change.
     **/
    virtual void update() { ;}
};

} // namespace STK

#endif /* STK_IRUNNERPTR2D_H_ */
