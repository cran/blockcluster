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
 * created on: 29 juil. 2011
 * Purpose:  main interface base class for running method.
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_IRunner.h
 *  @brief In this file we define the Interface base class for all the
 *  running classes.
 **/

#ifndef STK_IRUNNER_H
#define STK_IRUNNER_H

#include "../../STKernel/include/STK_String.h"

namespace STK
{
/** @ingroup Sdk
 *  @brief Abstract base class for all classes having a
 *  @code bool run(); @endcode method.
 *  All classes deriving from this class should implement the @c run method
 *  using this kind of code
 *  @code
 *  bool MyClass::run()
 *  {
 *    try
 *    {
 *    // do something
 *
 *    }
 *    catch(const Exception& e)
 *    { msg_error_ = e.error();
 *      return false;
 *    }
 *    return true;
 *  }
 *  @endcode
 **/
class IRunnerBase
{
  protected:
    /** default constructor */
    IRunnerBase();

    /** copy constructor
     * @param runner the runner to copy
     **/
    IRunnerBase( IRunnerBase const& runner);

  public:
    /** destructor*/
    virtual ~IRunnerBase();

    /** get the last error message.
     * @return the last error message
     **/
    inline String const& error() const { return msg_error_;}

    /** run the computations.
     * @return @c true if no error occur during the running process, @c false
     * otherwise
     **/
    virtual bool run() =0;

  protected:
    /** String with the last error message. */
    String msg_error_;
};

/** @ingroup Sdk
 *  @brief Abstract class for all classes making unsupervised learning.
 *
 *  This Interface is designed for unsupervised learning purpose. In a
 *  supervised learning setting, use IRunnerConstRefReg.
 *  The data set y is not copied. There is just a pointer on the data set
 *  stored internally.
 *
 *  The pure virtual methods to implement are
 *  @code
 *    bool run();
 *    bool run(weights);
 *  @endcode
 **/
template < typename TY>
class IRunnerConstRef : virtual public IRunnerBase
{
  protected:
    /** default constructor
     *  @param y The data set to run
     **/
    IRunnerConstRef( TY const& y)
                   : p_y_(&y)
    { }
    /** copy constructor */
    IRunnerConstRef( IRunnerConstRef const& runner)
                : IRunnerBase(runner)
                , p_y_(runner.p_y_)
    { }

  public:
    /** destructor*/
    virtual ~IRunnerConstRef() { }

    /** get the data set
     * @return a constant reference on the data set.
     **/
    inline TY const& y() const { return *p_y_;}

    /** set the data set. If the state of the runner change when a new data
     *  set is set, the user of this class have to overload the udpateY(
     *  method.
     *  @param y The data set to run
     */
    virtual void setY( TY const& y)
    {
      p_y_ = &y;
      update();
    }

    /** run the weighted computations.
     *  @return @c true if no error occur during the running process, @c false
     *  otherwise
     **/
    virtual bool run() =0;
    /** run the weighted computations.
     *  @param weights the weights of the samples
     *  @return @c true if no error occur during the running process, @c false
     *  otherwise
     **/
    virtual bool run( typename TY::TContainerVe const& weights) =0;

  protected:
    /** A pointer on the original data set. */
    TY const* p_y_;
    /** update the runner.
     * This virtual method will be called when the state of the runner will
     * change, i.e. when a new data set is set. By default do nothing.
     **/
    virtual void update() {}
};

/** @ingroup Sdk
 *  @brief Abstract class for all classes making supervised learning.
 *
 *  This Interface is designed for supervised learning purpose. In an
 *  unsupervised learning setting, use IRunnerConstRef.
 *  The data sets x and y are not copied. There is just pointers on the data
 *  sets stored internally.
 *
 *  The pure virtual method to implement are
 *  @code
 *    bool run();
 *    bool run(weights);
 *  @endcode
 *
 **/
template < typename TY, typename TX>
class IRunnerConstRefReg : virtual public IRunnerBase
{
  protected:
    /** default constructor
     *  @param x The x data set to run
     *  @param y The y data set to run
     **/
    IRunnerConstRefReg( TY const& y, TX const& x)
                      : p_y_(&y)
                      , p_x_(&x)
    { }
    /** copy constructor
     * @param runner the runner to copy
     **/
    IRunnerConstRefReg( IRunnerConstRefReg const& runner)
                      : IRunnerBase(runner)
                      , p_y_(runner.p_y_)
                      , p_x_(runner.p_x_)
    { }

  public:
    /** destructor*/
    virtual ~IRunnerConstRefReg() { }

    /** set the x data set (predictors). If the state of the runner change when
     *  a new x data set is set, the user of this class have to overload the
     *  udpate() method.
     *  @param x The x data set to run
     */
    virtual void setX( TX const& x)
    {
      p_x_ = &x;
      updateX();
    }

    /** set the data set. If the state of the runner change when a new data
     *  set is set, the user of this class have to overload the udpate() method.
     *  @param y The y data set to run
     */
    virtual void setY( TY const& y)
    {
      p_y_ = &y;
      updateY();
    }

    /** set the data set. If the state of the runner change when a new data
     *  set is set, the user of this class have to overload the udpate() method.
     *  @param y The y data set to run
     *  @param x The x data set to run
     */
    virtual void setData( TX const& y, TY const& x)
    {
      p_y_ = &y;
      p_x_ = &x;
      update();
    }

    /** run the weighted computations.
     *  @return @c true if no error occur during the running process, @c false
     *  otherwise
     **/
    virtual bool run() =0;
    /** run the weighted computations.
     *  @param weights the weights of the samples
     *  @return @c true if no error occur during the running process, @c false
     *  otherwise
     **/
    virtual bool run( typename TY::TContainerVe const& weights) =0;

  protected:
    /** A pointer on the original data set. */
    TY const* p_y_;
    /** A pointer on the original data set. */
    TY const* p_x_;
    /** @brief update the runner when y data set is set.
     * This virtual method will be called when the state of the runner will
     * change, i.e. when a new y data is set is set. By default do nothing.
     **/
    virtual void updateY() {}
    /** @brief update the runner when x data set is set.
     *  This virtual method will be called when the state of the runner will
     *  change, i.e. when a new x data set is set. By default do nothing.
     **/
    virtual void updateX() {}
    /** update the runner.
     * This virtual method will be called when the state of the runner will
     * change, i.e. when new x and y data sets are set. By default do nothing.
     **/
    virtual void update() { updateX(); updateY();}
};

} // namespace STK

#endif /* STK_IRUNNER_H */
