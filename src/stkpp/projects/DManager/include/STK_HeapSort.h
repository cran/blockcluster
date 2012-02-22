/*--------------------------------------------------------------------*/
/*     Copyright (C) 2006-2007  Serge Iovleff

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
 * Project: DManager
 * Purpose: sorting method acting on Containers
 * Author:   Serge Iovleff, serge.iovleff@stkpp.org
 *
 **/

/** @file STK_HeapSort.h
 *  @brief In this file we define an implementation of the heapsort
 *  algorithm for ITContainer1D containers.
 **/

#ifndef STK_HEAPSORT_H
#define STK_HEAPSORT_H

#include "../../Sdk/include/STK_ITContainer1D.h"
#include "../../Sdk/include/STK_ITContainer2D.h"

#ifdef STK_HEAPSORT_DEBUG
#include "../../STKernel/include/STK_Display1D.h"
#endif

namespace STK
{

/** @ingroup DManager
 *  @brief Sort the container T in ascending order
 *  @param T the container to sort
 **/
template<class TYPE, class TContainer1D>
void heapSort( ITContainer1D<TYPE, TContainer1D>& T)
{
  // number of elements
  const Integer nb_elt = T.size();
  if (nb_elt < 2) return;

  // if the container is base one, shift0 = 0 and shift1 = 1
  Integer shift1 = T.first(), shift0 = T.first() - 1;

  // create heap
  for (Integer first = nb_elt/2; first > 0; -- first)
  {
    // the value value to insert in the heap
    TYPE value = T[shift0 + first];
    // organize the heap
    Integer i=first, j=2*first;
    while (j <= nb_elt)
    {
      // j+1 is greatest child
      if ( j < nb_elt && T[shift0 + j] < T[shift1 + j] ) j++;
      // we have find a child gt value
      if (value >= T[shift0 + j]) break;
      // else shift the inner value
      T[shift0 + i] = T[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    T[shift0 + i] = value;
  }
#ifdef STK_HEAPSORT_DEBUG
  std::cout << "T=\n" << T.asLeaf() << "\n";
#endif
  // sort T
  for (Integer last = nb_elt;;)
  { // the value to sort
    TYPE value = T[shift0 + last];
    // Put the top of the heap at the end
    T[shift0 + last] = T[shift1];
    // decrease last.  last==1 : we end the job
    if (--last == 1)
    { T[shift1] = value;
      break;
    }
    // organize the heap
    Integer i=1, j=2;
    while (j <= last)
    { // j+1 is greatest child
      if ( j < last && T[shift0 + j] < T[shift1 + j] ) j++;
      // we have find a child gt value
      if (value >= T[shift0 + j]) break;
      // else shift the inner value
      T[shift0 + i] = T[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    T[shift0 + i] = value;
  }
}

/** @ingroup DManager
 *  @brief Sort the container T in ascending order and return the result in
 *  the container Tsort
 *  @param T the container to sort
 *  @param Tsort the container with the result
 **/
template<class TYPE, class TContainer1D>
void heapSort( const ITContainer1D<TYPE, TContainer1D>& T
             ,  ITContainer1D<TYPE, TContainer1D>& Tsort
             )
{
  // copy T in Tsort
  Tsort.asLeaf() = T.asLeaf();
  // number of elements
  const Integer nb_elt = Tsort.size();
  if (nb_elt < 2) return;

  // if the container is base one, shift0 = 0 and shift1 = 1
  Integer shift1 = Tsort.first(), shift0 = Tsort.first() - 1;

  // create heap
  for (Integer first = nb_elt/2; first > 0; -- first)
  {
    // the value value to insert in the heap
    TYPE value = Tsort[shift0 + first];
    // organize the heap
    Integer i=first, j=2*first;
    while (j <= nb_elt)
    {
      // j+1 is greatest child
      if ( j < nb_elt && Tsort[shift0 + j] < Tsort[shift1 + j] ) j++;
      // we have find a child gt value
      if (value >= Tsort[shift0 + j]) break;
      // else shift the inner value
      Tsort[shift0 + i] = Tsort[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    Tsort[shift0 + i] = value;
  }
#ifdef STK_HEAPSORT_DEBUG
  std::cout << "T=\n" << Tsort.asLeaf() << "\n";
#endif
  // sort T
  for (Integer last = nb_elt;;)
  { // the value to sort
    TYPE value = Tsort[shift0 + last];
    // Put the top of the heap at the end
    Tsort[shift0 + last] = Tsort[shift1];
    // decrease last.  last==1 : we end the job
    if (--last == 1)
    { Tsort[shift1] = value;
      break;
    }
    // organize the heap
    Integer i=1, j=2;
    while (j <= last)
    { // j+1 is greatest child
      if ( j < last && Tsort[shift0 + j] < Tsort[shift1 + j] ) j++;
      // we have find a child gt value
      if (value >= Tsort[shift0 + j]) break;
      // else shift the inner value
      Tsort[shift0 + i] = Tsort[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    Tsort[shift0 + i] = value;
  }
}

/** @ingroup DManager
 * Sort the container T in ascending order using index array.
 *  T is not modified, I contain the indices of the elements of T
 *  in ascending order.
 *  @param I the index array sorting T
 *  @param T the container to sort
 **/
template<class TYPE, class TContainer1D, class TContainer1DInt>
void heapSort( ITContainer1D<Integer, TContainer1DInt> & I
             , ITContainer1D<TYPE, TContainer1D> const& T
             )
{
  // number of elements
  Integer nb_elt = T.size();

  // create index array
  I.resize(T.range());
  Integer first = I.first(), last = I.last();
  for (Integer i=first; i<=last; i++)
  { I[i] = i;}

  if (nb_elt < 2) return;

  // if the container is base one, shift0 = 0 and shift1 = 1
  Integer shift1 = T.first(), shift0 = T.first() - 1;

  // create heap
  for (first = nb_elt/2; first > 0; --first)
  {
    // the value value to insert in the heap
    TYPE value = T[I[shift0 + first]];
    // organize the heap
    Integer i=first, j=2*first;
    while (j <= nb_elt)
    {
      // j+1 is greatest child
      if ( j < nb_elt && T[I[shift0 + j]] < T[I[shift1 + j]] ) j++;
      // we have find a child lt value
      if (value >= T[I[shift0 + j]]) break;
      // else shift the inner values
      I[shift0 + i] = I[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    I[shift0 + i] = shift0 + first;
  }
#ifdef STK_HEAPSORT_DEBUG
  std::cout << "I=\n" << I <<"\n";
#endif
  // sort T
  for (Integer last = nb_elt;;)
  {
    // the value to sort
    Integer ivalue = I[shift0 + last];
    TYPE value = T[ivalue];
    // Put the top of the heap at the end
    //T[shift0 + last] = T[shift1];
    I[shift0 + last] = I[shift1];
    // decrease last.  last==1 : we end the job
    if (--last == 1)
    { //T[shift1] = value;
      I[shift1] = ivalue;
      break;
    }
    // organize the heap
    Integer i=1, j=2;
    while (j <= last)
    { // j+1 is greatest child
      if ( j < last && T[I[shift0 + j]] < T[I[shift1 + j]] ) j++;
      // we have find a child gt value
      if (value >= T[I[shift0 + j]]) break;
      // else shift the inner value
      // T[shift0 + i] = T[shift0 + j];
      I[shift0 + i] = I[shift0 + j];
      // go down in the tree
      i = j;
      j*= 2;
    }
    // plug value in its final location
    // T[shift0 + i] = value;
    I[shift0 + i] = ivalue;
  }
#ifdef STK_HEAPSORT_DEBUG
  std::cout << "I=\n" << I <<"\n";
#endif
}


/** @ingroup DManager
 * Apply a sorting index array to the 1D container T.
 *  @param I the index array sorting T
 *  @param T the container to sort
 **/
template<class TYPE, class TContainer1D, class TContainer1DInt>
void applySort( ITContainer1D<TYPE, TContainer1D>& T
              , ITContainer1D<Integer, TContainer1DInt> const& I
              )
{
#ifdef STK_DEBUG
  if (I.range() != T.range())
  { throw runtime_error("In applySort(T, I) "
                        "incompatible lengths\n");
  }
#endif
  TContainer1D A(T.range());
  const Integer first = I.first(), last = I.last();
  for (Integer i=first; i<= last; i++)
  {
    A[i] = T[I[i]];
  }
  T.asLeaf() = A.asLeaf();
}

/** @ingroup DManager
 * Apply a sorting index array to the 2D container T row by row.
 *  @param I the index array sorting T
 *  @param T the container to sort
 **/
template < class TYPE, class TContainer2D, class TContainer1DInt>
void applySort( ITContainer2D< TYPE, TContainer2D>& T
              , ITContainer1D<Integer, TContainer1DInt> const& I
              )
{
#ifdef STK_DEBUG
  if (I.range() != T.rangeVe())
  { throw runtime_error("In applySort(T, I) "
                        "incompatible lengths\n");
  }
#endif
  TContainer2D A(T.rangeVe(), T.rangeHo());
  const Integer first = I.first(), last = I.last();
  for (Integer i=first; i<= last; i++)
  {
    A(i) = T(I[i]);
  }
  T.asLeaf() = A.asLeaf();
}

} // namespace STK

#endif /*STK_HEAPSORT_H*/
