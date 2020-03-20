// array.h

/**
 * Copyright (c) 2000 - 2005  Cubic Imaging LLC
 *
 * Permission to copy, modify and distribute the software and its documentation 
 * for noncommercial use is hereby granted without fee, provided that the above 
 * copyright notice appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation. Cubic Imaging LLC 
 * makes no representations about the suitability of this software for any purpose. 
 * It is provided "as is" without express or implied warranty.
 * 
 * Permission to copy, modify, and sell the software and its documentation 
 * in commercial use for profits has to be approved by written notice from 
 * Cubic Imaging LLC
 * 
 * \author Jason You (jyou@cubic-imaging.com)
 * \date created at 12/01/2000
 *
 * This file defines a template class to handle the basic array operations.
 **/

#ifndef _CUBICIMAGING_ARRAY_OF_ANY_TYPE_H_
#define _CUBICIMAGING_ARRAY_OF_ANY_TYPE_H_

#include <malloc.h>

//! template class for any type of array
template <class TT>
class Array
{
public: // constructor and destructor
  //! default constructor
  Array() : m_buffer(0), m_length(0), m_ownership(true)	{  }
  
  //! parameterized constructor
  Array(int length, const TT* newData = 0) : m_buffer(0), m_length(0), m_ownership(true) {  init(length, newData); }
  
  //! default copy constructor
  Array(const Array<TT> & v) : m_buffer(0), m_length(0), m_ownership(true) { *this = v; }

  //! copy constructor for other type of array
  //template<class TT1>
  //Array(Array<TT1> *v) : m_buffer(0), m_length(0), m_ownership(true)	{ *this = (*v); }

  //! destructor
  ~Array() { if(m_buffer) empty(); }

public: // accessors
  //! return the number of elements in the array
  int   length() const { return m_length; }

  //! return the address of allocated array buffer
  TT*   buffer() const { return m_buffer; }

  //! return the memory ownership
  bool  ownership() const { return m_ownership; }

public: 
  //! pointer cast operation
  operator TT*() { return m_buffer; }

  //! check if the array is empty
  bool isEmpty() const { return ( m_buffer == 0 ); }

  //! make the array empty
  void empty()
  {
    if( m_buffer && m_ownership ) free(m_buffer);
    m_buffer = 0;
    m_length = 0;
    m_ownership = true;
  }

  //! resize the array
  void resize(int newLen, TT *newData = 0) { init(newLen, newData); }

  //! attach external data to array
  void attach(int newLen, TT *newData)
  {
    empty();
    if( newLen > 0 && newData )
    {
      m_ownership = false;
      m_length = newLen;
      m_buffer = newData;
    }
  }

public:
  //! assignment operator for the same type of array
  Array<TT> & operator= (const Array<TT> & v)
  {
    if( &v != this ) { init(v.length(), v.buffer()); }
    return *this;
  }

  //! assignment operator for other type of array
  template <class TT1>
  void reset( Array<TT1> & v )
  {
    if((void *)&v != (void *)this )
    {
      init(v.length());
      if(m_buffer == 0) return;

      TT * pt1 = m_buffer;
      TT1* pt2 = v.buffer();
      int i = 0;
      while(i < m_length)
      {
        *pt1 = TT(*pt2);
        pt1++;
        pt2++;
        i++;
      }
    }
  }

  //! assign a constant value to array
  Array<TT> & operator= ( const TT t )
  {
    if( m_length <= 0 || m_buffer == 0) return *this;

    TT *pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt = t;
      pt++;
      i++;
    }
    return *this;
  }

  //! read access to elements
  const TT & operator[](int index) const { return *(m_buffer+index); }

  //! write access to elements
  TT & operator[](int index) { return *(m_buffer+index); }

  //! return a sub-array
  Array<TT> operator()(int start, int len)
  {
    if( len <= 0 || start < 0 || start + len > m_length || isEmpty() ) { return Array<TT>(0); }
    else { return Array<TT>(len, m_buffer+start); }
	}

  //! reset the values of a segment in the array
  void reset(int start, const Array<TT> & ar)
  {
    if( start < 0 || start >= m_length ) return;
    int len = m_length-start;
    len = (len < ar.length()? len : ar.length());
    TT * te = ar.buffer();
    TT * te1 = m_buffer + start;
    int i = 0;
    while( i < len )
    {
      *te1 = *te;
      te++;
      te1++;
      i++;
    }
  }

  //! check if the elements are sorted in ascending order
  bool isAscend()
  {
    if( m_length <= 1 ) return true;
    int i = 1;
    while( i < m_length )
    {
      if(m_buffer[i] < m_buffer[i-1]) { return false; }
      i++;
    }
    return true;
  }

  //! check if the elements are sorted in descending order
  bool isDescend()
  {
    if( m_length <= 1 ) return true;
    int i = 1;
    while( i < m_length )
    {
      if(m_buffer[i] > m_buffer[i-1]) { return false; }
      i++;
    }
    return true;
  }

  //! sort the array in ascending order
  void ascend() { heapsort(m_buffer, m_length); }

  //! sort the array in descending order
  void descend() { heapsort(m_buffer, m_length, false); }

  //! periodically shift the elements
  void shift(int dist)
  {
    if( dist < 0 ) dist = ((-dist/m_length+1)*m_length+dist);
    Array<TT> temp(*this);
    TT * te = temp.buffer();
    for(int i=0; i<m_length; i++) m_buffer[i] = te[(i+dist)%m_length];
  }

  //! flip the elements
  void flip()
  {
    TT te;
    for(int i=0; i<m_length/2; i++)
    {
      te = m_buffer[i];
      m_buffer[i] = m_buffer[m_length-i-1];
      m_buffer[m_length-i-1] = te;
    }
  }

  //! find the indices with a specific value 
  Array<unsigned int> find(const TT e)
  {
    int index = 0;
    int i = 0;
    while( i < m_length ) { if(m_buffer[i] == e) index++; i++; }
    Array<unsigned int> indices(index);
    unsigned int * te = indices.buffer();
    i = 0;
    index = 0;
    while( i < m_length ) { if(m_buffer[i] == e) { te[index] = i; index++; } i++; }
    return indices;
  }

  //! find the indices with values within certain range
  Array<unsigned int> find(TT low, TT high)
  {
    int index = 0;
    int i = 0;
    while( i < m_length ) { if(m_buffer[i] >= low && m_buffer[i] <= high) index++; i++; }
    Array<unsigned int> indices(index);
    unsigned int * te = indices.buffer();
    i = 0;
    index = 0;
    while( i < m_length ) { if(m_buffer[i] >= low && m_buffer[i] <= high) { te[index] = i; index++; } i++; }
    return indices;
  }

protected: // attributes
  TT*   m_buffer;
  int   m_length;
  bool  m_ownership;
 
private:
  //! initialize the array attributes
  void init(int newLen, const TT *newData = 0) 
  {
    if( newLen <= 0 ) { return; }

    if( m_buffer == 0 ) { m_buffer = (TT*) calloc(newLen, sizeof(TT)); }
    else { m_buffer = (TT*) realloc(m_buffer, newLen*sizeof(TT)); }
    if( m_buffer == 0 ) { m_length = 0; return; }
		
    m_length = newLen;
    if( newData )
    {
      int i = 0;
      TT* pt1 = m_buffer;
      TT* pt2 = (TT*)newData;
      while( i < m_length )
      {
        *pt1 = *pt2;
        pt1++;
        pt2++;
        i++;
      }
    }
  }
};

//! check if two arrays are equal
template <class TT>
bool operator!= ( const Array<TT> & v1, const Array<TT> & v2 )
{
  if( v1.length() != v2.length() ) return true;
  if( &v1 == &v2 ) return false;
  if( v1.length() == 0 ) return false;

  TT* te1 = v1.buffer();
  TT* te2 = v2.buffer();
  int i = 0;
  int len = v1.length();
  while(i < len)
  {
    if( *te1 != *te2 ) return true;
    te1++;
    te2++;
    i++;
  }
  return false;
}

//! check if two arrays are diffferent
template <class TT>
bool operator== ( const Array<TT> & v1, const Array<TT> & v2 )
{
  return ( !( v1 != v2 ) );
}

#endif	// _CUBICIMAGING_ARRAY_TYPE_H