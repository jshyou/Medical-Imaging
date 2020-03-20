// vector.h

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
 * This file defines a template class to handle the most commonly used 
 * numerical operations for vector. It is specifically emphasized that
 * the class defined below is designed for numerical operations, which
 * is quite different from the container class "vector" in the C++ 
 * Standard Template Library. So, some functions might not be working for
 * non-numerical type
 **/

#ifndef _CUBICIMAGING_VECTOR_OF_NUMERIC_TYPE_H_
#define _CUBICIMAGING_VECTOR_OF_NUMERIC_TYPE_H_

#include "utilStruct.h"
#include "array.h"
#include <memory.h>

/////////////////////////////////////////////////////
//! template class for numeric data array
template <class TT>
class Vector : public Array<TT>
{
public:
  Vector() { }
  Vector(int length, const TT *newData = NULL) { init(length, newData); }
  Vector(const Array<TT> & v) { *this = v; }
  //template<class TT1>
  //Vector(Array<TT1> *v) { *this = (*v); }
  ~Vector() { }

public:
  //! retrieve the buffer address
  TT*   data() const { return m_buffer; }

  //! resize the vector
  void  resize(int newLen, TT *newData = 0) { init(newLen, newData); }

public:
  //! assignment operator for the same type of data
  Vector<TT> & operator= (const Array<TT> & v)
  {
    if( &v != this ) { init(v.length(), v.buffer()); }
    return *this;
  }

  //! assignment operator for the different type of data
  template<class TT1>
  void rest(Array<TT1> & v)
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

  //! assign a constant value to all elements 
  Vector<TT> & operator= (const TT v) { *((Array<TT>*)this) = v; return *this; }

  // add another vector
  Vector<TT> & operator+= (const Vector<TT> & v)
  {
    if(v.isEmpty() || isEmpty()) { return *this; }

    TT * pt1 = m_buffer;
    TT * pt2 = v.buffer();
    int i = 0;
    int dim = ::Min(m_length, v.length());
    while( i < dim )
    {
      *pt1 += *pt2;
      pt1++;
      pt2++;
      i++;
    }
    return *this;
  }

  //! add a constant
  Vector<TT> & operator+= (const TT v)
  {
    if(isEmpty()) { return *this; }

    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt += v;
      pt++;
      i++;
    }
    return *this;
  }

  //! subtract by another vector
  Vector<TT> & operator-= (const Vector<TT> & v)
  {
    if(v.isEmpty() || isEmpty()) { return *this; }

    TT * pt1 = m_buffer;
    TT * pt2 = v.buffer();
    int i = 0;
    int dim = ::Min(m_length, v.length());
    while( i < dim )
    {
      *pt1 -= *pt2;
      pt1++;
      pt2++;
      i++;
    }
    return *this;
  }

  //! substract the elements by a constant
  Vector<TT> & operator-= (const TT v)
  {
    if(isEmpty()) { return *this; }

    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt -= v;
      pt++;
      i++;
    }
    return *this;
  }

  //! set all elements to a constant
  Vector<TT> & operator*= (const TT a )
  {
    if(isEmpty()) { return *this; }

    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt *= a;
      pt++;
      i++;
    }
    return *this;
  }

  //! divide the elements by a constant
  Vector<TT> & operator/= (const TT a)
  {
    if( isEmpty() ) { return *this; }

    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt /= a;
      pt++;
      i++;
    }
    return *this;
  }

public:
  //! find the minimum value
  TT Min() 
  {
    if( m_length <= 0 ) return 0;
    TT tt = m_buffer[0];
    TT* pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      if( tt > *pt ) tt = *pt;
      i++;
      pt++;
    }
    return tt;
  }

  //! find the maximum value
  TT Max()
  {
    if( m_length <= 0 ) return 0;
    TT tt = m_buffer[0];
    TT* pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      if( tt < *pt ) tt = *pt;
      i++;
      pt++;
    }
    return tt;
  }

  //! find the minimum and maximum values
  void getMinMax(TT & min, TT & max)
  {
    if( m_length <= 0 ) return ;
    min = m_buffer[0];
    max = m_buffer[0];
    TT* pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      if( min > *pt ) min = *pt;
      if( max < *pt ) max = *pt;
      i++;
      pt++;
    }
  }

  //! find the median value in the array
  TT median()
  {
    if( m_length <= 0 ) return 0;
    Vector<TT> temp(*this);
    TT * da = temp.data();
    heapsort(da, m_length);
    TT median;
    if( m_length%2 ) median = da[m_length/2];
    else median = (da[m_length/2]+da[m_length/2-1])/2;

    return median;
  }

  //! calculate the total summation of all elements
  double sum() 
  {
    double t = 0;
    TT* pt = this->buffer();
    int i = 0;
    while( i < m_length )
    {
      t += *pt++;
      i++;
    }
    return t;
  }

  //! calculate the norm
  double norm() 
  {
    double t = 0.;
    int i = 0;
    TT *pt = m_buffer;
    while( i < m_length )
    {
      t += (*pt)*(*pt);
      pt++;
      i++;
    }
    return sqrt(t/m_length);
  }

  //! calculate the mean
  double mean() 
  {
    if( m_length <= 0 ) return 0;

    int i = 0;
    double t = 0.0;
    TT* pt = m_buffer;
    while( i < m_length )
    {
      t += *pt;
      pt++;
      i++;
    }
    return (t/m_length);
  }

  //! calculate the standard deviation
  double dev() 
  {
    if( m_length <= 1 ) return 0.0;

    int i = 0;
    double t1 = 0.0, t2 = 0.0;
    TT* pt = m_buffer;
    while( i < m_length )
    {
      t1 += *pt;
      t2 += (*pt)*(*pt);
      pt++;
      i++;
    }

    t2 -= (t1/m_length*t1);
    if( t2 < 0.0 ) return 0.0;
    return sqrt(t2/(m_length-1));
  }

  //! calculate the mean and standard deviation
  void getMeanDev(double & mean, double & dev)
  {
    if( m_length <= 1 ) return ;

    int i = 0;
    mean = 0.;
    dev = 0.;
    TT* pt = m_buffer;
    while( i < m_length )
    {
      mean += *pt;
      dev += (*pt)*(*pt);
      pt++;
      i++;
    }
    mean = mean/m_length;
    dev = (dev-mean*mean*m_length)/(m_length-1);
    if( dev < 0.0 ) { dev = 0.0; return; }
    dev = sqrt(dev);
  }

  //! convert the negative elements to zero
  void positive()
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      if( *pt < 0 ) *pt = 0;
      pt++;
      i++;
    }
  }

  //! convert the positive elements to zero
  void negative()
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      if( *pt > 0 ) *pt = 0;
      pt++;
      i++;
    }
  }

  //! clip the elements to the prescribed value range
  void clip(TT low, TT high)
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      if( *pt > high ) *pt = high;
      else if( *pt < low ) *pt = low;
      pt++;
      i++;
    }
  }

  //! calculate the absolute matrix
  void abs()
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      *pt = (*pt>=0)?(*pt):(-*pt);
      pt++;
      i++;
    }
  }

  //! calculate the sqare matrix
  void square()
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      *pt *= *pt;
      pt++;
      i++;
    }
  }

  //! calculate the sqare-root matrix
  void sqroot()
  {
    TT * pt = data();
    int i = 0;
    while( i < m_length )
    {
      *pt = sqrt(*pt);
      pt++;
      i++;
    }
  }
  
  //! calculate the difference between neighboring elements
  Vector<TT> diff()
  {
    Vector<TT> temp(m_length);
    TT * pt = temp.data();
    for(int i=m_length-1; i>0; i--) pt[i] = m_buffer[i] - m_buffer[i-1];
    pt[0] = m_buffer[0];
    return temp;
  }

  //! calculate the integral
  Vector<TT> integral()
  {
    Vector<TT> temp(m_length);
    TT * pt = temp.data();
    pt[0] = m_buffer[0];

    for(int i=1; i<m_length; i++) pt[i] = pt[i-1] + m_buffer[i];

    return temp;
  }
  
  //! calculate the central difference
  Vector<double> differential_1()
  {
    Vector<double> diff(m_length);
    double* dif = diff.data();

    dif[0] = m_buffer[1] - m_buffer[0];
    dif[m_length-1] = m_buffer[m_length-1] - m_buffer[m_length-2];
    for( int i=1; i<m_length-1; i++) dif[i] = (m_buffer[i+1] - m_buffer[i-1])/2.0;

    return diff;
  }

  //! calculate the second-order difference
  Vector<double> differential_2()
  {
    Vector<double> diff(m_length);
    double* dif = diff.data();

    for( int i=1; i<m_length-1; i++) dif[i] = m_buffer[i+1] + m_buffer[i-1] - m_buffer[i] - m_buffer[i];
    dif[0] = diff[1];
    dif[m_length-1] = diff[m_length-2];

    return diff;
  }

private:
  //! initialize the vector attributes
  void init(int newLen, const TT *newData = 0) 
  {
    if( newLen <= 0 ) { return; }

    if( m_buffer == 0 ) { m_buffer = (TT*) calloc(newLen, sizeof(TT)); }
    else { m_buffer = (TT*) realloc(m_buffer, newLen*sizeof(TT)); }
    if( m_buffer == 0 ) { m_length = 0; return; }
    m_length = newLen;

    if( newData ) { memcpy(m_buffer, newData, newLen*sizeof(TT)); }
  }
};

// not logic operation
template<class TT>
Vector<TT> operator~ (const Vector<TT> & ve)
{
  int length = ve.length();
  Vector<TT> temp(length);
  TT * te1 = ve.buffer();
  TT * te2 = temp.buffer();

  int i = 0;
  while( i < length )
  {
    if(*te1) *te2 = 0;
    else *te2 = 1;
    te1++;
    te2++;
    i++;
  }
  return temp;
}

//////////////////////////////////////
// Arithmetic operations of vectors //
//////////////////////////////////////

//! addition operation between two vectors
template<class TT>
Vector<TT> operator+ (const Vector<TT> & v1, const Vector<TT> & v2)
{
  int len = ( v1.length() <= v2.length() ) ? v1.length() : v2.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = v2.buffer();
  TT *te3 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te3 = *te1 + *te2;
    te3++;
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! addition operation between a constant and a vector
template<class TT>
Vector<TT> operator+ (const Vector<TT> & v1, const TT v2)
{
  int len = v1.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te2 = (*te1 + v2);
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! addition operation between a constant and a vector
template<class TT>
Vector<TT> operator+ (const TT v2, const Vector<TT> & v1)
{
  int len = v1.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te2 = (*te1 + v2);
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! minus operation to a vector
template <class TT>
Vector<TT> operator- (const Vector<TT> & b)
{
  int len = b.length();
  Vector<TT> v(len);
  TT * te = b.buffer();
  TT * te1 = v.buffer();
  int i = 0;
  while( i < len )
  {
    *te1 = -*te;
    te++;
    te1++;
    i++;
  }
  return v;
}

//! subtraction operation between two vectors
template<class TT>
Vector<TT> operator- (const Vector<TT> & v1, const Vector<TT> & v2)
{
  int len = ( v1.length() <= v2.length() ) ? v1.length() : v2.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = v2.buffer();
  TT *te3 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te3 = *te1 - *te2;
    te3++;
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! subtraction operation between a constant and a vector
template<class TT>
Vector<TT> operator- (const Vector<TT> & v1, const TT v2)
{
  int len = v1.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te2 = (*te1 - v2);
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! subtraction operation between a constant and a vector
template<class TT>
Vector<TT> operator- (const TT v2, Vector<TT> v1)
{
  int len = v1.length();
  Vector<TT> vRet( len );
  TT *te1 = v1.buffer();
  TT *te2 = vRet.buffer();

  int i = 0;
  while( i < len )
  {
    *te2 = (v2 - *te1);
    te2++;
    te1++;
    i++;
  }
  return vRet;
}

//! multiply a constant to a vector
template <class TT>
Vector<TT> operator* (const TT a, const Vector<TT> & b)
{
  int len = b.length();
  Vector<TT> v(len);
  TT * t = b.buffer();
  TT * te = v.buffer();
  int i = 0;
  while( i < len )
  {
    *te = ((*t)*a);
    te++;
    t++;
    i++;
  }
  return v;
}

//! multiply a constant to a vector
template <class TT>
Vector<TT> operator* (const Vector<TT> & b, const TT a)
{
  int len = b.length();
  Vector<TT> v(len);

  TT * t = b.buffer();
  TT * te = v.buffer();
  int i = 0;
  while( i < len )
  {
    *te = ((*t)*a);
    te++;
    t++;
    i++;
  }
  return v;
}

//! divide a vector by a constant
template <class TT>
Vector<TT> operator/ (const Vector<TT> & b, const TT a)
{
  int len = b.length();
  Vector<TT> v(len);
  TT * t = b.buffer();
  TT * te = v.buffer();
  int i = 0;
  while( i < len )
  {
    *te = ((*t)/a);
    te++;
    t++;
    i++;
  }
  return v;
}

//! dot-product between two vectors
template<class TT>
double DotProduct(const Vector<TT> & v1, const Vector<TT> & v2)
{
  if( v1.length() != v2.length() ) { return 0.; }
  double tt = 0;
  TT* pt1 = v1.buffer();
  TT* pt2 = v2.buffer();
  int i = 0;
  int len = v1.length();
  while(i < len)
  {
    tt += (*pt1)*(*pt2);
    pt1++;
    pt2++;
    i++;
  }
  return tt;
}

//! cross product between two 3-dimensional vectors
template<class TT>
Vector<double> CrossProduct(const Vector<TT> & v1, const Vector<TT> & v2)
{
  Vector<double> v;

  int len = v1.length();
  if( len != 3 || v2.length() != 3 ) { return v; }
  v.resize(3);

  v[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v[1] = v2[0]*v1[2] - v2[2]*v1[0];
  v[2] = v1[0]*v2[1] - v1[1]*v2[0];

  return v;
}

// if you are tired of using <...>, try this
typedef Vector<char>            CHARVec;
typedef Vector<short>           SHORTVec;
typedef Vector<int>             INTVec;
typedef Vector<long>            LONGVec;
typedef Vector<float>           FLOATVec;
typedef Vector<double>          DOUBLEVec;
typedef Vector<unsigned char>   UCHARVec;
typedef Vector<unsigned short>  USHORTVec;
typedef Vector<unsigned int>    UINTVec;
typedef Vector<unsigned long>   ULONGVec;

#endif // _VECTOR_OF_NUMBER_TYPE_H_
