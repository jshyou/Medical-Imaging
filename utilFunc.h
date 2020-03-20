// utilfunc.h

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
 * \date created at 12/26/2000
 *
 * This file collects some math constants and utility functions in numerical operations.
 **/

#ifndef _CUBICIMAGING_UTILITY_FUNCTION_H_
#define _CUBICIMAGING_UTILITY_FUNCTION_H_

#include <math.h>

//////////////////////////////////////////////////
static const double pi  = 3.14159265358979324;  // ratio of the circumference of a circle to its diameter
static const double lg2 = 0.30102999566398120;	// log(2)/log(10)
static const double ln2 = 0.69314718055994529;  // log(2.)
static const double eps = 0.00000000000000001;  // accuracy for numerical operations
//////////////////////////////////////////////////

//! minimum of two variables
template<class TT>
inline TT Min(const TT & a, const TT & b) { return (a<b? a:b); }

//! maximum of two variables
template<class TT>
inline TT Max(const TT & a, const TT & b) { return (a<b? b:a); }

//! median of three variables
template<class TT>
inline TT Median(const TT & a, const TT & b, const TT & c)
{ 
  if( a < c ) { if(c < b ) return c; 
  else if(b < a ) return a; else return b; }
  else { if(a < b ) return a; else if(b < c ) return c; else  return b; }
}

//! absolute value of a number
template<class TT>
inline TT Abs(const TT & a) { return (a>0? (a): (-a)); }

//! the sign of a number
template<class TT>
inline int Sign(const TT & a) { return (a<0? (-1): (1)); }

//! square of a number
template<class TT>
inline TT Sqr(const TT & a) { return (a*a); }

//! norm of two numbers (precision-keeping method)
inline double Norm(const double c1, const double c2) 
{
  double a = (c1>=0.? c1 : (-c1));
  double b = (c2>=0.? c2 : (-c2));
  if( a >= b && a > eps) { b = b/a; return (a*sqrt(1+b*b)); }
  else if( b > a && b > eps) { a = a/b; return (b*sqrt(1+a*a)); }
  else return 0.0;
}

//! calculate the power of a numeric number
template<class TT>
inline TT Power(const TT & a, const unsigned short & e)
{
  TT v = a;
  int index = 1;
  while( index < e ) { v *= a; index ++; }
  return v;
}

//! Calculate the radix 2 of a positive number
inline double log2( const double x ) { return (log(x)/ln2); }

//! the nearest integer of a floating number
inline int Nearest( const double & x)
{
  if( x < 0. ) return int(x-0.5);
  else return int(x+0.5);
}

//! the kernel function for linear interpolation
inline double Linear( const double & x ) 
{
  if( x >= 1. || x <= -1. ) return 0.;
  else if( x < 0. ) return (1.+x);
  else return (1.-x);
}

//! the kernel function for bicubic interpolation
inline double Bicubic( const double & x )
{
  if( x >= 2. || x <= -2. ) return 0.;
  else if( x <= 1. && x >= -1.) return ( x>0.? (1.-x*x*(2.-x)) : 1.-x*x*(2.+x));
  else return ( x>0.? (4.+x*(x*(5.-x)-8.)) : (4.+x*(x*(5.+x)+8.)));
}

//! generate random number between 0.0 and 1.0
inline double Rand(int * seed)
{
  // L'Ecuyer constants
  static int m1 = 2147483563, m2 = 2147483399;
  static int a1 = 40014, a2 = 40692;
  static int q1 = 53668, q2 = 52774;
  static int r1 = 12211, r2 = 3791;

  // for recursive call
  static int val, bdTable[32], seed2;
  int index, div;

  // create Bays-Durham shuffle table
  if(*seed <= 0)
  {
    *seed = Max(1, -(*seed));
    seed2 = *seed;
    index = 64;
    while( index >= 0 )
    {
      div = (*seed)/q1;
      *seed = a1*((*seed)%q1) - r1*div;
      if( *seed < 0) *seed += m1;
      if(index < 32)  bdTable[index] = *seed;
      index --;
    }
    val = bdTable[0];
  }

  div =(*seed)/q1;
  *seed = a1*((*seed)%q1) - r1*div;
  if(*seed < 0) *seed += m1;
  div = seed2/q2;
  seed2 = a2*(seed2%q2) - r2*div;
  if(seed2 < 0) seed2 += m2;

  index = val%32;
  val = bdTable[index] - seed2;
  bdTable[index] = *seed;
  if(val < 0) val += m1;
  
  return (double(val)/m1);
}

//! calculate the log of gamma function using Lanczos' method
inline double GammaLn56(const double x)
{
  if( x <= 0.0 ) return -1.;

  // Lanczos' constants
  static double c0 =  1.000000000190015;
  static double c1 =  76.18009172947146;
  static double c2 = -86.50532032941677;
  static double c3 =  24.01409824083091;
  static double c4 = -1.231739572450155;
  static double c5 =  1.208650973866179e-3;
  static double c6 = -5.395239384953000e-6;

  double xx, tmp;
  tmp = x + 5.5;
  tmp -= (x + 0.5)*log(tmp);
  xx = (c0 + c1/(x+1) + c2/(x+2) + c3/(x+3) + c4/(x+4) + c5/(x+5) + c6/(x+6));

  if( x < 1.e-6 ) return (log(2.5066282746310005*xx) - tmp - log(x));
  else return (log(2.5066282746310005*xx/x) - tmp);
} 

//! swap the values for two variables with externally defined swap variable
template<class TT>
inline void swap(TT & a1, TT & a2, TT & tp)
{
  tp = a1; a1 = a2; a2 = tp;
}

//! sort two variables in ascending order with externally defined swap variable
template<class TT>
inline void sort2(TT & a1, TT & a2, TT & tp)
{
  if( a2 < a1 ) { tp = a2; a2 = a1; a1 = tp; }
}

//! sort three variables in ascending order with externally defined swap variable
template<class TT>
inline void sort3(TT & a1, TT & a2, TT & a3, TT & tp)
{
  if( a3 < a1 ) { tp = a3; a3 = a1; a1 = tp; }
  if( a2 < a1 ) { tp = a2; a2 = a1; a1 = tp; }
  else if( a3 < a2 ) { tp = a3; a3 = a2; a2 = tp; }
}

//! search the rank-th element in the array, the elements are sorted before the founded element
template<class TT>
inline TT ranksearch(int rank, TT * data, const int dim, bool isAscending = true)
{
  int index1, index2, mid, start = 0, end = dim-1;
  TT ref, temp;

  if( isAscending )
  {
    while( true )
    {
      if( end <= start+1 )
      {
        if(end == start+1 ) sort2(data[start], data[end], temp);
        return data[rank];
      }
      else
      {
        mid = (start + end) >> 1;
        sort2(data[start+1], data[mid], temp);
        sort3(data[start], data[start+1], data[end], temp);
        ref = data[start+1];
        index1 = start+2;
        index2 = end;
        while( true )
        {
          while(data[index1] < ref) { index1 ++; }
          while(data[index2] > ref) { index2 --; }
          if(index2 < index1 ) break;
          sort2(data[index1], data[index2], temp);
        }
        data[start+1] = data[index2];
        data[index2] = ref;
        if( index2 >= rank ) end = index2-1;
        if( index2 <= rank ) start = index1;
      }
    }
  }
  else
  {
    while( true )
    {
      if( end <= start+1 )
      {
        if(end == start+1 ) sort2(data[end], data[start], temp);
        return data[rank];
      }
      else
      {
        mid = (start + end) >> 1;
        sort2(data[mid], data[start+1], temp);
        sort3(data[end], data[start+1], data[start], temp);
        ref = data[start+1];
        index1 = start+2;
        index2 = end;
        while( true )
        {
          while(data[index1] > ref) { index1 ++; }
          while(data[index2] < ref) { index2 --; }
          if(index2 < index1 ) break;
          sort2(data[index2], data[index1], temp);
        }
        data[start+1] = data[index2];
        data[index2] = ref;
        if( index2 >= rank ) end = index2-1;
        if( index2 <= rank ) start = index1;
      }
    }
  }
}

//! perform heap sorting for any type of array if the operator "<" or ">" is defined
template<class TT>
void heapsort(TT * data, const int dim, bool isAscending = true)
{
  if( dim <= 1 )
    return;

  int index1, index2;
  int end = dim;
  int mid = (dim >> 1) + 1;
  TT temp;

  if( isAscending )
  {
    while( true )
    {
      if (mid > 1) { temp = data[--mid-1]; } 
      else 
      {
        temp = data[end-1];
        data[end-1] = data[0];
        if (--end == 1) { data[0] = temp; return; }
      }
      index1 = mid;
      index2 = mid << 1;
      while (index2 <= end) 
      {
        if (index2 < end && data[index2-1] < data[index2]) { ++index2; }
        if (temp < data[index2-1]) 
        {
          data[index1-1] = data[index2-1];
          index2 += (index1 = index2);
        } 
        else { index2 = end + 1; }
      }
      data[index1-1] = temp;
    }
  }
  else
  {
    while( true )
    {
      if (mid > 1) { temp = data[--mid-1]; } 
      else 
      {
        temp = data[end-1];
        data[end-1] = data[0];
        if (--end == 1) { data[0] = temp; return; }
      }
      index1 = mid;
      index2 = mid << 1;
      while (index2 <= end) 
      {
        if (index2 < end && data[index2-1] > data[index2]) { ++index2; }
        if (temp > data[index2-1]) 
        {
          data[index1-1] = data[index2-1];
          index2 += (index1 = index2);
        } 
        else { index2 = end + 1; }
      }
      data[index1-1] = temp;
    }
  }
}

//! inverse 2x2 matrix
inline void InverseMatrix(double mat[2][2], double inversion[2][2])
{
  double det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  if( Abs(det) < eps ) return;

  inversion[0][0] =  mat[1][1]/det;
  inversion[0][1] = -mat[0][1]/det;
  inversion[1][0] = -mat[1][0]/det;
  inversion[1][1] =  mat[0][0]/det;
}

//! inverse 3x3 matrix
inline void InverseMatrix(double mat[3][3], double inversion[3][3])
{
  double det = mat[0][0]*mat[1][1]*mat[2][2] + 
               mat[0][1]*mat[1][2]*mat[2][0] + 
               mat[0][2]*mat[1][0]*mat[2][1] -
               mat[0][2]*mat[1][1]*mat[2][0] -
               mat[0][0]*mat[1][2]*mat[2][1] -
               mat[0][1]*mat[1][0]*mat[2][2]    ;
  if( Abs(det) < eps ) return;

  inversion[0][0] = (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1])/det;
  inversion[0][1] = (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2])/det;
  inversion[0][2] = (mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2])/det;
  inversion[1][0] = (mat[2][0]*mat[1][2] - mat[1][0]*mat[2][2])/det;
  inversion[1][1] = (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0])/det;
  inversion[1][2] = (mat[1][0]*mat[0][2] - mat[0][0]*mat[1][2])/det;
  inversion[2][0] = (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0])/det;
  inversion[2][1] = (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1])/det;
  inversion[2][2] = (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1])/det;
}

///////////////////////////////////////////////////////////////////////////////
// endian exchange operations for data conversion between UNIX and Windows
//
//	swap( x )
//		transform between big-endian and lillte-endian 
//		representation of 16-, 32- or 64-bit numbers.
//
inline unsigned short swap(const unsigned short & source)
{
	return (((source & 0x00FF) << 8) | ((source & 0xFF00) >> 8));
}

inline short swap(const short & source)
{
	return (((source & 0x00FF) << 8) | ((source & 0xFF00) >> 8));
}

inline unsigned int swap(const unsigned int & source)
{
  return  ( ((source & 0x000000FFL) << 24) 
          | ((source & 0x0000FF00L) <<  8)
          | ((source & 0x00FF0000L) >>  8)
          | ((source & 0xFF000000L) >> 24) );
}

inline int swap(const int & source)
{
  return  ( ((source & 0x000000FFL) << 24) 
          | ((source & 0x0000FF00L) <<  8)
          | ((source & 0x00FF0000L) >>  8)
          | ((source & 0xFF000000L) >> 24) );
}

inline unsigned long swap(const unsigned long & source)
{
  return  ( ((source & 0x000000FFL) << 24) 
          | ((source & 0x0000FF00L) <<  8)
          | ((source & 0x00FF0000L) >>  8)
          | ((source & 0xFF000000L) >> 24) );
}

inline long swap(const long & source)
{
  return  ( ((source & 0x000000FFL) << 24) 
          | ((source & 0x0000FF00L) <<  8)
          | ((source & 0x00FF0000L) >>  8)
          | ((source & 0xFF000000L) >> 24) );
}

inline __int64 swap(const __int64 & source )
{
  return  ( ((source & 0x00000000000000FFL) << 56)
          | ((source & 0xFF00000000000000L) >> 56)
          | ((source & 0x000000000000FF00L) << 40)
          | ((source & 0x00FF000000000000L) >> 40) 
          | ((source & 0x0000000000FF0000L) << 24)
          | ((source & 0x0000FF0000000000L) >> 24) 
          | ((source & 0x00000000FF000000L) <<  8)
          | ((source & 0x000000FF00000000L) >>  8) );
}

inline float swap(const float & source)
{
  int source1 = *((int*)&source);
  source1 =   ( ((source1 & 0x000000FFL) << 24) 
              | ((source1 & 0x0000FF00L) <<  8)
              | ((source1 & 0x00FF0000L) >>  8)
              | ((source1 & 0xFF000000L) >> 24) );

  return *(float*)&source1;
}

inline double swap(const double & source )
{
  __int64 source1 = *(__int64*)&source;
  source1 =   ( ((source1 & 0x00000000000000FFL) << 56) 
              | ((source1 & 0xFF00000000000000L) >> 56)
              | ((source1 & 0x000000000000FF00L) << 40)
              | ((source1 & 0x00FF000000000000L) >> 40) 
              | ((source1 & 0x0000000000FF0000L) << 24)
              | ((source1 & 0x0000FF0000000000L) >> 24) 
              | ((source1 & 0x00000000FF000000L) <<  8)
              | ((source1 & 0x000000FF00000000L) >>  8) );

  return *(double*)&source1;
}

////////////////////////////////////////////////////////////////////
#endif	// _CUBICIMAGING_UTILITY_FUNCTION_H