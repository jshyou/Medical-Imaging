// utilstruct.h

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
 * This file collects several fundamental data types such as complex, Cartesian and polar coordinate points.
 **/

#ifndef _CUBICIMAGING_UTILITY_DATA_STRUCT_H_
#define _CUBICIMAGING_UTILITY_DATA_STRUCT_H_

#include "utilfunc.h"

///////////////////////////////////////////////////////////////////////
//!	Class to handle complex number operation
/*
  There is a template class for complex number in the STL. However, the complex number 
  is meaningful only for pair of floating numbers. So, I reimplement it as a regular structure.
*/
struct Complex
{
  double re;
  double im;

  Complex( double real = 0.0, double imag = 0.0 ) :	re( real ), im( imag ) { }
  Complex & operator = ( const double real ) { re = real; im = 0.0; return *this; }
  Complex & operator = ( const Complex & c) { if( &c != this) { re = c.re; im = c.im; } return *this; }

  //! Calculates the conjugate
  Complex conj() { return Complex(this->re, -this->im); }

  //! Calculates the norm of a complex number
	double norm()
  {
    double a = (re>=0.? re : (-re));
    double b = (im>=0.? im : (-im));
    if( a >= b && a > eps) { b = b/a; return (a*sqrt(1+b*b)); }
    else if( b > a && b > eps) { a = a/b; return (b*sqrt(1+a*a)); }
    else return 0.0;
  }

  //! Calculates the phase angle
  double angle()
  {
    double temp = atan2(im, re);
    if(temp >= 0) return temp;
    else return (pi+pi+temp);
  }

  operator double () { return re; }
};

//////////////////////////////////////
//	Complex arithmetic operation
//
inline Complex operator +( const Complex & c, const double a )
{
  return Complex((c.re+a), c.im);
}

inline Complex operator +( const double a, const Complex & c )
{
  return Complex((c.re+a), c.im);
}

inline Complex operator +( const Complex & c1, const Complex & c2 )
{
  return Complex(c1.re + c2.re, c1.im + c2.im);
}

inline Complex operator -( const Complex & c )
{
  return Complex(-c.re, -c.im);
}

inline Complex operator -( const Complex & c, const double a )
{
  return Complex((c.re-a), c.im);
}

inline Complex operator -( const double a, const Complex & c )
{
  return Complex((a-c.re), c.im);
}

inline Complex operator -( const Complex & c1, const Complex & c2 )
{
  return Complex(c1.re - c2.re, c1.im - c2.im);
}

inline Complex operator *( const double a ,const Complex & c )
{
  return Complex((a*c.re), (a*c.im));
}

inline Complex operator *( const Complex & c, const double a )
{
  return Complex((a*c.re), (a*c.im));
}

inline Complex operator *( const Complex & c1, const Complex & c2 )
{
  return Complex(c1.re*c2.re-c1.im * c2.im, c1.re*c2.im+c1.im * c2.re);
}

inline Complex operator /( const Complex & c, const double a )
{
  return Complex(c.re/a, c.im/a);
}

inline Complex operator/ ( const double a, Complex & c )
{
  double aa = c.norm();
  aa *= aa;
  return Complex(a*c.re/aa, a*c.im/aa);
}

inline Complex operator/ ( const Complex & c1, Complex & c2 )
{
  double aa, bb;
  if( Abs(c2.re) >= Abs(c2.im) )
  {
    aa = c2.im/c2.re;
    bb = c2.re+c2.im*aa;
    return Complex((c1.re+c1.im*aa)/bb, (c1.im-c1.re*aa)/bb);
  }
  else
  {
    aa = c2.re/c2.im;
    bb = c2.re*aa+c2.im;
    return Complex((c1.re*aa+c1.im)/bb, (c1.im*aa-c1.re)/bb);
  }
}

inline bool operator== (const Complex & c1, const Complex & c2)
{
  return (c1.re == c2.re && c1.im == c2.im);
}

inline bool operator!= (const Complex & c1, const Complex & c2)
{
  return (c1.re != c2.re || c1.im != c2.im);
}

//! structures for planar points in Cartesian or Polar coordinates
struct PolarPt;
struct CartePt;

struct PolarPt
{
  double radius;
  double theta;   // in degree, (not radian)

  PolarPt() : radius(0.0), theta(0.0) { }
  PolarPt(const double r, const double t) : radius(r), theta(t) { }
  PolarPt(const PolarPt & p) : radius(p.radius), theta(p.theta) { }

  double x() { return radius*cos(theta*pi/180.); }
  double y() { return radius*sin(theta*pi/180.); }

  operator CartePt();
};

struct CartePt
{
  double x;
  double y;

  CartePt() : x (0.), y (0.) { }
  CartePt(const double x0, const double y0) : x(x0), y(y0) { }
  CartePt(const CartePt & p) : x(p.x), y(p.y) { }

  double radius() { return sqrt(x*x+y*y); }
  double theta() { double sl = atan2(y, x); if(sl < 0) sl += (pi+pi); return (sl*180/pi); }

  operator PolarPt();
};

inline CartePt::operator PolarPt()
{
  return PolarPt(this->radius(), this->theta());
}

inline PolarPt::operator CartePt()
{
  return CartePt(this->x(), this->y());
}

inline CartePt operator+ (const CartePt & a, const CartePt & b)
{
  return CartePt(a.x + b.x, a.y + b.y);
}

inline CartePt operator- (const CartePt & a, const CartePt & b)
{
  return CartePt(a.x - b.x, a.y - b.y);
}

inline CartePt operator* (const double a, const CartePt & b)
{
  return CartePt(a*b.x, a*b.y);
}

inline CartePt operator* (const CartePt & b, const double a)
{
  return CartePt(a*b.x, a*b.y);
}

inline CartePt operator/ (const CartePt & b, const double a)
{
  return CartePt(b.x/a, b.y/a);
}

inline double DotProduct(const CartePt & v1, const CartePt & v2)
{
  return (v1.x*v2.x+v1.y*v2.y);
}

inline CartePt LinearTransform(const CartePt & in, const double mat[2][2])
{
  CartePt out;
  out.x = mat[0][0]*in.x + mat[0][1]*in.y;
  out.y = mat[1][0]*in.x + mat[1][1]*in.y;
  return out;
}

inline CartePt LinearTransform(const CartePt & in, const double mat[2][2], const CartePt & shift)
{
  CartePt out;
  out.x = mat[0][0]*in.x + mat[0][1]*in.y + shift.x;
  out.y = mat[1][0]*in.x + mat[1][1]*in.y + shift.y;
  return out;
}

inline bool operator== (const CartePt & a, const CartePt & b)
{
  return (a.x == b.x && a.y == b.y);
}

inline bool operator!= (const CartePt & a, const CartePt & b)
{
  return !(a == b);
}

//! structure for Cartesian coordinates as grid nodes
struct PPoint
{
  int x, y; 

  PPoint() : x(0), y(0) { }
  PPoint(const int x0, const int y0) : x(x0), y(y0) { }
  PPoint(const PPoint & a) : x(a.x), y(a.y) { }

  double norm()	{ return sqrt(double(x*x+y*y)); }
  double angle() { double temp = atan2((double)y, (double)x); if(temp < 0) temp += 2*pi; return temp; }
};

inline PPoint operator+ (const PPoint & a, const PPoint & b)
{
  return PPoint(a.x + b.x, a.y + b.y);
}

inline PPoint operator- (const PPoint & a, const PPoint & b)
{
  return PPoint(a.x - b.x, a.y - b.y);
}

inline PPoint operator* (const int a, const PPoint & b)
{
  return PPoint(a*b.x, a*b.y);
}

inline PPoint operator* (const PPoint & b, const int a)
{
  return PPoint(a*b.x, a*b.y);
}

inline PPoint operator/ (const PPoint & b, const int a)
{
  return PPoint(b.x/a, b.y/a);
}

inline bool operator== (const PPoint & a, const PPoint & b)
{
  return (a.x == b.x && a.y == b.y);
}

inline bool operator!= (const PPoint & a, const PPoint & b)
{
  return !(a == b);
}

//! structures for volumetric points in Cartesian or Sphere coordinates
struct SphereVt;
struct CartesVt;

struct SphereVt
{
  double radius;
  double theta;   // in degree
  double phi;     // in degree

  SphereVt() : radius(0.0), theta(0.0), phi(0.0) { }
  SphereVt(const double r, const double t, const double p): radius(r), theta(t), phi(p) { }
  SphereVt(const SphereVt & vt) : radius(vt.radius), theta(vt.theta), phi(vt.phi) { }

  double x() { double ratio = pi/180; return radius*cos(theta*ratio)*cos(phi*ratio); }
  double y() { double ratio = pi/180; return radius*cos(theta*ratio)*sin(phi*ratio); }
  double z() { return radius*sin(theta*pi/180); }
  operator CartesVt();
};

struct CartesVt
{
  double x, y, z;
  
  CartesVt() : x(0.0), y(0.0), z(0.0) { }
  CartesVt(const double x0, const double y0, const double z0) : x(x0), y(y0), z(z0) { }
  CartesVt(const CartesVt & vt) : x(vt.x), y(vt.y), z(vt.z) { }

  double radius() { return sqrt(x*x+y*y+z*z); }
  double theta();
  double phi();
  operator SphereVt();
};

inline SphereVt::operator CartesVt()
{
  return CartesVt(x(), y(), z());
}

inline CartesVt::operator SphereVt()
{
  return SphereVt(radius(), theta(), phi());
}

inline double CartesVt::theta() 
{ 
  double t = sqrt(x*x+y*y);
  if( t > 0. ) { return (atan(z/t)*180/pi); }
  else
  {
    if(z > 0.) return 90.0;
    else if(z < 0.) return -90.0;
    else return 0.0;
  }
}

inline double CartesVt::phi() 
{ 
  double p = atan2(y, x); 
  if( p < 0. ) p += (pi+pi); 
  return p*180/pi; 
}

inline CartesVt operator+ (const CartesVt & v1, const CartesVt & v2)
{
  return CartesVt(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

inline CartesVt operator- (const CartesVt & v1, const CartesVt & v2)
{
  return CartesVt(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

inline CartesVt operator* (const double s, const CartesVt & v)
{
  return CartesVt(v.x*s, v.y*s, v.z*s);
}

inline CartesVt operator* (const CartesVt & v, const double s)
{
  return CartesVt(v.x*s, v.y*s, v.z*s);
}

inline CartesVt operator/ (const CartesVt & v, const double s)
{
  return CartesVt(v.x/s, v.y/s, v.z/s);
}

inline CartesVt LinearTransform(const CartesVt & in, const double mat[3][3])
{
  CartesVt out;
  out.x = mat[0][0]*in.x + mat[0][1]*in.y + mat[0][2]*in.z;
  out.y = mat[1][0]*in.x + mat[1][1]*in.y + mat[1][2]*in.z;
  out.z = mat[2][0]*in.x + mat[2][1]*in.y + mat[2][2]*in.z;
  return out;
}

inline CartesVt LinearTransform(const CartesVt & in, const double mat[3][3], const CartesVt & shift)
{
  CartesVt out;
  out.x = mat[0][0]*in.x + mat[0][1]*in.y + mat[0][2]*in.z + shift.x;
  out.y = mat[1][0]*in.x + mat[1][1]*in.y + mat[1][2]*in.z + shift.y;
  out.z = mat[2][0]*in.x + mat[2][1]*in.y + mat[2][2]*in.z + shift.z;
  return out;
}

inline double DotProduct(const CartesVt & v1, const CartesVt & v2)
{
  return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
}

inline CartesVt CrossProduct(const CartesVt& v1, const CartesVt& v2)
{
  return CartesVt(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}

inline bool operator== (const CartesVt & a, const CartesVt & b)
{
  return (a.x == b.x && a.y == b.y && a.z == b.z);
}

inline bool operator!= (const CartesVt & a, const CartesVt & b)
{
  return !(a == b);
}

//! structure for Cartesian coordinates as triplet of integer numbers
struct VPoint
{
  int x, y, z; 

  VPoint() : x (0), y (0), z (0)	{ }
  VPoint(const int x0, const int y0, const int z0) : x (x0), y (y0), z (z0)	{ }
  VPoint(const VPoint& vP) : x(vP.x), y(vP.y), z(vP.z) { }

  double norm() { return sqrt((double)x*x+y*y+z*z); }
  double theta();
  double phi();
  int sqr() { return (x*x+y*y+z*z); }
};

inline double VPoint::theta() 
{ 
  double t = sqrt((double)x*x+y*y);
  if( t > 0. ) { return (atan(z/t)*180/pi); }
  else { if(z > 0.) return 90.0; else if(z < 0.) return -90.0; else return 0.0; }
}

inline double VPoint::phi() 
{ 
  double p = atan2((double)y, (double)x); 
  if( p < 0. ) p += (pi+pi); 
  return p*180/pi; 
}

inline VPoint operator+ (const VPoint& v1, const VPoint& v2)
{
  return VPoint(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

inline VPoint operator- (const VPoint& v1, const VPoint& v2)
{
  return VPoint(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

inline VPoint operator* (const int s, const VPoint& v)
{
  return VPoint(v.x*s, v.y*s, v.z*s);
}

inline VPoint operator* (const VPoint& v, const int s)
{
  return VPoint(v.x*s, v.y*s, v.z*s);
}

inline VPoint operator/ (const VPoint& v, const int s)
{
  return VPoint(v.x/s, v.y/s, v.z/s);
}

inline bool operator== (const VPoint & a, const VPoint & b)
{
  return (a.x == b.x && a.y == b.y && a.z == b.z);
}

inline bool operator!= (const VPoint & a, const VPoint & b)
{
  return !(a == b);
}

inline int DotProduct(const VPoint & v1, const VPoint & v2)
{
  return (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
}

inline VPoint CrossProduct(const VPoint & v1, const VPoint & v2)
{
  return VPoint(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}

//! ultility structure for monitoring the duration of a specific operation
#include <time.h>
struct CTimer
{
  clock_t t1, t2;
  double	duration;

  CTimer() : t1(0), t2(0), duration(0.0) { }
  inline void start()	{ t1 = clock(); }
  inline void stop() { t2 = clock(); duration = double(t2-t1)/CLOCKS_PER_SEC; }
};

#endif // _UTILITY_TYPE_H
