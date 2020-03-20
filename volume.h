// volumn.h

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
 * \date created at 12/10/2000
 *
 * This file defines a template class to manage a three-dimensional dataset.
 */

#ifndef _CUBICIMAGING_VOLUME_OF_NUMERIC_TYPE_H_
#define _CUBICIMAGING_VOLUME_OF_NUMERIC_TYPE_H_

#include "matrix.h"

//! template class for volumetric data
template <class TT>
class Volume : public Vector<TT>
{
public:
  Volume() : m_depth(0), m_height(0), m_width(0), m_data(NULL) { }
  Volume(int dep, int hei, int wid, TT*** newData = NULL)
    : Vector<TT>(dep*hei*wid), m_depth(0), m_height(0), m_width(0), m_data(NULL)
  { 
    init(dep, hei, wid, newData); 
  }
  Volume(const TT* newData, int dep, int hei, int wid )
    : Vector<TT>(dep*hei*wid), m_depth(0), m_height(0), m_width(0), m_data(NULL)
  {
    init(dep, hei, wid);
    if(m_buffer != NULL)
      memcpy(m_buffer, newData, sizeof(TT)*m_length);
  }
  Volume( const Volume<TT> & v ) : m_depth(0), m_height(0), m_width(0), m_data(NULL) { *this = v; }
  //template<class TT1>
  //Volume( Volume<TT1> & v ) : m_depth(0), m_height(0), m_width(0), m_data(NULL) { *this = v; }

  ~Volume() { empty(); }

public:
  int depth() const { return m_depth; }
  int height() const { return m_height; }
  int width() const { return m_width; }
  TT*** data() const { return m_data; }
  bool isEmpty() {return (m_data == NULL); }
  void empty()
  {
    if( m_data != NULL )
    {
      delete[] m_data[0];
      delete[] m_data;
      m_data = NULL;
      m_depth = 0;
      m_height = 0;
      m_width = 0;
    }

    Array<TT>::empty();
  }

  void resize(int dep, int hei, int wid, TT *** newData)
  {
    Array<TT>::resize(dep*hei*wid);
    init(dep, hei, wid, newData);
  }

  void resize(const TT * newData, int dep, int hei, int wid)
  {
    Array<TT>::resize(dep*hei*wid);
    init(dep, hei, wid);
    memcpy(m_buffer, newData, sizeof(TT)*m_length);
  }

public:
  operator TT***() { return m_data; }

  Volume<TT> & operator= (const Volume<TT> & vl)
  {
    if(vl.isEmpty()) { empty(); return *this; }
    if( this != &vl )
    {
      if(m_length != vl.length())
      {
        Array<TT>::resize(vl.length());
        init(vl.depth(), vl.height(), vl.width());
      }
      else if(m_depth != vl.depth() || m_height != vl.height() || m_width != vl.width() )
      {
        init(vl.depth(), vl.height(), vl.width());
      }
      memcpy(m_buffer, vl.buffer(), sizeof(TT)*m_length);
    }
    return *this;
  }
/*
  template<class TT1>
  Volume<TT> & operator= (Volume<TT1> & vl)
  {
    if(vl.isEmpty()) { empty(); return *this; }
    if( (void*)this != (void*)&vl )
    {
      if(m_length != vl.length())
      {
        Array<TT>::resize(vl.length());
        init(vl.depth(), vl.height(), vl.width());
      }
      else if(m_depth != vl.depth() || m_height != vl.height() || m_width != vl.width() )
      {
        init(vl.depth(), vl.height(), vl.width());
      }
      
      int i = 0;
      TT * pt1 = m_buffer;
      TT1* pt2 = vl.buffer();
      while( i < m_length )
      {
        *pt1 = TT(*pt2);
        pt1++;
        pt2++;
        i++;
      }
    }
    return *this;
  }
*/  
  Volume<TT> & operator= (const TT v)
  {
    int i = 0;
    TT * pt = m_buffer;
    while( i < m_length )
    {
      *pt = v;
      pt++;
      i++;
    }
    return *this;
  }

public:
  Matrix<TT> getXY(int dep)
  {
    if(dep < 0 || dep >= m_depth) return Matrix<TT>();

    Matrix<TT> temp(m_height, m_width);
    memcpy(temp.buffer(), m_buffer+dep*m_height*m_width, m_height*m_width*sizeof(TT));
    return temp;
  }
  
  Matrix<TT> getYZ(int width)
  {
    if(width < 0 || width >= m_width) return Matrix<TT>();

    Matrix<TT> temp(m_depth, m_height);
    TT **pt = temp.data();
    for(int i=0; i<m_depth; i++)
    {
      for(int j=0; j<m_height; j++)
        pt[i][j] = m_data[i][j][width];
    }
    return temp;
  }
  
  Matrix<TT> getXZ(int height)
  {
    if(height < 0 || height >= m_height) return Matrix<TT>();

    Matrix<TT> temp(m_depth, m_width);
    TT **pt = temp.data();
    for(int i=0; i<m_depth; i++)
    {
      for(int j=0; j<m_width; j++)
        pt[i][j] = m_data[i][height][j];
    }
    return temp;
  }

protected:
  TT***	m_data;
  int		m_depth;
  int		m_height;
  int		m_width;

private:
  void init(int dep, int hei, int wid, TT*** data = NULL)
  {
    if(dep<=0 || hei<=0 || wid<=0) return;
    if( m_buffer == NULL ) return;

    if( m_data != NULL )
    {
      delete[] m_data[0];
      delete[] m_data;
      m_data = NULL;
    }
    m_data = new TT**[dep];
    if( m_data == NULL ) return;

    m_data[0] = new TT*[dep*hei];
    if(m_data[0] == NULL)
    {
      delete[] m_data;
      m_data = NULL;
      return;
    }
    
    int i, j;
    for(i=1; i<dep; i++)
      m_data[i] = m_data[0] + hei*i;
    for(i=0; i<dep; i++)
      for(j=0; j<hei; j++)
        m_data[i][j] = m_buffer + i*hei*wid + j*wid;

    m_depth = dep;
    m_height = hei;
    m_width = wid;

    if( data != NULL )
      for(i=0; i<dep; i++)
        for(j=0; j<hei; j++)
          memcpy(m_data[i][j], data[i][j], sizeof(TT)*wid);
  }
};

// if you are tired of using <...>, try this
typedef Volume<unsigned char>	  UCHARVol;
typedef Volume<unsigned short>  USHORTVol;
typedef Volume<unsigned int>    UINTVol;
typedef Volume<unsigned long>	  ULONGVol;
typedef Volume<char>            CHARVol;
typedef Volume<short>           SHORTVol;
typedef Volume<int>             INTVol;
typedef Volume<long>            LONGVol;
typedef Volume<float>           FLOATVol;
typedef Volume<double>          DOUBLEVol;

#endif //_VOLUMN_OF_NUMERIC_TYPE_H_