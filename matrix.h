// matrix.h

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
 * This file defines a template class to manage two-dimensional matrix data together with 
 * the commonly used numerical operations. Note, some functions might not work for non-number type.
 **/

#ifndef _CUBICIMAGING_MATRIX_OF_NUMERIC_TYPE_H_
#define _CUBICIMAGING_MATRIX_OF_NUMERIC_TYPE_H_

#include "vector.h"

///////////////////////////////////////////////////////////////////////////
//
// Class Matrix
// 
// The following figure illustrates the Carteisan coordinate representation of 
// matrix elements. This rule applies to all coordinate-related operattions.

/*
y^(m_rows)
 |
r|*******************************************(matrix elements)
o|*******************************************
w|*******************************************
 |*******************************************
0---------------------------------------------------->x(m_cols)
      column direction
*/

//! template class for matrix of numberic type
template <class TT>
class Matrix : public Vector<TT>
{
public:
  Matrix() : m_data( 0 ), m_rows( 0 ), m_cols( 0 ) { }
  Matrix(int rows, int cols, TT **newData = 0) : m_data( 0 ), m_rows( 0 ), m_cols( 0 ) { init(rows, cols, newData); }
  Matrix(const TT * newData, int rows, int cols) : m_data( 0 ), m_rows( 0 ), m_cols( 0 ) { init(rows, cols, (TT*)newData); }
  Matrix(const Matrix<TT> & m) : m_data( 0 ), m_rows( 0 ), m_cols( 0 ) { *this = m; }
  //template<class TT1>
  //Matrix(Matrix<TT1> *m) : m_data( 0 ), m_rows( 0 ), m_cols( 0 ) { *this = (*m); }

  ~Matrix() { empty(); }

public:
  int   rows() const { return m_rows; }
  int   cols() const { return m_cols; }
  TT**  data() const { return m_data; }

  operator TT**() { return m_data; }
  bool isEmpty() const { return ( m_data == 0 ); }
  void empty()
  {
    if( m_data ) delete [] m_data;
    m_data = 0;
    m_rows = 0;
    m_cols = 0;

    Array<TT>::empty();
  }

public:
  void resize(int rows, int cols, TT **newData = 0 ) { init(rows, cols, newData); }
  void resize(int rows, int cols, TT *newData ) { init(rows, cols, newData); }
  void attach(int rows, int cols, TT *newData)
  {
    empty();
    m_ownership = false;
    m_rows = rows;
    m_cols = cols;
    m_data = new TT*[rows]; if(m_data == 0) { empty(); return; }
    m_length = rows*cols;
    m_buffer = newData;
    for( int i = 0; i < m_rows; i++ ) m_data[i] = newData + m_cols*i;
  }

public:
  // assignment operator
  Matrix<TT> & operator= (const Matrix<TT> & m)
  {
    if( m.isEmpty() ) { empty(); return *this; }
    if( &m != this )
    {
      if(this->length() != m.length() || m_rows != m.rows() || m_cols != m.cols()) { init( m.rows(), m.cols(), m.buffer()); }
      else { memcpy(m_buffer, m.buffer(), sizeof(TT)*m.length()); }
    }
    return *this;
  }

  template <class TT1>
  void reset( Matrix<TT1> & m )
  {
    if( m.isEmpty() ) { this->empty(); return ; }

    if(this->length() != m.length() || m_rows != m.rows() || m_cols != m.cols()) { init( m.rows(), m.cols()); }

    int len = this->length();
    TT * pt1 = m_buffer;
    TT1 * pt2 = m.buffer();
    int i = 0;
    while( i < len )
    {
      *pt1 = TT(*pt2);
      pt1++;
      pt2++;
      i++;
    }
  }

  Matrix<TT> & operator= ( const TT t )
  {
    TT* pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt = t;
      pt++;
      i++;
    }
    return *this;
  }

  // read/write access to elements
  TT* operator[]( int row )
  {
    if( isEmpty() || row < 0 || row >= m_rows) return 0;
    return *(m_data+row);
  }

  const TT* operator[]( int row ) const
  {
    if( isEmpty() || row < 0 || row >= m_rows) return 0;
    return *(m_data+row);
  }

  Vector<TT> getRow( int row )
  {
    if( isEmpty() || row < 0 || row >= m_rows ) { return Vector<TT>(); }    
    return Vector<TT>(m_cols, m_data[row] );
  }

  Vector<TT> getCol( int col )
  {
    if( isEmpty() || col < 0 || col >= m_cols ) { return Vector<TT>(); }
    Vector<TT> v(m_rows);
    TT * buff = v.buffer();
    for( int i = 0; i < m_rows; i++ ) buff[i] = m_data[i][col];
    return v;
  }

  Matrix<TT> operator()( int startRow, int startCol, int rows, int cols )
  {
    if( isEmpty() || startRow < 0 || rows < 1 || startRow + rows > m_rows 
      || startCol < 0 || cols < 1 || startCol + cols > m_cols )
    {
      return Matrix<TT>();
    }
    Matrix<TT> m( rows, cols );
    TT ** data = m.data();
    for( int i = 0; i < rows; i++ )
    {
      for( int j = 0; j < cols; j++ ) data[i][j] = m_data[i + startRow][j + startCol];
    }
    return m;
  }

  void setRow( int row, const Vector<TT> & ve)
  {
    if( row<0 || row>m_rows-1 || ve.length() <= 0 ) return;

    int len = (m_cols<ve.length()? m_cols : ve.length());
    int i;
    TT * da = ve.buffer();
    for(i=0; i<len; i++) m_data[row][i] = da[i];
  }

  void setCol(int col, const Vector<TT> & ve)
  {
    if( col<0 || col>m_cols-1 || ve.length() <= 0 ) return;
    int len = (m_rows<ve.length()? m_rows : ve.length());
    int i;
    TT *da = ve.buffer();
    for(i=0; i<len; i++) m_data[i][col] = da[i];
  }

  void setBlock(int row, int col, const Matrix<TT> & ma)
  {
    if( row < 0 || col < 0 || row > m_rows-1 || col > m_cols-1 || 
      ma.rows() <= 0 || ma.cols() <= 0)
      return;
    int rowLen = ((m_rows-row)<ma.rows()? (m_rows-row):ma.rows());
    int colLen = ((m_cols-col)<ma.cols()? (m_cols-col):ma.cols());
    if( rowLen <= 0 || colLen <= 0 ) return;
    int i, j;
    TT ** temp = ma.data();
    for(i=0; i<rowLen; i++)
      for(j=0; j<colLen; j++)
        m_data[i+row][j+col] = temp[i][j];
  }

  // arithmatic operations
  Matrix<TT> & operator+= ( const Matrix<TT> & m )
  {
    if( m.isEmpty() )
    {
      return *this;
    }
    int row = (m.rows()<m_rows? m.rows() : m_rows);
    int col = (m.cols()<m_cols? m.cols() : m_cols);
    if( row < 0 || col < 0 ) { return *this; }
    TT ** pt1 = m_data;
    TT ** pt2 = m.data();
    
    int i = 0;
    int j = 0;
    while( i < row )
    {
      j = 0;
      while( j < col )
      {
        pt1[i][j] += pt2[i][j];
        j++;
      }
      i++;
    }
    return *this;
  }

  Matrix<TT> & operator+= ( const TT m )
  {

    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt += m;
      pt++;
      i++;
    }
    return *this;
  }

  Matrix<TT> & operator-= ( const Matrix<TT> & m )
  {
    if( m.isEmpty() ) { return *this; }
    int row = (m.rows()<m_rows? m.rows() : m_rows);
    int col = (m.cols()<m_cols? m.cols() : m_cols);
    if( row < 0 || col < 0 ) { return *this; }
    TT ** pt1 = m_data;
    TT ** pt2 = m.data();
    int i = 0;
    int j = 0;
    while( i < row )
    {
      j = 0;
      while( j < col )
      {
        pt1[i][j] -= pt2[i][j];
        j++;
      }
      i++;
    }
    return *this;
  }

  Matrix<TT> & operator-= ( const TT m )
  {
    TT * pt = m_buffer;
    int i = 0;
    while( i < m_length )
    {
      *pt -= m;
      pt++;
      i++;
    }
    return *this;
  }

  Matrix<TT> & operator*= ( const TT a )
  {
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

  Matrix<TT> & operator/= ( const TT a )
  {
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
  // utility functions
  Matrix<TT> transpose()
  {
    int i, j;
    Matrix<TT> mRet(m_cols, m_rows);
    TT ** te = mRet.data();
    for( i = 0; i < m_rows; i++ )
      for( j = 0; j < m_cols; j++ )
        te[j][i] = m_data[i][j];
    return mRet;
  }

  // Periodically shift matrix elements
  void shift(int dist_x, int dist_y)
  {
    if( dist_x < 0 ) dist_x = ((-dist_x/m_cols+1)*m_cols+dist_x);
    if( dist_y < 0 ) dist_y = ((-dist_y/m_rows+1)*m_rows+dist_y);

    Matrix<TT> temp = *this;
    TT ** te = temp.data();
    int i, j;
    for(i=0; i<m_rows; i++)
      for(j=0; j<m_cols; j++)
        m_data[i][j] = te[(i+dist_y)%m_rows][(j+dist_x)%m_cols];
  }

  void flipRow()
  {
    TT te;
    int i, j;
    for(i=0; i<m_rows/2; i++)
      for(j=0; j<m_cols; j++)
      {
        te = m_data[i][j];
        m_data[i][j] = m_data[m_rows-i-1][j];
        m_data[m_rows-i-1][j] = te;
      }
  }

  void flipCol()
  {
    TT te;
    int i, j;
    for(i=0; i<m_rows; i++)
      for(j=0; j<m_cols/2; j++)
      {
        te = m_data[i][j];
        m_data[i][j] = m_data[i][m_cols-1-j];
        m_data[i][m_cols-1-j] = te;
      }
  }

  void swapRows(int row1, int row2)
  {
    if( row1<0 || row1>=m_rows || row2<0 || row2>=m_rows ) return;

    TT temp;
    for(int i=0; i<m_cols; i++)
    {
      temp = m_data[row1][i];
      m_data[row1][i] = m_data[row2][i];
      m_data[row2][i] = temp;
    }
  }

  void swapCols(int col1, int col2)
  {
    if( col1<0 || col1>=m_cols || col2<0 || col2>=m_cols ) return;

    TT temp;
    for(int i=0; i<m_rows; i++)
    {
      temp = m_data[i][col1];
      m_data[i][col1] = m_data[i][col2];
      m_data[i][col2] = temp;
    }
  }

  void invert() // only works for double square matrix
  {
    if( this->isEmpty() || m_rows != m_cols ) return;

    int n = m_rows;
    double** mI = data();

    // Allocate index arrays and set to zero
    int *pivotFlag = new int[n];
    int *swapRow = new int[n];
    int *swapCol = new int[n];
    if( pivotFlag == 0 || swapRow == 0 || swapCol == 0 ) return;
    memset(pivotFlag, 0, n*sizeof(int));
    memset(swapRow, 0, n*sizeof(int));
    memset(swapCol, 0, n*sizeof(int));

    // Pivoting n iterations
    int i, row, irow, col, icol, swap;
    double big;
    double abs_element;
    double temp;
    double pivotInverse;
    for( i = 0; i < n; i++ )
    {
      big = 0.0;
      for( row = 0; row < n; row++ )
      {
        if( pivotFlag[row] == 0 )
        {
          for( col = 0; col < n; col++ )
          {
            if( pivotFlag[col] == 0 )
            {
              abs_element = fabs( mI[row][col] );
              if( abs_element >= big )
              {
                big = abs_element;
                irow = row;
                icol = col;
              }
            }
          }
        }
      }

      pivotFlag[icol]++;

      if( irow != icol )
      {
        for( col = 0; col < n; col++ )
        {
          temp = mI[irow][col];
          mI[irow][col] = mI[icol][col];
          mI[icol][col] = temp;
        }
      }

      swapRow[i] = irow;
      swapCol[i] = icol;

      // Bad news if the pivot is zero
      if( mI[icol][icol] == 0.0 )
        return;

      pivotInverse = 1.0 / mI[icol][icol];

      mI[icol][icol] = 1.0;
      for( col = 0; col < n; col++ )
        mI[icol][col] *= pivotInverse;

      for( row = 0; row < n; row++ )
      {
        if( row != icol )
        {
          temp = mI[row][icol];
          mI[row][icol] = 0.0;
          for( col = 0; col < n; col++)
            mI[row][col] -= mI[icol][col] * temp;
        }
      }
    }

    for( swap = n-1; swap >= 0; swap-- )
    {
      if( swapRow[swap] != swapCol[swap] )
      {
        for( row = 0; row < n; row++ )
        {
          temp = mI[row][swapRow[swap]];
          mI[row][swapRow[swap]] = mI[row][swapCol[swap]];
          mI[row][swapCol[swap]] = temp;
        }
      }
    }

    delete [] pivotFlag;
    delete [] swapRow;
    delete [] swapCol;
  }

  // Calculates the determinant of a square matrix
  double det()
  {
    if( this->isEmpty() || m_rows != m_cols ) return 0.;

    int n = m_rows;
    Matrix<TT> tempMat = *this;
    TT** mDet = tempMat.data();

    // Initialize the answer
    int swapRow;
    double big;
    double abs_element;
    double pivotInverse;
    double temp;
    int row, col, pivot;
    double det = 1.0;

    for( pivot = 0; pivot < n-1; pivot++ )
    {
      big = fabs( mDet[pivot][pivot] );
      swapRow = 0;
      for( row = pivot + 1; row < n; row++ )
      {
        abs_element = fabs( mDet[row][pivot] );
        if( abs_element > big )
        {
          swapRow = row;
          big = abs_element;
        }
      }
      if( swapRow != 0 )
      {
        swapRows( pivot, swapRow );
        det *= -mDet[pivot][pivot];
      }
      else
      {
        det *= mDet[pivot][pivot];
      }

      if( fabs( det ) < 1.0e-50 )
        return det;

      pivotInverse = 1.0 / mDet[pivot][pivot];
      for( col = pivot + 1; col < n; col++ )
        mDet[pivot][col] = mDet[pivot][col] * pivotInverse;

      for( row = pivot + 1; row < n; row++ )
      {
        temp = mDet[row][pivot];
        for( col = pivot + 1; col < n; col++ )
        {
          mDet[row][col] = mDet[row][col] - mDet[pivot][col] * temp;
        }
      }
    }

    det *= mDet[n-1][n-1];
    return det;  
  }

  double trace()
  {
    if( this->isEmpty() || m_rows != m_cols ) return (-1.e300);

    double sum = 0.0;
    int i = 0;

    while( i < m_rows )
    {
      sum += m_data[i][i];
      i++;
    }

    return sum;
  }

//////////////////////////////////////////////////////////////////
// Diffrential operation as a two-dimensional discrete function //
//////////////////////////////////////////////////////////////////

  Matrix<TT> diff_col()
  {
    int i, j;
    Matrix<TT> fx(m_rows, m_cols);
    TT** fx_data = fx.data();

    for( i=0; i<m_rows; i++ )
    {
      for( j=1; j<m_cols-1; j++ )
      {
        fx_data[i][j] = ((m_data[i][j+1] - m_data[i][j-1])/2);
      }
    }

    for( i=0; i<m_rows; i++)
    {
      fx_data[i][0] = m_data[i][1]-m_data[i][0];
      fx_data[i][m_cols-1] = m_data[i][m_cols-1] - m_data[i][m_cols-2];
    }

    return fx;
  }

  Matrix<TT> diff_row()
  {
    int i, j;
    Matrix<TT> fy(m_rows, m_cols);
    TT** fy_data = fy.data();
    
    for( j=0; j<m_cols; j++ )
    {
      for( i=1; i<m_rows-1; i++ )
      {
        fy_data[i][j] = ((m_data[i+1][j] - m_data[i-1][j])/2);
      }
    }

    for( j=0; j<m_cols; j++)
    {
      fy_data[0][j] = m_data[1][j] - m_data[0][j];
      fy_data[m_rows-1][j] = m_data[m_rows-1][j] - m_data[m_rows-2][j];
    }

    return fy;
  }

  Matrix<TT> integral()
  {
    int i, j;
    Matrix<TT> sum(m_rows, m_cols);
    TT** sum_data = sum.data();

    for( i=0; i<m_rows; i++ )
    {
      sum_data[i][0] = m_data[i][0];
      for( j=1; j<m_cols; j++ )
      {
        sum_data[i][j] = sum_data[i][j-1] + m_data[i][j];
      }
    }
    for( j=0; j<m_cols; j++ )
    {
      for( i=1; i<m_rows; i++)
      {
        sum_data[i][j] = sum_data[i-1][j] + sum_data[i][j];
      }
    }

    return sum;
  }

  // Calculate the first partial gradient using central difference
  Matrix<double> differential_col1()
  {
    int i, j;
    Matrix<double> fx(m_rows, m_cols);
    double** fx_data = fx.data();

    for( i=0; i<m_rows; i++ )
    {
      for( j=1; j<m_cols-1; j++ )
      {
        fx_data[i][j] = ((m_data[i][j+1] - m_data[i][j-1])/2.0);
      }
    }

    for( i=0; i<m_rows; i++)
    {
      fx_data[i][0] = m_data[i][1]-m_data[i][0];
      fx_data[i][m_cols-1] = m_data[i][m_cols-1] - m_data[i][m_cols-2];
    }

    return fx;
  }

  Matrix<double> differential_row1()
  {
    int i, j;
    Matrix<double> fy(m_rows, m_cols);
    double** fy_data = fy.data();
    
    for( j=0; j<m_cols; j++ )
    {
      for( i=1; i<m_rows-1; i++ )
      {
        fy_data[i][j] = ((m_data[i+1][j] - m_data[i-1][j])/2.0);
      }
    }

    for( j=0; j<m_cols; j++)
    {
      fy_data[0][j] = m_data[1][j] - m_data[0][j];
      fy_data[m_rows-1][j] = m_data[m_rows-1][j] - m_data[m_rows-2][j];
    }

    return fy;
  }

  // Calculate the second derivative using three point-difference
  Matrix<double> differential_col2()
  {
    int i, j;
    Matrix<double> fxx(m_rows, m_cols);
    double** fxx_data = fxx.data();

    for( i=0; i<m_rows; i++ )
    {
      for( j=1; j<m_cols-1; j++ )
      {
        fxx_data[i][j] = m_data[i][j+1] + m_data[i][j-1] -m_data[i][j] - m_data[i][j];
      }
    }

    for( i=0; i<m_rows; i++)
    {
      fxx_data[i][0] = fxx.data()[i][1];
      fxx_data[i][m_cols-1] = fxx.data()[i][m_cols-2];
    }

    return fxx;
  }

  Matrix<double> differential_row2()
  {
    int i, j;
    Matrix<double> fyy(m_rows, m_cols);
    double** fyy_data = fyy.data();

    for( j=0; j<m_cols; j++ )
    {
      for( i=1; i<m_rows-1; i++ )      
      {
        fyy.data()[i][j] = m_data[i+1][j] + m_data[i-1][j] -m_data[i][j] - m_data[i][j];
      }
    }

    for( j=0; j<m_cols; j++ )
    {
      fyy.data()[0][j] = fyy.data()[1][j];
      fyy.data()[m_rows-1][j] = fyy.data()[m_rows-2][j];
    }

    return fyy;
  }

  Matrix<double> differential_cross()
  {
    int i, j;
    Matrix<double> fxy(m_rows, m_cols);
    double fxy_data = fxy.data();

    for( i=1; i<m_rows-1; i++ )
    {
      for( j=1; j<m_cols-1; j++ )            
      {
        fxy_data[i][j] = (m_data[i+1][j+1] + m_data[i-1][j-1] -m_data[i+1][j-1] - m_data[i-1][j+1])/4.0;
      }
    }
    for( i=0; i<m_rows; i++)
    {
      fxy_data[i][0] = fxy[i][1];
      fxy_data[i][m_cols-1] = fxy[i][m_cols-2];
    }
    for( j=0; j<m_cols; j++)
    {
      fxy_data[0][j] = fxy[1][j];
      fxy_data[m_rows-1][j] = fxy[m_rows-2][j];
    }

    return fxy;
  }

////////////////////////////////////////////////////////
// coordinate transform as a discrete planar function //
////////////////////////////////////////////////////////

  // assume that the original matrix is in Cartesian coordinate system
  // then transform into polar coordinate representation
  // the m_row of transformed matrix represents radius, m_width represents angle
  Matrix<TT> getPolar(const int numOfRays)
  {
    if( m_rows == 0 || m_cols == 0 || numOfRays <= 0) return Matrix<TT>();

    int radius = (int)(sqrt(rows()*rows()+cols()*cols())/2) + 1;

    int i, j;

    Matrix<TT> polar(radius, numOfRays);
    TT min = this->Min();

    Vector<double> cos_(numOfRays);
    Vector<double> sin_(numOfRays);

    for( i=0; i<numOfRays; i++)
    {
      cos_[i] = cos(2*pi*(i+0.5)/numOfRays);
      sin_[i] = sin(2*pi*(i+0.5)/numOfRays);
    }

    double x, y;
    double cx, cy;
    int xx, yy;

    for( j=0; j<radius; j++)
    {
      for( i=0; i<numOfRays; i++)
      {
        x = (j+0.5)*cos_[i] + m_cols/2-0.5;
        xx = (int)floor(x);
        cx = x - xx;
        y = (j+0.5)*sin_[i] + m_rows/2-0.5;
        yy = (int)floor(y);
        cy = y - yy;

        if( xx >= 0 && xx+1 < m_cols && yy >= 0 && yy+1 < m_rows )
          polar[j][i] = TT((m_data[yy][xx] + cx*(m_data[yy][xx+1]-m_data[yy][xx]) +
          cy*(m_data[yy+1][xx]-m_data[yy][xx] + 
          cx*(m_data[yy+1][xx+1]+m_data[yy][xx]-m_data[yy+1][xx]-m_data[yy][xx+1]))));
        else
          polar[j][i] = min;
      }
    }

    return polar;
  }

  // assume that the original matrix is in Polar coordinate system
  // then transform into Cartesian coordinate representation
  Matrix<TT> getCart()
  {
    // set the row and col for cartesian coodinate, the new Mat is square
    int row = 2*rows();
    int col = 2*rows();

    Matrix<TT> cart(row, col);
    TTmin = this->Min();

    TT ** ca = cart.data();

    int i, j;
    double rad;
    double ang;
    int rr, an;
    CPPoint te;

    double anGrid;
    if(m_cols) anGrid = 2*pi/m_cols;
    else anGrid = pi;

    for( i=0; i<row; i++)
    {
      for( j=0; j<col; j++)
      {
        te = CPPoint(j-col/2.+0.5, i-row/2.+0.5).polar();
        rad = te.x;
        ang = te.y/anGrid;

        rr = (int)floor(rad);
        rad = rad - rr;

        an = (int)floor(ang);
        ang = ang - an;

        if(an == m_cols-1)
        {
          if( rr < m_rows-1 )
            ca[i][j] = (TT)(m_data[rr][an] + rad*(m_data[rr+1][an]-m_data[rr][an])+
              ang*((m_data[rr][0]-m_data[rr][an]) +
              rad*(m_data[rr+1][0]+m_data[rr][an]-m_data[rr+1][an]-m_data[rr][0])));
          else if( rr == m_rows-1 )
            ca[i][j] = (TT)(m_data[rr][an] + ang*(m_data[rr][0]-m_data[rr][an]));
          else
            ca[i][j] = min;
        }
        else
        {
          if( rr < m_rows-1 )
            ca[i][j] = (TT)(m_data[rr][an] + rad*(m_data[rr+1][an]-m_data[rr][an])+
              ang*((m_data[rr][an+1]-m_data[rr][an]) +
              rad*(m_data[rr+1][an+1]+m_data[rr][an]-m_data[rr+1][an]-m_data[rr][an+1])));
          else if( rr == m_rows-1 )
            ca[i][j] = (TT)(m_data[rr][an] + ang*(m_data[rr][an+1]-m_data[rr][an]));
          else
            ca[i][j] = min;
        }
      }
    }

    return cart;
  }

  Vector<double> getRowProj()
  {
    Vector<double> rProj(m_rows);
    double* pt = rProj.buffer();

    int i, j;
    for(i=0; i<m_rows; i++)
    {
      for(j=0; j<m_cols; j++)
        pt[i] += m_data[i][j];
    }
    return rProj;
  }

  Vector<double> getColProj()
  {
    Vector<double> cProj(m_cols);
    double* pt = cProj.buffer();

    int i, j;
    for(i=0; i<m_rows; i++)
    {
      for(j=0; j<m_cols; j++)
        pt[j] += m_data[i][j];
    }
    return cProj;
  }

  Matrix<double> getProjection(const int numOfRays)
  {
    double step = pi/numOfRays;
    Matrix<double> result(numOfRays, m_rows);
    double** pt = result.data();

    Matrix<TT> rotated(m_rows, m_cols);
    TT ** rot = rotated.data();

    double arc;
    double cos_, sin_;
    CPPoint original, transformed;
    int xx, yy;
    double alpha_x, alpha_y;
    double half_row = m_rows/2.;
    double half_col = m_cols/2.;
    double aa, bb;
    int m, i, j;
    for(int m=0; m<numOfRays; m++)
    {
      arc = m*step;
      cos_ = cos(-arc);
      sin_ = sin(-arc);

      for( i = 0; i < m_rows; i++)
      {
        for( j = 0; j < m_cols; j++)
        {
          transformed.x = j - half_col + 0.5;
          transformed.y = i - half_row + 0.5;
          original.x = transformed.x*cos_ - transformed.y*sin_;
          original.y = transformed.x*sin_ + transformed.y*cos_;

          aa = original.x + half_col - 0.5;
          xx = (int)floor(aa);
          alpha_x = aa - xx;

          bb = original.y + half_row - 0.5;
          yy = (int)floor(bb);
          alpha_y = bb - yy;

          if( yy < m_rows-1 && yy >=0 && xx < m_cols-1 && xx >= 0)
            rot[i][j] = TT(st[yy][xx] + alpha_x*(st[yy][xx+1]-st[yy][xx])
                    +alpha_y*((st[yy+1][xx]-st[yy][xx]) + alpha_x
                    *(st[yy+1][xx+1]+st[yy][xx]-st[yy][xx+1]-st[yy+1][xx])));
          pt[m][i] += rot[i][j];
        }
      }
    }

    return result;
  }

  // shift and rotate the matrix in Cartesian coordinate system
  Matrix<TT> getShiftRotate( int shift_x, 
                             int shift_y,
                             double rotateAngle,        // in degrees
                             int intp = 1               // 0 for the nearest neibooughboring, 1 for linear interpolation
                           )
  {
    if( shift_y <= -m_rows || shift_y >= m_rows || 
        shift_x <= -m_cols || shift_x >= m_cols    )
        return Matrix<TT>();

    int i, j;
    TT min_value = this->Min();

    // Perform the shift
    Matrix<TT> shift(m_rows, m_cols);
    TT** st = shift.data();
    for( i=0; i<m_rows; i++)
    {
      for(j=0; j<m_cols; j++)
      {
        if( (i>=shift_y) && (i-shift_y)<m_rows && 
            (j>=shift_x) && (j-shift_x)<m_cols )
          st[i][j] = m_data[i-shift_y][j-shift_x];
        else
          st[i][j] = min_value;
      }
    }

    // Check if rotation is necessary
    double arc = rotateAngle*pi/180;
    if( fabs(sin(arc))*::Max(m_rows, m_cols) < 1 )
      return shift;

    // Perform rotation based on the shifted image
    Matrix<TT> rotated(m_rows, m_cols);
    rotated = min_value;
    TT ** rot = rotated.data();

    double cos_, sin_;
    cos_ = cos(-arc);
    sin_ = sin(-arc);

    CPPoint original, transformed;
    int xx, yy;
    double alpha_x, alpha_y;

    double half_row = m_rows/2.;
    double half_col = m_cols/2.;

    double aa, bb;

    for( i = 0; i < m_rows; i++)
    {
      for( j = 0; j < m_cols; j++)
      {
        transformed.x = j - half_col + 0.5;
        transformed.y = i - half_row + 0.5;
        original.x = transformed.x*cos_ - transformed.y*sin_;
        original.y = transformed.x*sin_ + transformed.y*cos_;

        switch( intp )
        {
        case 0:
          xx = Nearest(original.x + half_col - 0.5);
          yy = Nearest(original.y + half_row - 0.5);
          if( yy <= m_rows-1 && yy >=0 && xx <= m_cols-1 && xx >= 0)
            rot[i][j] = st[yy][xx];
          break;
        case 1:
          aa = original.x + half_col - 0.5;
          xx = (int)floor(aa);
          alpha_x = aa - xx;

          bb = original.y + half_row - 0.5;
          yy = (int)floor(bb);
          alpha_y = bb - yy;

          if( yy < m_rows-1 && yy >=0 && xx < m_cols-1 && xx >= 0)
            rot[i][j] = TT(st[yy][xx] + alpha_x*(st[yy][xx+1]-st[yy][xx])
                    +alpha_y*((st[yy+1][xx]-st[yy][xx]) + alpha_x
                    *(st[yy+1][xx+1]+st[yy][xx]-st[yy][xx+1]-st[yy+1][xx])));
          break;
        default:
          break;
        }
      }
    }

    return rotated;
  }

protected:
  TT**  m_data;
  int   m_rows;
  int   m_cols;

private:
  void init(int rows, int cols, TT **newData = 0)
  {
    Vector<TT>::resize(rows*cols);
    if(rows <= 0 || cols <= 0 || buffer() == 0) { empty(); return; }
    m_data = new TT*[rows];
    if(m_data == 0) { empty(); return; }

    m_rows = rows;
    m_cols = cols;
	int i;
    for( i = 0; i < m_rows; i++ ) m_data[i] = m_buffer + m_cols*i;

    if( newData != 0 )
      for( i = 0; i < m_rows; i++ ) memcpy( m_data[i], newData[i], m_cols * sizeof( TT ) );
  }

  void init(int rows, int cols, TT *newData)
  {
    Vector<TT>::resize(rows*cols, newData);
    if(rows <= 0 || cols <= 0 || buffer() == 0) { empty(); return; }
    m_data = new TT*[rows];
    if(m_data == 0) { empty(); return; }

    m_rows = rows;
    m_cols = cols;
    for( int i = 0; i < m_rows; i++ ) m_data[i] = m_buffer + m_cols*i;
  }
};

// not logic operation
template<class TT>
Matrix<TT> operator~ (const Matrix<TT> & ve)
{
  int length = ve.length();
  Matrix<TT> temp(ve.rows(), ve.cols());
  TT * te1 = ve.buffer();
  TT * te2 = temp.buffer();
  int i = 0;

  while( i++ < length )
  {
    if(*te1) *te2 = 0;
    else *te2 = 1;
    te1++;
    te2++;
  }
  return temp;
}

///////////////////////////////////////////////////////////////////////
//
// Matrix arithmetic operation
//
///////////////////////////////////////////////////////////////////////

// Add two matrices element-by-element on the intersection part
template<class TT>
Matrix<TT> operator+ (const Matrix<TT> & m1, const Matrix<TT> & m2)
{
  int row = (m1.rows() <= m2.rows()) ? m1.rows() : m2.rows();
  int col = (m1.cols() <= m2.cols()) ? m1.cols() : m2.cols();

  Matrix<TT> m3; // initilize an empty matrix
    
  if( row <= 0 || col <= 0 ) return m3;
    
  m3.resize( row, col ); // reset the size of matrix

  // using temporal pointers to speed up the access to element
  TT ** temp1 = m1.data();
  TT ** temp2 = m2.data();
  TT ** temp3 = m3.data();

  for( int i = 0; i < row; i++ )
    for( int j = 0; j < col; j++ )
      temp3[i][j] = temp1[i][j] + temp2[i][j];

  return m3;
}

// Add a constant to a matrix
template<class TT>
Matrix<TT> operator+ (const Matrix<TT> & m1, const TT m2)
{
  int row = m1.rows();
  int col = m1.cols();
    
  Matrix<TT> m3( row, col );

  TT * temp3 = m3.buffer();
  TT * temp1 = m1.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp3 = *temp1 + m2;

    temp3++;
    temp1++;
    i++;
  }

  return m3;
}

template<class TT>
Matrix<TT> operator+ (const TT m2, const Matrix<TT> & m1)
{
  int row = m1.rows();
  int col = m1.cols();
    
  Matrix<TT> m3( row, col );

  TT * temp3 = m3.buffer();
  TT * temp1 = m1.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp3 = *temp1 + m2;

    temp3++;
    temp1++;
    i++;
  }
  return m3;
}

// minus of a matrix
template <class TT>
Matrix<TT> operator- (const Matrix<TT> & b)
{
    
  Matrix<TT> mRet( b.rows(), b.cols() );
  TT *temp = mRet.buffer();
  TT *temp1 = b.buffer();
  int size = b.rows()*b.cols();

  int i = 0;
  while( i < size)
  {
    *temp = -(*temp1);

    temp++;
    temp1++;
    i++;
  }

  return mRet;
}

// Substract a matrix from another on the intersection part
template<class TT>
Matrix<TT> operator- (const Matrix<TT> & m1, const Matrix<TT> & m2)
{
  int row = (m1.rows() <= m2.rows()) ? m1.rows() : m2.rows();
  int col = (m1.cols() <= m2.cols()) ? m1.cols() : m2.cols();

  Matrix<TT> m3; // initilize an empty matrix
    
  if( row <= 0 || col <= 0 ) return m3;
    
  m3.resize( row, col ); // reset the size of matrix

  // using temporal pointers to speed up the access to element
  TT ** temp1 = m1.data();
  TT ** temp2 = m2.data();
  TT ** temp3 = m3.data();

  for( int i = 0; i < row; i++ )
    for( int j = 0; j < col; j++ )
      temp3[i][j] = temp1[i][j] - temp2[i][j];

  return m3;
}

template<class TT>
Matrix<TT> operator- (const Matrix<TT> & m1, const TT m2)
{
  int row = m1.rows();
  int col = m1.cols();
    
  Matrix<TT> m3( row, col );

  TT * temp3 = m3.buffer();
  TT * temp1 = m1.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp3 = *temp1 - m2;

    temp3++;
    temp1++;
    i++;
  }

  return m3;
}

template<class TT>
Matrix<TT> operator- ( const TT m2, const Matrix<TT> & m1)
{
  int row = m1.rows();
  int col = m1.cols();
    
  Matrix<TT> m3( row, col );

  TT * temp3 = m3.buffer();
  TT * temp1 = m1.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp3 = m2 - *temp1;

    temp3++;
    temp1++;
    i++;
  }

  return m3;
}

// Mutiply two matrices
template<class TT>
Matrix<TT> operator* (const Matrix<TT> & m1, const Matrix<TT> & m2)
{
  int rows = m1.rows();
  int cols = m2.cols();

  if( rows == 0 || cols == 0 || m1.cols() != m2.rows() ) return m1;

  TT tsum;
  TT** t1 = m1.data();
  TT** t2 = m2.data();

  Matrix<TT> mRet( rows, cols );
  TT ** temp = mRet.data();

  for( int i = 0; i < rows; i++ )
  {
    for( int j = 0; j < cols; j++ )
    {
      tsum = 0;
      for( int k = 0; k < m1.cols(); k++ )
        tsum += t1[i][k] * t2[k][j];
      temp[i][j] = tsum;
    }
  }
  return mRet;
}

// Mutiply two vectors
template<class TT>
Matrix<TT> operator* (const Vector<TT> & v1, const Vector<TT> & v2)
{
  int rows = v1.length();
  int cols = v2.length();

  TT* t1 = v1.buffer();
  TT* t2 = v2.buffer();

  Matrix<TT> mRet( rows, cols );
  TT ** temp = mRet.data();

  for( int i = 0; i < rows; i++ )
  {
    for( int j = 0; j < cols; j++ )
    {
      temp[i][j] = t1[i]*t2[j];
    }
  }
  return mRet;
}

/////////////////////////////////////////////////////////////////
// Mutiply matrix by a constant
/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
template <class TT>
Matrix<TT> operator* (const TT c, const Matrix<TT> & ma)
{
  int row = ma.rows();
  int col = ma.cols();
    
  Matrix<TT> mRet( row, col );
  TT * temp = mRet.buffer();
  TT * temp1 = ma.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp = (c*(*temp1));

    temp++;
    temp1++;
    i++;
  }

  return mRet;
}

template <class TT>
Matrix<TT> operator* (const Matrix<TT> & ma, const TT c)
{
  int row = ma.rows();
  int col = ma.cols();
    
  Matrix<TT> mRet( row, col );
  TT * temp = mRet.buffer();
  TT * temp1 = ma.buffer();

  int i = 0;
  int size = row*col;

  while( i < size )
  {
    *temp = (c*(*temp1));

    temp++;
    temp1++;
    i++;
  }

  return mRet;
}

template <class TT>
Matrix<TT> operator/ (const Matrix<TT> & ma, const TT c)
{
  int row = ma.rows();
  int col = ma.cols();

  if( c < 1.e-20 && c > -1.e-20)
  {
    return ma;
  }

  Matrix<TT> mRet( row, col );
  TT * temp = mRet.buffer();
  TT * temp1 = ma.buffer();

  int size = row*col;
  int i = 0;

  while( i < size )
  {
    *temp = TT( *temp1 / c);

    temp++;
    temp1++;
    i++;
  }

  return mRet;
}

//////////////////////////////////////////////////////////
//  multiplication between vector and matrix
//
template <class TT>
Vector<TT> operator* ( const Vector<TT> & vt, const Matrix<TT> & ma)
{
  if( vt.length()<=0 || ma.rows()<=0 || ma.cols()<=0 ) { return Vector<TT>(); }
  if( vt.length() != ma.rows() ) { return Vector<TT>(); }

  Vector<TT> vRet(ma.cols());
  TT * temp = vRet.data();

  TT * temp1 = vt.data();
  TT ** temp2 = ma.data();

  int row = ma.rows();
  int col = ma.cols();

  for( int i=0; i< col; i++)
    for( int j=0; j<row; j++)
      temp[i] += (temp1[j]*temp2[j][i]);

  return vRet;
}

template <class TT>
Vector<TT> operator* ( const Matrix<TT> & ma, const Vector<TT> & vt)
{
  if( vt.length()<=0 || ma.rows()<=0 || ma.cols()<=0 ) { return Vector<TT>(); }
  if( vt.length() != ma.cols() ) { return Vector<TT>(); }

  Vector<TT> vRet(ma.rows());
  TT * temp = vRet.data();

  TT * temp1 = vt.data();
  TT ** temp2 = ma.data();

  int row = ma.rows();
  int col = ma.cols();

  for( int i=0; i< row; i++)
    for( int j=0; j<col; j++)
      temp[i] += (temp2[i][j]*temp1[j]);

  return vRet;
}

// if you are tired of using <...>, try this
typedef Matrix<char>            CHARMat;
typedef Matrix<short>           SHORTMat;
typedef Matrix<int>             INTMat;
typedef Matrix<long>            LONGMat;
typedef Matrix<float>           FLOATMat;
typedef Matrix<double>          DOUBLEMat;
typedef Matrix<unsigned char>   UCHARMat;
typedef Matrix<unsigned short>  USHORTMat;
typedef Matrix<unsigned int>    UINTMat;
typedef Matrix<unsigned long>   ULONGMat;

#endif //_MATRIX_OF_NUMBER_TYPE_H
