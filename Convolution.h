// convolution.h

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
 * The general window averaging is implemented for one- and two-dimensional data.
 * The stochastic Savitzky-Golay filter is listed here.
 */

#ifndef _CUBICIMAGING_CONVOLUTION_OPERATION_H_
#define _CUBICIMAGING_CONVOLUTION_OPERATION_H_

#include "matrix.h"

//! 5-point Savitzky-Golay filter
static const double SavGol2[5] = { -0.086, 0.343, 0.486, 0.343, -0.086 };

//! 9-point Savitzky-Golay filter
static const double SavGol4[9] = { 0.035, -0.128,  0.070, 0.315, 0.417, 
                                   0.315,  0.070, -0.128, 0.035 };
//! 11-point Savitzky-Golay filter
static const double SavGol5[11] = { 0.042, -0.105, -0.023,  0.140, 0.280, 0.333,
                                    0.280,  0.140, -0.023, -0.105, 0.042 };

//! one-dimensional Gaussian filter generator
void make_gaussian_kernel(const double sigma, DOUBLEVec & kernel);

//! two-dimensional Hanning filter
static const double HanningKernel[25] = {
	0.0069444444,    0.0208333333,    0.0277777778,    0.0208333333,    0.0069444444,
	0.0208333333,    0.0625000000,    0.0833333333,    0.0625000000,    0.0208333333,
	0.0277777778,    0.0833333333,    0.1111111111,    0.0833333333,    0.0277777778,
	0.0208333333,    0.0625000000,    0.0833333333,    0.0625000000,    0.0208333333,
	0.0069444444,    0.0208333333,    0.0277777778,    0.0208333333,    0.0069444444	};

//! two-dimensional Gaussian filter generator
void make_gaussian_kernel(const double sigma, DOUBLEMat & kernel);

//! 5-point smoothing
template<class TT>
void Smooth53(Vector<TT> & in)
{
  if(in.length() < 5)
    return;

  TT* da = in.data();
  int dim = in.length();
  Vector<TT> out(dim);
  TT* da1 = out.data();

  da1[0] = (69*da[0] + 4*(da[1]+da[3]) - 6*da[2] - da[4])/70;
  da1[1] = (2*(da[0]+da[4]) + 27*da[1] + 12*da[2] - 8*da[3])/35;
  da1[dim-2] = (2*(da[dim-1]+da[dim-5]) + 27*da[dim-2] + 12*da[dim-3] - 8*da[dim-4])/35;
  da1[dim-1] = (69*da[dim-1] + 4*(da[dim-2]+da[dim-4]) - 6*da[dim-3] - da[dim-5])/70;
  for(int i=2; i<dim-2; i++)
  {
    da1[i] = (17*da[i] + 12*(da[i-1]+da[i+1]) - 3*(da[i-2]+da[i+2]))/35;
  }

  in = out;
}

//! one-dimensional convolution with predefined filter
template<class TT>
void Convolution(Vector<TT> & data, const Vector<double> & kernel, bool symmetry = true)
{
	register int len = data.length();
	if( len <= 0 ) return;
	register int len_ = kernel.length(); // must be an odd number
	if( !(len_%2) ) return;
	register int len_k = kernel.length()/2;
	register int i, j;

	Vector<TT> temp(len+2*len_k);
	TT *te = temp.data();
	TT *da = data.data();
	double * ker = kernel.data();

	for( i=0; i<len_k; i++)
		te[i] = da[0];
	for( i=len; i<len+len_k; i++)
		te[i] = da[len-1];
	memcpy(te+len_k, da, sizeof(TT)*len);

	register double value;
  register int len_1 = len_-1;

  if( symmetry )
  {
	  for(i=0; i<len; i++)
	  {
		  value = 0;
		  for(j=0; j<len_k; j++)
			  value += ker[j]*(te[i+j]+te[i+len_1-j]);
      value += ker[len_k]*te[i+len_k];
		  data[i] = TT(value);
	  }
  }
  else
  {
	  for(i=0; i<len; i++)
	  {
		  value = 0;
		  for(j=0; j<len_; j++)
			  value += ker[j]*te[i+j];
		  data[i] = TT(value);
	  }
  }
}

//! one-dimensional convolution with predefined filter
template<class TT>
void Convolution(Vector<TT> & data, Vector<int> & kernel, bool symmetry = true)
{
	register int len = data.length();
	if( len <= 0 ) 
    return;
	register int len_ = kernel.length(); // must be an odd number
  register int sum = int(kernel.sum());
	if( !(len_%2) || sum <= 0 ) 
    return;
	register int len_k = kernel.length()/2;
	register int i, j;

	Vector<TT> temp(len+2*len_k);
	TT *te = temp.data();
	TT *da = data.data();
	int * ker = kernel.data();

	for( i=0; i<len_k; i++)
		te[i] = da[0];
	for( i=len; i<len+len_k; i++)
		te[i] = da[len-1];
	memcpy(te+len_k, da, sizeof(TT)*len);

	register int value;
  register int len_1 = len_-1;

  if( symmetry )
  {
	  for(i=0; i<len; i++)
	  {
		  value = 0;
		  for(j=0; j<len_k; j++)
			  value += ker[j]*(te[i+j]+te[i+len_1-j]);
      value += ker[len_k]*te[i+len_k];
		  data[i] = TT(value/sum);
	  }
  }
  else
  {
	  for(i=0; i<len; i++)
	  {
		  value = 0;
		  for(j=0; j<len_; j++)
			  value += ker[j]*te[i+j];
		  data[i] = TT(value/sum);
	  }
  }
}

//! two-dimensional convolution with predefined filter
template<class TT>
void Convolution(Matrix<TT> & da, const Matrix<double> & ker, bool symmetry = true)
{
	register int row = da.rows();
	register int col = da.cols();
	if( row <= 0 || col <= 0 ) 
    return;
	register int row_ = ker.rows(); // must be an odd number
	register int col_ = ker.cols(); // must be an odd number
	if( !(row_%2 && col_%2) ) 
    return;

	register int row_k = ker.rows()/2;
	register int col_k = ker.cols()/2;

	int i, j, m, n;

	Matrix<TT> te(row+2*row_k, col+2*col_k);
	TT** temp = te.data();
	TT** data = da.data();
	double** kernel = ker.data();

	for( i=0; i<row; i++)
    memcpy(temp[i+row_k]+col_k, data[i], sizeof(TT)*col);

  for( i=row_k; i<row+row_k; i++)
	{
		for( j=0; j<col_k; j++)
			temp[i][j] = data[i-row_k][0];
		for( j=col; j<col+col_k; j++)
			temp[i][j] = data[i-row_k][col-1];
	}
	for( i=0; i<row_k; i++)
		for( j=0; j<col+2*col_k; j++)
			temp[i][j] = temp[row_k][j];
	for( i=row+row_k; i<row+2*row_k; i++)
		for( j=0; j<col+2*col_k; j++)
			temp[i][j] = temp[row+row_k-1][j];

	register double value;

  register int row_1 = row_ - 1;
  register int col_1 = col_ - 1;

  if( symmetry )
  {
	  for( i=0; i<row; i++)
    {
		  for( j=0; j<col; j++)
		  {
			  value = 0;
			  for(m=0; m<row_k; m++)
        {
				  for(n=0; n<col_k; n++)
					  value += kernel[m][n]*(temp[i+m][j+n]+temp[i+row_1-m][j+n]+temp[i+m][j+col_1-n]+temp[i+row_1-m][j+col_1-n]);
        }
        for(n=0; n<col_k; n++)
          value += kernel[row_k][n]*(temp[i+row_k][j+n]+temp[i+row_k][j+col_1-n]);
        for(m=0; m<row_k; m++)
          value += kernel[m][col_k]*(temp[i+m][j+col_k]+temp[i+row_1-m][j+col_k]);
        value += kernel[row_k][col_k]*temp[i+row_k][j+col_k];

			  data[i][j] = TT(value);
		  }
    }
  }
  else
  {
    for( i=0; i<row; i++)
    {
		  for( j=0; j<col; j++)
		  {
			  value = 0;
			  for(m=0; m<row_; m++)
        {
				  for(n=0; n<col_; n++)
					  value += kernel[m][n]*temp[i+m][j+n];
        }
			  data[i][j] = TT(value);
		  }
    }
  }
}

//! two-dimensional convolution with predefined filter
template<class TT>
void Convolution(Matrix<TT> & da, Matrix<int> & ker, bool symmetry = true)
{
	register int row = da.rows();
	register int col = da.cols();
	if( row <= 0 || col <= 0 ) 
    return;
	register int row_ = ker.rows(); // must be an odd number
	register int col_ = ker.cols(); // must be an odd number
  register int sum = ker.sum();
	if( !(row_%2 && col_%2)  || sum <= 0 ) 
    return;

	register int row_k = ker.rows()/2;
	register int col_k = ker.cols()/2;

	register int i, j, m, n;

	Matrix<TT> te(row+2*row_k, col+2*col_k);
	TT** temp = te.data();
	TT** data = da.data();
	int** kernel = ker.data();

	for( i=0; i<row; i++)
    memcpy(temp[i+row_k]+col_k, data[i], sizeof(TT)*col);

	for( i=row_k; i<row+row_k; i++)
	{
		for( j=0; j<col_k; j++)
			temp[i][j] = data[i-row_k][0];
		for( j=col; j<col+col_k; j++)
			temp[i][j] = data[i-row_k][col-1];
	}
	for( i=0; i<row_k; i++)
		for( j=0; j<col+2*col_k; j++)
			temp[i][j] = temp[row_k][j];
	for( i=row+row_k; i<row+2*row_k; i++)
		for( j=0; j<col+2*col_k; j++)
			temp[i][j] = temp[row+row_k-1][j];

	register int value;
  register int row_1 = row_ - 1;
  register int col_1 = col_ - 1;

  if( symmetry )
  {
	  for( i=0; i<row; i++)
    {
		  for( j=0; j<col; j++)
		  {
			  value = 0;
			  for(m=0; m<row_k; m++)
        {
				  for(n=0; n<col_k; n++)
					  value += kernel[m][n]*(temp[i+m][j+n]+temp[i+row_1-m][j+n]+temp[i+m][j+col_1-n]+temp[i+row_1-m][j+col_1-n]);
        }
        for(n=0; n<col_k; n++)
          value += kernel[row_k][n]*(temp[i+row_k][j+n]+temp[i+row_k][j+col_1-n]);
        for(m=0; m<row_k; m++)
          value += kernel[m][col_k]*(temp[i+m][j+col_k]+temp[i+row_1-m][j+col_k]);
        value += kernel[row_k][col_k]*temp[i+row_k][j+col_k];

			  data[i][j] = TT(value/sum);
		  }
    }
  }
  else
  {
	  for( i=0; i<row; i++)
    {
		  for( j=0; j<col; j++)
		  {
			  value = 0;
			  for(m=0; m<row_; m++)
        {
				  for(n=0; n<col_; n++)
					  value += kernel[m][n]*temp[i+m][j+n];
        }
			  data[i][j] = TT(value/sum);
		  }
    }
  }
}

#endif	_APPIMGTECH_CONVOLUTION_OPERATION_H