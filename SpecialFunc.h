// specialFunc.h

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
 * \date Created at 12/26/2000

 * This file collects a lot of commonly used special functions for statistic analysis.
 */

#ifndef _CUBICIMAGING_SPECIAL_FUNCTION_H_
#define _CUBICIMAGING_SPECIAL_FUNCTION_H_

//////////////////////////////////////////////////////
// commonly used special functions
double  GammaDeviates(int event, int *seed);
int     PoissonDeviates(const double mean, int * seed);
double  GaussDeviates(int * seed);
double  LegendrePoly(int l, int m, const double x);
double  BetaIncomplete(double a, double b, double x);

double  StudentT_Test(double chisqr, int fDegree);
double  F_Test(double var1, double var2, double fDegree);
double  Spearman(double* arr1, int dim1, double* arr2, int dim2);
double  KendallTau(double* arr1, int dim1, double* arr2, int dim2);
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// utility functions for calculating the special functions
double  GammaP(const double a, const double x);
double  GammaQ(const double a, const double x);
double  GammaSeriesExt(double a, double x);
double  GammaFractionExt(double a, double x);
double  BetaFraction(double a, double b, double x);

#endif
