
#ifndef _CUBICIMAGING_UTILITY_DISCRETE_TRANSFORM_H_
#define _CUBICIMAGING_UTILITY_DISCRETE_TRANSFORM_H_

void FFT(double *x, double *y, int n, int dir);

void DFT(double *x, double *y, int n, int dir);

void DST1(double *x, int n);

void DCT3(double *x, int n, int dir);

#endif