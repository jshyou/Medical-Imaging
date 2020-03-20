// convolution.cpp

#include "convolution.h"

void make_gaussian_kernel(double sigma, Vector<double> & ker)
{
   int i, center;
   double x, fx, sum = 0.0;

   int windowsize = (int)(1 + 2 * ceil(2.5 * sigma));
   center = (windowsize) / 2;

   double *kernel = new double[windowsize];

   for(i=0; i<windowsize; i++)
   {
     x = (double)(i - center);
     fx = exp(-0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
     kernel[i] = fx;
     sum += fx;
   }

   for(i=0; i<windowsize; i++)
   {
	   kernel[i] /= sum;
   }

   ker = Vector<double>(windowsize, kernel);

   delete [] kernel;
}

void make_gaussian_kernel(double sigma, Matrix<double> & ker)
{
   int i, j, center;
   double x, y, f, sum = 0.0;

   int windowsize = (int)(1 + 2 * ceil(2.5 * sigma));
   center = (windowsize) / 2;

   ker.resize(windowsize, windowsize);
   double ** kernel = ker.data();

   for(i=0; i<windowsize; i++)
   {
	   for(j=0; j<windowsize; j++)
	   {
		   x = (i - center);
		   y = (j - center);
		   f = exp(-0.5*(x*x+y*y)/(sigma*sigma)) / (sigma * sqrt(pi+pi));
		   kernel[i][j] = f;
		   sum += f;
	   }
   }

   for(i=0; i<windowsize*windowsize; i++)
   {
	   kernel[0][i] /= sum;
   }
}
