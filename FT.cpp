
#include <math.h>
#include "FT.h"

void FFT(double *x, double *y, int n, int dir)
{
   long m, i, i1, j, k, i2, l, l1, l2;
   double c1, c2, tx, ty, t1, t2, u1, u2, z;

   // Calculate the exponentail number
   m = 0;
   while(n/(2<<m)){ m++; }
   if(n%(1<<m) != 0) return;

   // Do the bit reversal
   i2 = n >> 1;
   j = 0;
   for (i=0; i<n-1; i++) 
   {
      if (i < j) 
      {
         tx = x[i]; ty = y[i];
         x[i] = x[j]; y[i] = y[j];
         x[j] = tx; y[j] = ty;
      }
      k = i2;
      while (k <= j) 
      {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   // Compute the FFT
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0; l<m; l++) 
   {
      l1 = l2; l2 <<= 1;
      u1 = 1.0; u2 = 0.0;
      for (j=0; j<l1; j++) 
      {
         for (i=j; i<n; i+=l2) 
         {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   // Normalizing transform
  double a = sqrt(1.0/n);
  for (i=0; i<n; i++) { x[i] *= a; y[i] *= a; }
}

void DFT(double *x, double *y, int n, int dir)
{
   double *mem = new double[4*n]; if(mem == 0) return;
   double *x2 = mem;
   double *y2 = mem + n;
   double *cs = mem + n*2;
   double *sn = mem + n*3;

   long i, k, m;
   double arg= -dir*2.0*3.141592653589793;
   double cosarg = cos(arg/n);
   double sinarg = sin(arg/n);
   cs[0] = 1.0; sn[0] = 0.0;
   for (i=1; i<n; i++)
   {
     cs[i] = cs[i-1]*cosarg - sn[i-1]*sinarg;
     sn[i] = cs[i-1]*sinarg + sn[i-1]*cosarg;
   }

   for (i=0; i<n; i++) 
   {
      x2[i] = 0; y2[i] = 0;
      for (k=0; k<n; k++) 
      {
        m = (i*k)%n;
        x2[i] += (x[k] * cs[m] - y[k] * sn[m]);
        y2[i] += (x[k] * sn[m] + y[k] * cs[m]);
      }
   }

   // Copy back the normalized data
   double a = sqrt(1.0/n);
    for (i=0; i<n; i++) 
    {
        x[i] = x2[i] * a;
        y[i] = y2[i] * a;
    }

    delete[] mem;
}

void DST1(double *x, int n)
{
   double *mem = new double[3*n]; if(mem == 0) return;
   double *y = mem;
   double *sn = mem + n*1;
   double *cs = mem + n*2;

   int i, k, m, r;
   double cosarg = cos(3.141592653589793/n);
   double sinarg = sin(3.141592653589793/n);
   cs[0] = 1.0; sn[0] = 0.0;
   for (i=1; i<n; i++)
   {
     cs[i] = cs[i-1]*cosarg - sn[i-1]*sinarg;
     sn[i] = cs[i-1]*sinarg + sn[i-1]*cosarg;
   }

   for (i=0; i<n; i++) 
   {
      y[i] = 0;
      for (k=0; k<n; k++) 
      {
        m = i*k; r = m%n; m = m/n;
        if(m%2 == 0) 
          y[i] += x[k] * sn[r];
        else 
          y[i] -= x[k] * sn[r];
      }
   }

   // Copy back the normalized data
   double a = sqrt(2.0/n);
   for (i=0; i<n; i++) { x[i] = y[i]*a; }
   
   delete[] mem;
}

void DCT3(double *x, int n, int dir)
{
   double *mem = new double[5*n]; if(mem == 0) return;
   double *y = mem;
   double *cs = mem + n*1;
   double *sn = mem + n*2;
   double *ch = mem + n*3;
   double *sh = mem + n*4;

   int i, k, m, r;
   double cosarg, sinarg, sqrt2 = sqrt(2.0); 
   cosarg = cos(3.141592653589793/n);
   sinarg = sin(3.141592653589793/n);
   cs[0] = 1.0; sn[0] = 0.0;
   for (i=1; i<n; i++)
   {
     cs[i] = cs[i-1]*cosarg - sn[i-1]*sinarg;
     sn[i] = cs[i-1]*sinarg + sn[i-1]*cosarg;
   }
   cosarg = cos(0.5*3.141592653589793/n);
   sinarg = sin(0.5*3.141592653589793/n);
   ch[0] = 1.0; sh[0] = 0.0;
   for (i=1; i<n; i++)
   {
     ch[i] = ch[i-1]*cosarg - sh[i-1]*sinarg;
     sh[i] = ch[i-1]*sinarg + sh[i-1]*cosarg;
   }

   if(dir == -1)
   {
     for (i=0; i<n; i++) 
     {
        y[i] = x[0]/sqrt2;
        for (k=1; k<n; k++) 
        {
          m = i*k; r = m%n; m = m/n;
          if(m%2 == 0) 
            y[i] += (x[k]*cs[r]*ch[k] - x[k]*sn[r]*sh[k]);
          else 
            y[i] -= (x[k]*cs[r]*ch[k] - x[k]*sn[r]*sh[k]);
        }
     }
   }
   else
   {
     y[0] = 0.0;
     for (k=0; k<n; k++) y[0] += x[k];
     y[0] /= sqrt2;
     for (i=1; i<n; i++) 
     {
       y[i] = 0.0;
       for (k=0; k<n; k++) 
       {
         m = i*k; r = m%n; m = m/n;
         if(m%2 == 0) 
           y[i] += (x[k]*cs[r]*ch[i] - x[k]*sn[r]*sh[i]);
         else 
           y[i] -= (x[k]*cs[r]*ch[i] - x[k]*sn[r]*sh[i]);
       }
     }
   }

   // Copy back the normalized data
   double a = sqrt(2.0/n);
   for (i=0; i<n; i++) { x[i] = y[i]*a; }

   delete[] mem;
}
