
#include "specialFunc.h"
#include "utilfunc.h"

double GammaDeviates(int evt, int *seed)
{
	if (evt < 1) return 0.0;

	double am, e, s, v1, v2, x, y;
  int j;

	if (evt < 6) 
  {
		x = 1.0;
		for (j=1; j<=evt; j++) x *= Rand(seed);
		x = -log(x);
	} 
  else 
  {
		while( true ) 
    {
			while( true ) 
      {
				while( true ) 
        {
					v1 = Rand(seed);
					v2 = 2.0*Rand(seed)-1.0;
          if(v1*v1+v2*v2 <= 1.0) break;
				}
				y = v2/v1;
				am = evt-1;
				s = sqrt(2.0*am+1.0);
				x = s*y+am;
        if(x > 0.0) break;
			}
			e = (1.0+y*y)*exp(am*log(x/am)-s*y);
      if(Rand(seed) <= e) break;
		}
	}
	return x;
}

int PoissonDeviates(const double mean, int * seed)
{
  if( mean <= 1.e-8 ) return 0;

	static double sq, almean, g, oldm = -1.0;
	double t, em, y;
  int poisson;

	if(mean < 12.0) 
  {
		if(mean != oldm) 
    {
			oldm = mean;
			g = exp(-mean);
		}
		poisson = 0;
		t = Rand(seed); 
		while(t > g) 
    {
			poisson ++;
			t *= Rand(seed);
		}
	} 
  else 
  {
		if(mean != oldm) 
    {
			oldm = mean;
			sq = sqrt(2.0*mean);
			almean = log(mean);
			g = mean*almean - GammaLn56(mean+1.0);
		}
		while( true ) 
    {
			while( true )
      {
				y = tan(pi*Rand(seed));
				em = sq*y + mean;
        if(em >= 0.0) break;
			}
			poisson = int(em);
			t = 0.9*(1.0+y*y)*exp(poisson*almean-GammaLn56(poisson+1.0)-g);
      if(Rand(seed) <= t) break;
		}
	}

	return poisson;
}

double GaussDeviates(int * seed)
{
	static bool switcher = true;
	static double gaussValue;
	double fac, v1, v2, rsq = 2.;

	if( switcher ) 
  {
		while(rsq >= 1.0 || rsq < eps)
    {
			v1 = 2.0*Rand(seed) - 1.0;
			v2 = 2.0*Rand(seed) - 1.0;
			rsq = v1*v1 + v2*v2;
		} 

    fac = sqrt(-2.0*log(rsq)/rsq);
		gaussValue = v1*fac;
		switcher = !switcher;
		return (v2*fac);
	} 
  else 
  {
		switcher = !switcher;
		return gaussValue;
	}
}

double LegendrePoly(int l, int m, const double x)
{
	double fact, pll, pmm, pmmp1, somx2;
	int i, ll;

	if (m < 0 || m > l || Abs(x) > 1.0) return 1.e100;
	
  pmm = 1.0;
	if(m > 0) 
  {
		somx2 = sqrt((1.0-x)*(1.0+x));
		fact = 1.0;
		for(i=1; i<=m; i++)
    {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if(l > m)
  {
		pmmp1 = x*(2*m+1)*pmm;
		if(l == (m+1)) return pmmp1;
		else 
    {
			for(ll=m+2; ll<=l; ll++) 
      {
				pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm = pmmp1;
				pmmp1 = pll;
			}
			return pll;
		}
	}
	else return pmm;
}

double  BetaIncomplete(double a, double b, double x)
{
	double bt;

	if (x < 0.0 || x > 1.0) return -1e30;;
	if (x == 0.0 || x == 1.0) return 0.0;

	bt = exp(GammaLn56(a+b)-GammaLn56(a)-GammaLn56(b)+a*log(x)+b*log(1.0-x));

	if(x < (a+1.0)/(a+b+2.0)) return bt*BetaFraction(a,b,x)/a;
	else return 1.0-bt*BetaFraction(b,a,1.0-x)/b;
}

/////////////////////////////////////////////////////////////////////////////////////
// utility functions
//
double GammaP(const double a, const double x)
{
	if(x < 0.0 || a <= 0.0) return -1.;

	if(x < (a+1.0)) return GammaSeriesExt(a, x);
  else return (1.0 - GammaFractionExt(a, x));
}

double GammaQ(const double a, const double x)
{
	if(x < 0.0 || a <= 0.0)  return -1.;

	if(x < (a+1.0)) return (1. - GammaSeriesExt(a, x));
  else return GammaFractionExt(a, x);
}

double GammaSeriesExt(double a, double x)
{
	int iter, itmax = 100;
	double sum, del, ap, gln, eps = 1.e-8;
	gln = GammaLn56(a);

  if(x > 0.0 && a > 0.0) 
  {
		ap = a;
		del = sum = 1.0/a;
    iter = 0;
		while(iter < itmax) 
    {
			++ap;
			del *= x/ap;
			sum += del;
			if (Abs(del) < Abs(sum)*eps) return (sum*exp(-x+a*log(x)-(gln)));
      iter ++;
		}
	}

  return 0.0;
}

double GammaFractionExt(double a, double x)
{
  if(x <= 0.0 || a <= 0.0) return 0.0;

	int iter, itmax = 100;
	double an, b, c, d, del, h, gln;
  double fpmin = 1.e-30, eps = 1.e-8;
	gln = GammaLn56(a);

  b = x+1.0-a;
	c = 1.0/fpmin;
	d = 1.0/b;
	h = d;
  iter = 1;
	while(iter <= itmax) 
  {
		an = -iter*(iter-a);
		b += 2.0;
		d = an*d+b;
		if(Abs(d) < fpmin) d = fpmin;
		c = b+an/c;
		if(Abs(c) < fpmin) c = fpmin;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if(Abs(del-1.0) < eps) return (exp(-x+a*log(x)-gln)*h);

    iter ++;
	}

  return 0.0;
}

double  BetaFraction(double a, double b, double x)
{
	int m, m2;
	double aa, c, d, del, h, qab, qam, qap;
  const int itmax = 100;
  const double eps = 1.0e-8;
  const double fpmin = 1.0e-30;

	qab = a+b;
	qap = a+1.0;
	qam = a-1.0;
	c = 1.0;
	d = 1.0-qab*x/qap;
	if(Abs(d) < fpmin) d = fpmin;
	d = 1.0/d;
	h = d;

  m = 1;
	while(m <= itmax) 
  {
		m2 = 2*m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0+aa*d;
		if(Abs(d) < fpmin) d = fpmin;

    c = 1.0+aa/c;
		if(Abs(c) < fpmin) c = fpmin;

		d = 1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0+aa*d;
		if(Abs(d) < fpmin) d = fpmin;

		c = 1.0+aa/c;
		if(Abs(c) < fpmin) c = fpmin;

		d = 1.0/d;
		del = d*c;
		h *= del;

		if(Abs(del-1.0) < eps) break;
    m++;
	}

  return h;
}