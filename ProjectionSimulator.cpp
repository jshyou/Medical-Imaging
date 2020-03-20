// ProjectionSimulator.cpp: implementation of the ProjectionSimulator class.
//
//////////////////////////////////////////////////////////////////////

#include "ProjectionSimulator.h"
#include <math.h>
#include <memory.h>

const double pi = 3.1415926535897932385;
const double rad = 48.5; // for normal simulation 91.5;// for noise analysis# 54.0; //for iterative FBP# 
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ProjectionSimulator::ProjectionSimulator()
{
}

ProjectionSimulator::~ProjectionSimulator()
{
}

double ProjectionSimulator::ellip(double x, double y, double x0, double y0, double phi, double a0, double b0)
{
	double value = 1.0;
	double s0, t0, s, t, thresh, a, b;

	a = cos(phi); b = sin(phi);
	s0 = x0*a + y0*b; t0 = -x0*b + y0*a;
	s = x*a + y*b; t = -x*b + y*a;

  thresh = (s-s0)*(s-s0)/(a0*a0) + (t-t0)*(t-t0)/(b0*b0);
	if(thresh > 1.0) value = 0.;

	return value;
}

double ProjectionSimulator::ellip(double x, double y, double z, double x0, double y0, double z0, double theta, double phi, double a0, double b0, double c0)
{
	double value = 1.0;
	double s0, t0, u0, s, t, u, thresh, a, b, c1, c2;

	a = cos(phi); b = sin(phi); 
  c1 = cos(theta); c2 = sin(theta);
	s0 = x0*a + y0*b; t0 = -x0*b + y0*a; u0 = z0;
	s = x*a + y*b; t = -x*b + y*a; u = z;

  thresh = (s-s0)*(s-s0)/(a0*a0) + (t-t0)*(t-t0)/(b0*b0) + (u-u0)*(u-u0)/(c0*c0);
	if(thresh > 1.0) value = 0.;

	return value;
}

double ProjectionSimulator::ellipProj(double x0, double y0, double phi, double a, double b, double l, double theta, double mu)
{
	double s0, t0, psi, a0, b0, c0, t1, t2, sqr, eradon;
	double eps = 1.0e-7;

	s0 =  x0*cos(phi) + y0*sin(phi);
	t0 = -x0*sin(phi) + y0*cos(phi);
	psi = pi/2 + theta - phi;
	a0 = sin(psi)*sin(psi)/(a*a) + cos(psi)*cos(psi)/(b*b);
	b0 = sin(psi)*(cos(psi)*l-s0)/(a*a) + cos(psi)*(t0-l*sin(psi))/(b*b);
	c0 = (s0-l*cos(psi))*(s0-l*cos(psi))/(a*a)+(t0-l*sin(psi))*(t0-l*sin(psi))/(b*b)-1.0;

  sqr = b0*b0 - a0*c0;
	if(sqr > eps)
	{
    sqr = sqrt(sqr);
		t1 = (b0 - sqr)/a0;
		t2 = (b0 + sqr)/a0;
		if(mu < eps)
			eradon = t2 - t1;
		else
			eradon = (exp(mu*t2) - exp(mu*t1))/mu;
	}
	else
		eradon = 0.0;

	return eradon;
}

void ProjectionSimulator::GenerateSLPhantom(double **SLImage, int dim)
{
  int m, n;
  double a, b;
	double HALFGRID = dim/2. - 0.5;
  double Delta = 1.0/HALFGRID;

  for(m=0; m<dim; m++)
	{
		for(n=0; n<dim; n++)
		{
			b = (m-HALFGRID)*Delta;
			a = (n-HALFGRID)*Delta;

      SLImage[m][n]
				 = 9.*ellip(a, b, 0.0, 0.0, 0.0, 35.0/rad, 46.0/rad)
				 - 4.*ellip(a, b, 0.0, -0.5/rad, 0.0, 33.0/rad, 43.50/rad)
				 +    ellip(a, b, 0.0, 17.5/rad, 0.0, 10.5/rad, 12.50/rad)
				 +    ellip(a, b, 0.0, -5.0/rad, 0.0, 2.5/rad, 2.5/rad)
				 + 2.*ellip(a, b, 0.0, 5.0/rad, 0.0, 2.5/rad, 2.5/rad)
				 +    ellip(a, b, 0.0, -30.0/rad, 0.0, 1.5/rad, 1.5/rad)
				 +    ellip(a, b, -6.0/rad, -30.0/rad, 0.0, 2.5/rad, 1.5/rad)
				 +    ellip(a, b, 5.0/rad, -30.0/rad, 0.0, 1.5/rad, 2.5/rad)
				 - 2.*ellip(a, b, 11.0/rad, 0.0, -18.0*pi/180.0, 5.5/rad, 15.5/rad)
				 - 2.*ellip(a, b, -11.0/rad, 0.0, 18.0*pi/180.0, 8.0/rad, 20.5/rad);
		}
	}
}

void ProjectionSimulator::GenerateSLPhantom(double ***SLVol, int dim, int slice)
{
  // to be implemented
}

void ProjectionSimulator::GenerateAttenMap(double ** attMap, int dim, bool unif)
{
  int m, n;
  double a, b;
	double HALFGRID = dim/2. - 0.5;
  double Delta = 1.0/HALFGRID;
  double step = RADIUS_IMAGE*2./dim;

  if( unif )
  {
	  for(m=0; m<dim; m++)
	  {
		  for(n=0; n<dim; n++)
		  {
			  b = (m-HALFGRID)*Delta;
			  a = (n-HALFGRID)*Delta;
			  attMap[m][n] = step*ellip(a, b, 0.0, 0.0, 0.0, 38.0/rad, 48.0/rad);
		  }
	  }
    return;
  }

	for(m=0; m<dim; m++)
	{
		for(n=0; n<dim; n++)
		{
			b = (m-HALFGRID)*Delta;
			a = (n-HALFGRID)*Delta;
			attMap[m][n] 
				 = 0.75*ellip(a, b, 0.0, 0.0, 0.0, 38.0/rad, 48.0/rad)
         - 0.50*ellip(a, b, 0.0,-21.5/rad, 0.0, 20./rad, 15.0/rad)
				 - 0.50*ellip(a, b, 0.0, 21.5/rad, 0.0, 20./rad, 15.0/rad)
         + 0.25*ellip(a, b,-25./rad, 0.0, 0.0, 5./rad, 5.0/rad)
         + 0.25*ellip(a, b, 25./rad, 0.0, 0.0, 5./rad, 5.0/rad);
      attMap[m][n] *= step;
		}
	}
}

void ProjectionSimulator::ProjectorAPar(double **Radon, int angle, int bin)
{
	int i, k;
	double s, theta, mu = 0.;
  const double pi = 3.141592654;
	double ANGLEGRID = pi/angle;

  for(k=0; k<angle; k++)
	{
		theta = k*ANGLEGRID;
		for(i=0; i<bin; i++)
		{				
			s = (i+i-bin+1.0)/bin;
			Radon[k][i] = 
				 9*ellipProj(0.0, 0.0, 0.0, 35.0/rad, 46.0/rad, s, theta, mu)
				-4*ellipProj(0.0, -0.5/rad, 0.0, 33.0/rad, 43.50/rad, s, theta, mu)
				+  ellipProj(0.0, 17.5/rad, 0.0, 10.5/rad, 12.50/rad, s, theta, mu)
				+  ellipProj(0.0, -5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+2*ellipProj(0.0, 5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+  ellipProj(0.0, -30.0/rad, 0.0, 1.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(-6.0/rad, -30.0/rad, 0.0, 2.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(5.0/rad, -30.0/rad, 0.0, 1.5/rad, 2.5/rad, s, theta, mu)
				-2*ellipProj(11.0/rad, 0.0, -18.0*pi/180.0, 5.5/rad, 15.5/rad, s, theta, mu)
				-2*ellipProj(-11.0/rad, 0.0, 18.0*pi/180.0, 8.0/rad, 20.5/rad, s, theta, mu);	
		}
	}
}

void ProjectionSimulator::ProjectorAPar(double **Radon, int angle, int bin, double mu)
{
	int i, k;
	double s, theta;
  const double pi = 3.141592654;
	double AngleStep = 2*pi/angle;

  for(k=0; k<angle; k++)
	{
		theta = k*AngleStep;
		for(i=0; i<bin; i++)
		{				
			s = (i+i-bin+1.0)/bin;
			Radon[k][i] = 
				 9*ellipProj(0.0, 0.0, 0.0, 35.0/rad, 46.0/rad, s, theta, mu)
				-4*ellipProj(0.0, -0.5/rad, 0.0, 33.0/rad, 43.50/rad, s, theta, mu)
				+  ellipProj(0.0, 17.5/rad, 0.0, 10.5/rad, 12.50/rad, s, theta, mu)
				+  ellipProj(0.0, -5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+2*ellipProj(0.0, 5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+  ellipProj(0.0, -30.0/rad, 0.0, 1.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(-6.0/rad, -30.0/rad, 0.0, 2.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(5.0/rad, -30.0/rad, 0.0, 1.5/rad, 2.5/rad, s, theta, mu)
				-2*ellipProj(11.0/rad, 0.0, -18.0*pi/180.0, 5.5/rad, 15.5/rad, s, theta, mu)
				-2*ellipProj(-11.0/rad, 0.0, 18.0*pi/180.0, 8.0/rad, 20.5/rad, s, theta, mu);
		}
	}
}

void ProjectionSimulator::ProjectorAFanA(double **Radon, int ANGLE, int RAY, double mu)
{
	int i, k;
	double s, theta;
  const double pi = 3.141592654;
	double projection_step = RAY_VIEW/RAY;
	double rotat_step = (pi+pi)/ANGLE;
 
	for(k=0; k<ANGLE; k++)
	{
		for(i=0; i<RAY; i++)
		{				
			s = RADIUS_FOCA*sin((i-RAY/2.+0.5)*projection_step);
			theta = k*rotat_step+(i-RAY/2.+0.5)*projection_step;

			Radon[k][i]
				=9*ellipProj(0.0, 0.0, 0.0, 35.0/rad, 46.0/rad, s, theta, mu)
				-4*ellipProj(0.0, -0.5/rad, 0.0, 33.0/rad, 43.50/rad, s, theta, mu)
				+  ellipProj(0.0, 17.5/rad, 0.0, 10.5/rad, 12.50/rad, s, theta, mu)
				+  ellipProj(0.0, -5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+2*ellipProj(0.0, 5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+  ellipProj(0.0, -30.0/rad, 0.0, 1.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(-6.0/rad, -30.0/rad, 0.0, 2.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(5.0/rad, -30.0/rad, 0.0, 1.5/rad, 2.5/rad, s, theta, mu)
				-2*ellipProj(11.0/rad, 0.0, -18.0*pi/180.0, 5.5/rad, 15.5/rad, s, theta, mu)
				-2*ellipProj(-11.0/rad, 0.0, 18.0*pi/180.0, 8.0/rad, 20.5/rad, s, theta, mu);
		}			
	}
}

void ProjectionSimulator::ProjectorAFanS(double **Radon, int ANGLE, int RAY, double mu)
{
	int i, k;
	double s, theta;
  const double pi = 3.141592654;
	double projection_step = RAY_FOV/RAY;
	double rotat_step = (pi+pi)/ANGLE;

	double ds = (RADIUS_FOCA+RADIUS_DETC)*(RADIUS_FOCA+RADIUS_DETC);
	double step, tilt;

	for(k=0;k<ANGLE;k++)
	{
		for(i=0;i<RAY;i++)
		{				
			step = (i-RAY/2.+0.5)*projection_step;
			tilt = step/sqrt(step*step+ds);

			s = RADIUS_FOCA*tilt;
			theta = (k*rotat_step+asin(tilt));

			Radon[k][i]
				=9*ellipProj(0.0, 0.0, 0.0, 35.0/rad, 46.0/rad, s, theta, mu)
				-4*ellipProj(0.0, -0.5/rad, 0.0, 33.0/rad, 43.50/rad, s, theta, mu)
				+  ellipProj(0.0, 17.5/rad, 0.0, 10.5/rad, 12.50/rad, s, theta, mu)
				+  ellipProj(0.0, -5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+2*ellipProj(0.0, 5.0/rad, 0.0, 2.5/rad, 2.5/rad, s, theta, mu)
				+  ellipProj(0.0, -30.0/rad, 0.0, 1.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(-6.0/rad, -30.0/rad, 0.0, 2.5/rad, 1.5/rad, s, theta, mu)
				+  ellipProj(5.0/rad, -30.0/rad, 0.0, 1.5/rad, 2.5/rad, s, theta, mu)
				-2*ellipProj(11.0/rad, 0.0, -18.0*pi/180.0, 5.5/rad, 15.5/rad, s, theta, mu)
				-2*ellipProj(-11.0/rad, 0.0, 18.0*pi/180.0, 8.0/rad, 20.5/rad, s, theta, mu);	
		}			
	}
}

void ProjectionSimulator::ProjectorDPar(double **SLImage, int dim, double **Radon, int angle)
{
	int i, k, j, a, b;
	double x, y, c1, c2, ix, iy;
	double HALFGRID = dim/2.0-0.5;
	double ANGLE_STEP = pi/angle;
	
	double *cos_ = new double[angle];
	double *sin_ = new double[angle];
	for(k=0; k<angle; k++)
	{
		cos_[k] = cos(k*ANGLE_STEP);
		sin_[k] = sin(k*ANGLE_STEP);
	}

	for(k=0; k<angle; k++)
  {
		for(j=0, iy=-HALFGRID; j<dim; j++, iy+=1.0)
		{
			Radon[k][j] = 0.0;
		  for(i=0, ix=-HALFGRID; i<dim; i++, ix+=1.0)
			{
				x = ix*cos_[k] - iy*sin_[k];
				y = ix*sin_[k] + iy*cos_[k];
				if(x<=HALFGRID && x>=-HALFGRID && y<=HALFGRID && y>=-HALFGRID)
				{
          y += HALFGRID; a = int(y); c1 = y - a;
          x += HALFGRID; b = int(x); c2 = x - b;
          if(a < dim-1 && b < dim-1)
					  Radon[k][j] += SLImage[a][b] + c1*(SLImage[a+1][b]-SLImage[a][b])	+
                           c2*(SLImage[a][b+1]-SLImage[a][b] +
                           c1*(SLImage[a+1][b+1]-SLImage[a][b+1]-SLImage[a+1][b]+SLImage[a][b]));
          else
            Radon[k][j] += SLImage[a][b];
				}
			}
		}
  }
	delete[] cos_; delete[] sin_;
}

void ProjectionSimulator::ProjectorDPar(double **SLImage, int dim, double **Radon, int angle, double atten)
{
	int i, k, j;
	double ANGLE_STEP = 2*pi/angle;
  double **IMAGE1 = CreateMatrix<double>(dim, dim);
  double step = RADIUS_IMAGE*2./dim;
  atten *= step;

	for(k=0; k<angle; k++)
  {
    RotateCarte(SLImage, IMAGE1, dim, k*ANGLE_STEP);
    for(i=0; i<dim; i++) 
    {
      Radon[k][i] = 0.0;
      for(j=0; j<dim; j++) 
        Radon[k][i] += IMAGE1[i][j]*exp(-atten*(j-dim/2.+0.5));
    }
  }

	FreeMatrix(IMAGE1);
}

void ProjectionSimulator::ProjectorDPar(double **SLImage, int dim, double **Radon, int angle, double **ATTENUATOR)
{
	int i, k, j;
	double ANGLE_STEP = pi/angle*2;
	double **ATTENU = CreateMatrix<double>(dim, dim);
  double **IMAGE1 = CreateMatrix<double>(dim, dim);

	for(k=0; k<angle; k++)
  {
		RotateCarte(ATTENUATOR, ATTENU, dim, k*ANGLE_STEP);
    RotateCarte(SLImage, IMAGE1, dim, k*ANGLE_STEP);
		for(i=0; i<dim; i++) for(j=dim-2; j>=0; j--) ATTENU[i][j] += ATTENU[i][j+1];

    for(i=0; i<dim; i++)
    {
      Radon[k][i] = 0.;
      for(j=0; j<dim; j++) Radon[k][i] += IMAGE1[i][j]*exp(-ATTENU[i][j]);
    }
  }

	FreeMatrix(ATTENU);
	FreeMatrix(IMAGE1);
}

void ProjectionSimulator::ProjectorDFanA(double **Radon, int angle, int ray, double **SLImage, int dim, double **atten)
{
  int k, m, n;
  double focalPos;
  double angle_step = 3.14159265358979324*2/angle;
  double **pImage = CreateMatrix<double>(ray, dim);
  double **attImage = 0;
  if(atten) attImage = CreateMatrix<double>(ray, dim);

  for(k=0; k<angle; k++)
  {
    focalPos = k*angle_step;
    CarteToPolarA(SLImage, dim, pImage, ray, focalPos);

    if( atten ) 
    {
      CarteToPolarA(atten, dim, attImage, ray, focalPos);
      for(m=0; m<ray; m++)
      {
        for(n=dim-2; n>=0; n--) attImage[m][n] += attImage[m][n+1];
      }
      for(m=0; m<ray; m++)
      {
        Radon[k][m] = 0.;
        for(n=0; n<dim; n++) Radon[k][m] += pImage[m][n]*exp(-attImage[m][n]);
      }
    }
    else
    {
      for(m=0; m<ray; m++)
      {
        Radon[k][m] = 0.;
        for(n=0; n<dim; n++) Radon[k][m] += pImage[m][n];
      }
    }
  }

  FreeMatrix(pImage);
  if(atten) FreeMatrix(attImage);
}

void ProjectionSimulator::ProjectorDFanS(double **Radon, int angle, int ray, double **SLImage, int dim, double **atten)
{
  int k, m, n;
  double focalPos;
  double angle_step = 3.14159265358979324*2/angle;
  double **pImage = CreateMatrix<double>(ray, dim);
  double **attImage = 0;
  if(atten) attImage = CreateMatrix<double>(ray, dim);

  for(k=0; k<angle; k++)
  {
    focalPos = k*angle_step;
    CarteToPolarS(SLImage, dim, pImage, ray, focalPos);

    if( atten ) 
    {
      CarteToPolarS(atten, dim, attImage, ray, focalPos);
      for(m=0; m<ray; m++)
      {
        for(n=dim-2; n>=0; n--) attImage[m][n] += attImage[m][n+1];
      }
      for(m=0; m<ray; m++)
      {
        Radon[k][m] = 0.;
        for(n=0; n<dim; n++) Radon[k][m] += pImage[m][n]*exp(-attImage[m][n]);
      }
      continue;
    }

    for(m=0; m<ray; m++)
    {
      Radon[k][m] = 0.;
      for(n=0; n<dim; n++) Radon[k][m] += pImage[m][n];
    }
  }

  FreeMatrix(pImage);
  if(atten) FreeMatrix(attImage);
}

void ProjectionSimulator::Smooth3(double *data, int len, double *temp)
{
  memcpy(temp, data, sizeof(double)*len); //0.035, -0.128,  0.070, 0.315, 0.417
  for(int i=1; i<len-1; i++)
  {
    data[i] = (temp[i-1]+temp[i+1])*0.25 + temp[i]*0.5;
  }
}

void ProjectionSimulator::Smooth5(double **IMAGE, int row, int col)
{
	int i, j;
	double **IMAGE1 = CreateMatrix<double>(row, col);

  for(i=1; i<row-1; i++)
  {
		for(j=1; j<col-1; j++)
			IMAGE1[i][j] = (4.*IMAGE[i][j]+IMAGE[i-1][j]+IMAGE[i][j-1]+IMAGE[i+1][j]+IMAGE[i][j+1])/8.;
  }
	for(i=1; i<row-1; i++)
  {
		for(j=1; j<col-1; j++) IMAGE[i][j] = IMAGE1[i][j];
  }

  FreeMatrix(IMAGE1);
}

void ProjectionSimulator::Smooth9(double **IMAGE, int row, int col)
{
	int i, j;
	double **IMAGE1 = CreateMatrix<double>(row, col);

  for(i=1; i<row-1; i++)
  {
		for(j=1; j<col-1; j++)
			IMAGE1[i][j] = (4.*IMAGE[i][j]+2.*(IMAGE[i-1][j]+IMAGE[i][j-1]+IMAGE[i+1][j]+IMAGE[i][j+1])+
                      IMAGE[i-1][j-1]+IMAGE[i-1][j+1]+IMAGE[i+1][j-1]+IMAGE[i+1][j+1])/16.;
  }
	for(i=1; i<row-1; i++)
  {
		for(j=1; j<col-1; j++) IMAGE[i][j] = IMAGE1[i][j];
  }

  FreeMatrix(IMAGE1);
}

void ProjectionSimulator::Smooth53(double *da1, int dim, double *da)
{
  if(dim < 5) return;
  memcpy(da, da1, sizeof(double)*dim);

  da1[0] = (69*da[0] + 4*(da[1]+da[3]) - 6*da[2] - da[4])/70;
  da1[1] = (2*(da[0]+da[4]) + 27*da[1] + 12*da[2] - 8*da[3])/35;
  da1[dim-2] = (2*(da[dim-1]+da[dim-5]) + 27*da[dim-2] + 12*da[dim-3] - 8*da[dim-4])/35;
  da1[dim-1] = (69*da[dim-1] + 4*(da[dim-2]+da[dim-4]) - 6*da[dim-3] - da[dim-5])/70;

  for(int i=2; i<dim-2; i++)
    da1[i] = (17*da[i] + 12*(da[i-1]+da[i+1]) - 3*(da[i-2]+da[i+2]))/35;
}

void ProjectionSimulator::SavGol5(double *data, int len, double *temp)
{
  memcpy(temp, data, sizeof(double)*len); //-0.086, 0.343,  0.486
  for(int i=2; i<len-2; i++)
  {
    data[i] = -(temp[i-2]+temp[i+2])*0.086 + (temp[i-1]+temp[i+1])*0.343 + temp[i]*0.486;
  }
}

void ProjectionSimulator::SavGol9(double *data, int len, double *temp)
{
  memcpy(temp, data, sizeof(double)*len); //0.035, -0.128,  0.070, 0.315, 0.417
  for(int i=4; i<len-4; i++)
  {
    data[i] = (temp[i-4]+temp[i+4])*0.035 - (temp[i-3]+temp[i+3])*0.128 +
              (temp[i-2]+temp[i+2])*0.070 + (temp[i-1]+temp[i+1])*0.315 + temp[i]*0.417;
  }
}

void ProjectionSimulator::FilterMedian(double *arr, int len)
{
  int i;
  double *tm = new double[len];
  memcpy(tm, arr, sizeof(double)*len);

  for(i=1; i<len-1; i++)
  {
    if( tm[i-1] < tm[i+1] ) 
    { 
      if(tm[i+1] < tm[i] ) arr[i] = tm[i+1];
      else if(tm[i] < tm[i-1] ) arr[i] = tm[i-1]; 
      else arr[i] = tm[i]; 
    }
    else 
    { 
      if(tm[i-1] < tm[i] ) arr[i] = tm[i-1]; 
      else if(tm[i] < tm[i+1] ) arr[i] = tm[i+1]; 
      else  arr[i] = tm[i]; 
    }
  }

  delete[] tm;
}

void ProjectionSimulator::FilterMedian(double **IMAGE, int row, int col)
{
  int i, j;
  double tb[3];
  double * pt;
  double **tm = CreateMatrix<double>(row, col);
  memcpy(tm[0], IMAGE[0], row*col*sizeof(double));

  for(i=1; i<row-1; i++)
  {
    pt = IMAGE[i];
    for(j=1; j<col-1; j++)
    {
      if( tm[i-1][j] <= tm[i+1][j] )
      {
        if( tm[i][j] < tm[i-1][j] )
        {
          tb[0] = tm[i][j];
          tb[1] = tm[i-1][j];
          tb[2] = tm[i+1][j];
        }
        else if( tm[i][j] > tm[i+1][j] )
        {
          tb[0] = tm[i-1][j];
          tb[1] = tm[i+1][j];
          tb[2] = tm[i][j];
        }
        else
        {
          tb[0] = tm[i-1][j];
          tb[1] = tm[i][j];
          tb[2] = tm[i+1][j];
        }
      }
      else
      {
        if( tm[i][j] > tm[i-1][j] )
        {
          tb[2] = tm[i][j];
          tb[1] = tm[i-1][j];
          tb[0] = tm[i+1][j];
        }
        else if( tm[i][j] < tm[i+1][j] )
        {
          tb[0] = tm[i][j];
          tb[1] = tm[i+1][j];
          tb[2] = tm[i-1][j];
        }
        else
        {
          tb[0] = tm[i+1][j];
          tb[1] = tm[i][j];
          tb[2] = tm[i-1][j];
        }
      }
      if( tm[i][j-1] <= tm[i][j+1] )
      {
        if( tm[i][j-1] >= tb[2] )
        {
          pt[j] = tb[2];
        }
        else if( tm[i][j+1] <= tb[0] )
        {
          pt[j] = tb[0];
        }
        else if( tm[i][j-1] >= tb[1] )
        {
          pt[j] = tm[i][j-1];
        }
        else if( tm[i][j+1] <= tb[1] )
        {
          pt[j] = tm[i][j+1];
        }
        else
        {
          pt[j] = tb[1];
        }
      }
      else
      {
        if( tm[i][j+1] >= tb[2] )
        {
          pt[j] = tb[2];
        }
        else if( tm[i][j-1] <= tb[0] )
        {
          pt[j] = tb[0];
        }
        else if( tm[i][j+1] >= tb[1] )
        {
          pt[j] = tm[i][j+1];
        }
        else if( tm[i][j-1] <= tb[1] )
        {
          pt[j] = tm[i][j-1];
        }
        else
        {
          pt[j] = tb[1];
        }
      }
    }
  }

  FreeMatrix(tm);
}

void ProjectionSimulator::RotateCarte(double **IMAGE, double **OUTPUT, int SIZE, double ANGLE)
{
	int i, j, a, b;
	double x, y, ix, iy, c1, c2;
	double HALFGRID = SIZE/2.-0.5;
	double cos_ = cos(ANGLE);
  double sin_ = sin(ANGLE);

	for(i=0, iy=-HALFGRID; i<SIZE; i++, iy+=1.)
  {
		for(j=0, ix=-HALFGRID; j<SIZE; j++, ix+=1.)
		{
			x = ix*cos_ - iy*sin_;
			y = ix*sin_ + iy*cos_;
			if(x<=HALFGRID && x>=-HALFGRID && y<=HALFGRID && y>=-HALFGRID)
			{	
        y += HALFGRID; a = int(y); c1 = y - a; 
        x += HALFGRID; b = (int)(x); c2 = x - b; 
        if(a < SIZE-1 && b < SIZE-1)
        {
				  OUTPUT[i][j] = (IMAGE[a][b]+c1*(IMAGE[a+1][b]-IMAGE[a][b]) + 
	  				              c2*(IMAGE[a][b+1]-IMAGE[a][b] +
					                c1*(IMAGE[a+1][b+1]-IMAGE[a][b+1]-IMAGE[a+1][b]+IMAGE[a][b])));
        }
        else
          OUTPUT[i][j] = IMAGE[a][b];
			}
			else
				OUTPUT[i][j] = 0.;
		}
  }
}

void ProjectionSimulator::RotateCarte(double **IMAGE, float **OUTPUT, int SIZE, double ANGLE)
{
	int i, j, a, b;
	double x, y, ix, iy, c1, c2;
	double HALFGRID = SIZE/2.-0.5;
	double cos_ = cos(ANGLE);
  double sin_ = sin(ANGLE);

	for(i=0, iy=-HALFGRID; i<SIZE; i++, iy+=1.)
  {
		for(j=0, ix=-HALFGRID; j<SIZE; j++, ix+=1.)
		{
			x = ix*cos_ - iy*sin_;
			y = ix*sin_ + iy*cos_;
			if(x<=HALFGRID && x>=-HALFGRID && y<=HALFGRID && y>=-HALFGRID)
			{	
        y += HALFGRID; a = int(y); c1 = y - a; 
        x += HALFGRID; b = (int)(x); c2 = x - b; 
        if(a < SIZE-1 && b < SIZE-1)
        {
				  OUTPUT[i][j] = float(IMAGE[a][b]+c1*(IMAGE[a+1][b]-IMAGE[a][b]) + 
	  				              c2*(IMAGE[a][b+1]-IMAGE[a][b] +
					                c1*(IMAGE[a+1][b+1]-IMAGE[a][b+1]-IMAGE[a+1][b]+IMAGE[a][b])));
        }
        else
          OUTPUT[i][j] = float(IMAGE[a][b]);
			}
			else
				OUTPUT[i][j] = 0.;
		}
  }
}

void ProjectionSimulator::CarteToPolarA(double **cImage, int SIZE, double ** pImage, int ANGLE, double focalPos)
{
	int i, j, m, k;
  double HALFGRID = SIZE/2.-0.5;
	double ray_step = RAY_VIEW/ANGLE;
	double FOCA_SHIFT = SIZE*(RADIUS_FOCA-RADIUS_DETC)/(2*RADIUS_DETC);
  double focalX = SIZE*RADIUS_FOCA*cos(focalPos)/(2*RADIUS_DETC);
  double focalY = SIZE*RADIUS_FOCA*sin(focalPos)/(2*RADIUS_DETC);
  double c, d, cx, cy;

	double *cos_ = new double[ANGLE];
	double *sin_ = new double[ANGLE];
	for(k=0; k<ANGLE; k++)
	{
		cos_[k] = cos((k-ANGLE/2.+0.5)*ray_step + focalPos);
		sin_[k] = sin((k-ANGLE/2.+0.5)*ray_step + focalPos);
	}

	for(k=0; k<ANGLE; k++)
  {
		for(m=0; m<SIZE; m++)
		{
      cy = (FOCA_SHIFT+m)*sin_[k] - focalY;
      cx = (FOCA_SHIFT+m)*cos_[k] - focalX;
			if(cx<=HALFGRID && cx>=-HALFGRID && cy<=HALFGRID && cy>=-HALFGRID)
			{
        cy += HALFGRID; i = int(cy); c = cy - i;
        cx += HALFGRID; j = int(cx); d = cx - j;
        if(i < SIZE-1 && j < SIZE-1 )
				  pImage[k][m] = cImage[i][j]+c*(cImage[i+1][j]-cImage[i][j]) +
	  					           d*(cImage[i][j+1]-cImage[i][j] + 
						             c*(cImage[i+1][j+1]-cImage[i][j+1]-cImage[i+1][j]+cImage[i][j]));
        else
          pImage[k][m] = cImage[i][j];
			}
			else
				pImage[k][m] = 0.0;
		}
  }

  delete[] cos_; delete[] sin_;
}

void ProjectionSimulator::PolarAToCarte(double ** pImage, int ANGLE, double focalPos, double **cImage, int SIZE)
{
	int i, j, m, k;
  double HALFGRID = SIZE/2.-0.5;
  double HALFVIEW = RAY_VIEW/2.;
	double ray_step = RAY_VIEW/ANGLE;
	double FOCA_SHIFT = SIZE*(RADIUS_FOCA-RADIUS_DETC)/(2*RADIUS_DETC);
  double focalX = SIZE*RADIUS_FOCA*cos(focalPos+pi)/(2*RADIUS_DETC);
  double focalY = SIZE*RADIUS_FOCA*sin(focalPos+pi)/(2*RADIUS_DETC);
  double radius, phi, px, py, c, d;

	for(i=0; i<SIZE; i++)
  {
		for(j=0; j<SIZE; j++)
		{
      py = i-HALFGRID - focalY;
      px = j-HALFGRID - focalX;
			radius = sqrt(py*py + px*px);
      if(radius < FOCA_SHIFT)
      {
        cImage[i][j] = 0.0;
        continue;
      }
			phi = atan2(py, px); //if(phi < 0) phi += (pi+pi);
      phi -= focalPos;
      if(phi < -pi) phi += (pi+pi);

			if(phi <= HALFVIEW && phi >= -HALFVIEW)
			{
        phi = (phi+RAY_VIEW/2.)/ray_step; k = int(phi); c = phi-k;
        radius -= FOCA_SHIFT; m = int(radius); d = radius-m;
        if(k < ANGLE-1 && m < SIZE-1)
				  cImage[i][j] = pImage[k][m]+c*(pImage[k+1][m]-pImage[k][m]) +
	  					           d*(pImage[k][m+1]-pImage[k][m] +
						             c*(pImage[k+1][m+1]-pImage[k][m+1]-pImage[k+1][m]+pImage[k][m]));
        else
          cImage[i][j] = pImage[k][m];
			}
      else
        cImage[i][j] = 0.0;
		}
  }
}

void ProjectionSimulator::CarteToPolarS(double **cImage, int SIZE, double ** pImage, int ANGLE, double focalPos)
{
	int i, j, m, k;
  double HALFGRID = SIZE/2.-0.5;
	double ray_step = RAY_FOV/ANGLE;
	double FOCA_SHIFT = SIZE*(RADIUS_FOCA-RADIUS_DETC)/(2*RADIUS_DETC);
  double focalX = SIZE*RADIUS_FOCA*cos(focalPos)/(2*RADIUS_DETC);
  double focalY = SIZE*RADIUS_FOCA*sin(focalPos)/(2*RADIUS_DETC);
  double c, d, cx, cy, step, tilt, theta;
	double ds = (RADIUS_FOCA+RADIUS_DETC)*(RADIUS_FOCA+RADIUS_DETC);

	double *cos_ = new double[ANGLE];
	double *sin_ = new double[ANGLE];
	for(k=0; k<ANGLE; k++)
	{
		step = (k-ANGLE/2.+0.5)*ray_step;
		tilt = step/sqrt(step*step+ds);
		theta = asin(tilt);

		cos_[k] = cos(theta + focalPos);
		sin_[k] = sin(theta + focalPos);
	}

	for(k=0; k<ANGLE; k++)
  {
		for(m=0; m<SIZE; m++)
		{
      cy = (FOCA_SHIFT+m)*sin_[k] - focalY;
      cx = (FOCA_SHIFT+m)*cos_[k] - focalX;
			if(cx<=HALFGRID && cx>=-HALFGRID && cy<=HALFGRID && cy>=-HALFGRID)
			{
        cy += HALFGRID; i = int(cy); c = cy - i;
        cx += HALFGRID; j = int(cx); d = cx - j;
        if(i < SIZE-1 && j < SIZE-1 )
				  pImage[k][m] = cImage[i][j]+c*(cImage[i+1][j]-cImage[i][j]) +
	  					           d*(cImage[i][j+1]-cImage[i][j] + 
						             c*(cImage[i+1][j+1]-cImage[i][j+1]-cImage[i+1][j]+cImage[i][j]));
        else
          pImage[k][m] = cImage[i][j];
			}
			else
				pImage[k][m] = 0.0;
		}
  }

  delete[] cos_; delete[] sin_;
}

void ProjectionSimulator::PolarSToCarte(double **polImage, int angle, double focalPos, double **carImage, int dim)
{
}

void ProjectionSimulator::CarteToConeP(double ***carImage, int dimZ, int dimXY, double ***coneImage, int theta, int phi, double focalPos)
{
  // to be implemented
}

void ProjectionSimulator::ConePToCarte(double ***coneImage, int theta, int phi, double focalPos, double ***carImage, int dimZ, int dimXY)
{
  // to be implemented
}
