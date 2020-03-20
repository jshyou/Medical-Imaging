// ReconIterative.cpp: implementation of the ReconIterative class.
//
//////////////////////////////////////////////////////////////////////

#include "ReconIterative.h"
#include "ProjectionSimulator.h"
#include <math.h>
#include <memory.h>

const double pi  = 3.1415926535897932385;

ReconIterative::ReconIterative() 
 : angle_step(0.), angleNum(0), AA(0), BB(0) 
{ 
}

ReconIterative::~ReconIterative() 
{ 
  delete[] AA; 
  delete[] BB; 
}

void ReconIterative::Init(int ANGLE_NUM)
{
  delete[] AA; delete[] BB;
  AA = new double[ANGLE_NUM];
  BB = new double[ANGLE_NUM];

  int k;
	for(k=0; k<ANGLE_NUM; k++) BB[k] = cos(k*angle_step);
	for(k=0; k<ANGLE_NUM; k++) AA[k] = sin(-k*angle_step);

  angleNum = ANGLE_NUM;
}

/////////////////////////////////////////////////////////
// EM algorithm for parallel-beam acquisition geometry //
/////////////////////////////////////////////////////////

void ReconIterative::EM_Parall(double **RADON, int ANGLE, int BIN, double **IMAGE, int step, int iter)
{
  angle_step = pi/ANGLE;
  Init(ANGLE);

  int i, j, k, m, n, iteration;
	double **IMAGE1 = CreateMatrix<double>(BIN, BIN);
	double **Radon1 = CreateMatrix<double>(ANGLE, BIN);

  // total sum of projections
  double prjSum = 0.;
	for(k=0; k<ANGLE; k++) for(i=0; i<BIN; i++) prjSum += RADON[k][i];
	double eps = 1.e-3*prjSum/(ANGLE*BIN);

  // detect the contour of the image to be reconstructed
	for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) IMAGE[m][n] = 1.;
	for(k=0; k<ANGLE; k++)
	{
		backprj(RADON[k], k, BIN, IMAGE1);
		for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) if(IMAGE1[m][n] < eps) IMAGE[m][n] = 0.;
	}

  // perform iteration for subgrouped projections
  iteration = 0;
	while(iteration < iter)
	{
		for(j=0; j<step; j++)
		{
			project(IMAGE, BIN, Radon1, ANGLE, j, step);
      k = j;
      while(k < ANGLE)
      {
			  for(i=0; i<BIN; i++)
			  {
				  if(Radon1[k][i] > eps && RADON[k][i] > eps) Radon1[k][i] = RADON[k][i]/Radon1[k][i];
				  else Radon1[k][i] = 0.0;
			  }
        k += step;
      }
			backprj(Radon1, ANGLE, j, step, IMAGE1, BIN);
			for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) if( IMAGE[m][n] > 0.0 ) IMAGE[m][n] *= IMAGE1[m][n];
		}
		iteration++;
	}

  FreeMatrix(Radon1);
	FreeMatrix(IMAGE1);
}

void ReconIterative::EM_Parall(double **RADON, int ANGLE, int BIN, double **IMAGE, double **ATTENMap, int step, int iter, bool half)
{
  if( !half ) angle_step = 2*pi/ANGLE;
  else angle_step = pi/ANGLE;
  Init(ANGLE);

	int i, j, k, m, n, iteration;

	double **IMAGE1 = CreateMatrix<double>(BIN, BIN);
  double **Radon1 = CreateMatrix<double>(ANGLE, BIN);
  double **ATTENU = CreateMatrix<double>(BIN, BIN);

  double prjSum = 0.;
	for(k=0; k<ANGLE; k++) for(i=0; i<BIN; i++) prjSum += RADON[k][i];
	double eps = 1.e-3*prjSum/(ANGLE*BIN);

  // detect image contour
	for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) IMAGE[m][n] = 1.;
	for(k=0; k<ANGLE; k++)
	{
    backprj(RADON[k], k, BIN, IMAGE1);
		for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) if(IMAGE1[m][n] < eps) IMAGE[m][n] = 0.;
	}

  // perform iteration for subgrouped projections
  iteration = 0;
	while(iteration < iter)
	{
		for(j=0; j<step; j++)
		{
			project(IMAGE, ATTENMap, BIN, Radon1, ANGLE, j, step, ATTENU, IMAGE1);
      k = j;
      while(k < ANGLE)
      {
			  for(i=0; i<BIN; i++)
			  {
				  if(Radon1[k][i] > eps && RADON[k][i] > eps) Radon1[k][i] = RADON[k][i]/Radon1[k][i];
				  else Radon1[k][i] = 0.0;
			  }
        k += step;
      }
			backprj(Radon1, ANGLE, j, step, IMAGE1, BIN);
			for(m=0; m<BIN; m++) for(n=0; n<BIN; n++) if( IMAGE[m][n] > 0.0 ) IMAGE[m][n] *= IMAGE1[m][n];
		}
		iteration++;
	}

  FreeMatrix(Radon1);
	FreeMatrix(IMAGE1);
  FreeMatrix(ATTENU);
}

void ReconIterative::project(double **IMAGE, int SIZE, double **Radon, int views, int start, int step)
{
	int i, k, j, a, b;
	double x, y, c1, c2, ix, iy, HALFGRID = SIZE/2.-0.5;
	
  k = start;
	while(k < views)
  {
		for(i=0, ix=-HALFGRID; i<SIZE; i++, ix+=1.0)
		{
			Radon[k][i] = 0.0;
			for(j=0, iy=-HALFGRID; j<SIZE; j++, iy+=1.0)
			{
				x = ix*BB[k]-iy*AA[k]; y = ix*AA[k]+iy*BB[k];
				if(x<=HALFGRID && x>=-HALFGRID && y<=HALFGRID && y>=-HALFGRID)
				{
          x += HALFGRID; y += HALFGRID;
					a = int(x); c1 = x - a; if(a==SIZE-1) a--;
					b = int(y); c2 = y - b; if(b==SIZE-1) b--;
					
					Radon[k][i] += IMAGE[a][b]+c1*(IMAGE[a+1][b]-IMAGE[a][b]) +
		  				           c2*(IMAGE[a][b+1]-IMAGE[a][b] +
						             c1*(IMAGE[a+1][b+1]-IMAGE[a][b+1]-IMAGE[a+1][b]+IMAGE[a][b]));					
				}
			}
		}
    k += step;
  }
}

void ReconIterative::backprj(double **Radon, int views, int start, int step, double **IMAGE, int BIN)
{
	int k, m, n, j;
	double a, y, mm, nn, HALFGRID = BIN/2.-0.5;

	for(m=0, mm=-HALFGRID; m<BIN; m++, mm+=1.0)
  {
		for(n=0, nn=-HALFGRID; n<BIN; n++, nn+=1.0)
		{
			IMAGE[m][n] = 0.0;
      k = start;
			while(k < views)
			{
				y = BB[k]*mm + AA[k]*nn;
				if(y <= HALFGRID && y >= -HALFGRID)
				{
					y += HALFGRID; j = int(y); a = y-j;
					if(j < BIN-1) IMAGE[m][n] += Radon[k][j]+a*(Radon[k][j+1]-Radon[k][j]);
					else IMAGE[m][n] += Radon[k][j];
				}
        k += step;
			}
      IMAGE[m][n] /= (views/step);
		}
  }
}

void ReconIterative::backprj(double *Radon, int k, int BIN, double **IMAGE)
{
	int m, n, j;
	double a, x, iy, ix, HALFGRID = BIN/2.-0.5;
  double cos_ = cos(k*angle_step);
  double sin_ = sin(k*angle_step);

	for(m=0, iy=-HALFGRID; m<BIN; m++, iy+=1.0)
  {
		for(n=0, ix=-HALFGRID; n<BIN; n++, ix+=1.0)
		{
			IMAGE[m][n] = 0.0;
			x = cos_*iy - sin_*ix;
			if(x<=HALFGRID && x>=-HALFGRID)
			{
				x += HALFGRID; j = int(x); a = x-j;
				if(j < BIN-1)
					IMAGE[m][n] += Radon[j]+a*(Radon[j+1]-Radon[j]);
				else
					IMAGE[m][n] += Radon[j];
			}
		}
  }
}

void ReconIterative::project(double **IMAGE, double **ATTENMap, int SIZE, double **Radon, int views, int start, int step, double **ATTENU, double **IMAGE1)
{
	int i, j, k;

  k = start;
  while(k < views)
  {
    // rotate the image and attenuation map
    ProjectionSimulator::RotateCarte(IMAGE, IMAGE1, SIZE, k*angle_step);
    ProjectionSimulator::RotateCarte(ATTENMap, ATTENU, SIZE, k*angle_step);

    // line integral of attenuation map
	  for(i=0; i<SIZE; i++) for(j=SIZE-2; j>=0; j--) ATTENU[i][j] += ATTENU[i][j+1];

    // attenuated Radon transform
    for(i=0; i<SIZE; i++)
    {
      Radon[k][i] = 0.;
      for(j=0; j<SIZE; j++) Radon[k][i] += IMAGE1[i][j]*exp(-ATTENU[i][j]);
    }
    k += step;
  }
}

//////////////////////////////////////////////////////
// reconstruction for fan-beam acquisition geometry //
//////////////////////////////////////////////////////

void ReconIterative::EM_FanAng(double **RADON, int ANGLE, int RAY, double **IMAGE, int SIZE, int step, int OSEM_Iter)
{
  angle_step = 2*pi/ANGLE;
  double ray_step = RAY_VIEW/RAY;
	double eps =1.e-8;

	int i, j, k, m, n, iteration;

  double **fanPrj = CreateMatrix<double>(ANGLE, RAY);
	double **cImage = CreateMatrix<double>(SIZE, SIZE);
  double **sImage = CreateMatrix<double>(SIZE, SIZE);
  double **pImage = CreateMatrix<double>(RAY, SIZE);

  // detect image contour
	for(m=0; m<SIZE; m++)
  {
		for(n=0; n<SIZE; n++) IMAGE[m][n] = 1.;
  }
	for(k=0; k<ANGLE; k++)
  {
    backprj_fanA(RADON[k], RAY, cImage, SIZE, k*angle_step, pImage);
	  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) if(cImage[m][n] < eps) IMAGE[m][n] = 0.;
  }

  iteration = 0;
	while(iteration < OSEM_Iter)
  {
		for(j=0; j<step; j++)
		{
      k = j;
      while(k < ANGLE)
      {
			  project_fanA(fanPrj[k], RAY, IMAGE, SIZE, k*angle_step, pImage);
			  for(i=0; i<RAY; i++)
			  {
				  if(fanPrj[k][i] > eps && RADON[k][i] > eps) fanPrj[k][i] = RADON[k][i]/fanPrj[k][i];
				  else fanPrj[k][i] = 0.0;
			  }
        k += step;
      }

      k = j;
      for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) sImage[m][n] = 0.;
      while(k < ANGLE)
      {
			  backprj_fanA(fanPrj[k], RAY, cImage, SIZE, k*angle_step, pImage);
        for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) sImage[m][n] += cImage[m][n];
        k += step;
      }

			for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= (sImage[m][n]*step/ANGLE);
		}

    ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
		iteration++;
  }

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

	FreeMatrix(fanPrj);
	FreeMatrix(cImage);
  FreeMatrix(pImage);
  FreeMatrix(sImage);
}

void ReconIterative::EM_FanAng(double **RADON, int ANGLE, int RAY, double **IMAGE, int SIZE, double **atten, int step, int OSEM_Iter)
{
  angle_step = 2*pi/ANGLE;
  double ray_step = RAY_VIEW/RAY;
	double eps =1.e-8;

	int i, j, k, m, n, iteration;

  double **fanPrj = CreateMatrix<double>(ANGLE, RAY);
	double **cImage = CreateMatrix<double>(SIZE, SIZE);
  double **pImage = CreateMatrix<double>(RAY, SIZE);
  double **sImage = CreateMatrix<double>(SIZE, SIZE);
  double **pAtten = CreateMatrix<double>(RAY, SIZE);

  // detect image contour
	for(m=0; m<SIZE; m++)
  {
		for(n=0; n<SIZE; n++) IMAGE[m][n] = 1.;
  }
	for(k=0; k<ANGLE; k++)
  {
    backprj_fanA(RADON[k], RAY, cImage, SIZE, k*angle_step, pImage);
	  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) if(cImage[m][n] < eps) IMAGE[m][n] = 0.;
  }

  iteration = 0;
	while(iteration < OSEM_Iter)
  {
		for(j=0; j<step; j++)
		{
      k = j;
      while(k < ANGLE)
      {
        project_fanA(fanPrj[k], RAY, atten, IMAGE, SIZE, k*angle_step, pImage, pAtten);
			  for(i=0; i<RAY; i++)
			  {
				  if(fanPrj[k][i] > eps && RADON[k][i] > eps) fanPrj[k][i] = RADON[k][i]/fanPrj[k][i];
				  else fanPrj[k][i] = 0.0;
			  }
        k += step;
      }

      k = j;
      for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) sImage[m][n] = 0.;
      while(k < ANGLE)
      {
			  backprj_fanA(fanPrj[k], RAY, cImage, SIZE, k*angle_step, pImage);
        for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) sImage[m][n] += cImage[m][n];
        k += step;
      }

			for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= (sImage[m][n]*step/ANGLE);
		}

    ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
		iteration++;
  }

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

	FreeMatrix(fanPrj);
	FreeMatrix(cImage);
  FreeMatrix(pImage);
  FreeMatrix(sImage);
  FreeMatrix(pAtten);
}

void ReconIterative::project_fanA(double *Radon, int RAY, double **IMAGE, int SIZE, double pos, double **pImage)
{
  ProjectionSimulator::CarteToPolarA(IMAGE, SIZE, pImage, RAY, pos);
  int i, j;
  for(i=0; i<RAY; i++)
  {
    Radon[i] = 0.;
    for(j=0; j<SIZE; j++) Radon[i] += pImage[i][j];
  }
}

void ReconIterative::project_fanA(double *Radon, int RAY, double **atten, double **IMAGE, int SIZE, double pos, double **pImage, double **pAtten)
{
  ProjectionSimulator::CarteToPolarA(IMAGE, SIZE, pImage, RAY, pos);
  ProjectionSimulator::CarteToPolarA(atten, SIZE, pAtten, RAY, pos);
  int i, j;
  for(i=0; i<RAY; i++)
  {
    for(j=SIZE-2; j>=0; j--) pAtten[i][j] += pAtten[i][j+1];
    
    Radon[i] = 0.;
    for(j=0; j<SIZE; j++) Radon[i] += pImage[i][j]*exp(-pAtten[i][j]);
  }
}

void ReconIterative::backprj_fanA(double *Radon, int RAY, double **cImage, int SIZE, double pos, double **pImage)
{
  int i, j;
  for(i=0; i<RAY; i++)
  {
    for(j=0; j<SIZE; j++) pImage[i][j] = Radon[i];
  }
  ProjectionSimulator::PolarAToCarte(pImage, RAY, pos, cImage, SIZE);
}

