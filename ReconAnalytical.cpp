// ReconAnalytical.cpp: implementation of the ReconAnalytical class.
//
//////////////////////////////////////////////////////////////////////

#include "ReconAnalytical.h"
#include "ProjectionSimulator.h"
#include <math.h>
#include <memory.h>

const double pi  = 3.1415926535897932385;
const double ln2 = 0.69314718055994529;
const double eps = 1.e-8;

ReconAnalytical::ReconAnalytical()
{
}

ReconAnalytical::~ReconAnalytical()
{
}

void ReconAnalytical::FBP_ParSL(double **RADON,
                                int ANGLE,
                                int RAYS,
                                double **IMAGE,
                                int SIZE,
                                int filter)
{
  ///////////////////////////
  // filtering projections //
  ///////////////////////////
  int m, n, i, k, j;
	double *Filter = new double[RAYS*2+1];
  double cc = 0.5/(pi*ANGLE);
	for(i=0; i<=RAYS*2; i++) Filter[i] = cc/(0.25-(i-RAYS)*(i-RAYS)); //Shepp-Logan filter

  if(filter >= 2) ProjectionSimulator::FilterMedian(RADON, ANGLE, RAYS);

	double *R = new double[RAYS];
	for(k=0; k<ANGLE; k++)
	{
		memcpy(R, RADON[k], sizeof(double)*RAYS);
		for(m=0; m<RAYS; m++)
		{
			RADON[k][m] = 0.0;
			for(n=0; n<RAYS; n++) RADON[k][m] += R[n]*Filter[n-m+RAYS];
		}

    if(filter >= 1) ProjectionSimulator::SavGol5(RADON[k], RAYS, R);
	}
	delete[] R; delete[] Filter;

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

	double HALFGRID = SIZE/2.-0.5;
	double ANGLE_STEP = pi/ANGLE;
  double ratio = double(RAYS)/SIZE;
	double r, s, c;
  double half_m, half_m2, half_n;
  double halfRat = HALFGRID*ratio;
  double sqrR = HALFGRID*HALFGRID;
  double* cos_ = new double[ANGLE];
  double* sin_ = new double[ANGLE];
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP)*ratio;
    sin_[i] = sin(i*ANGLE_STEP)*ratio;
  }

	for(m=0; m<SIZE; m++)
  {
    half_m = m - HALFGRID;
    half_m2 = half_m*half_m;
		for(n=0; n<SIZE; n++)
		{
      half_n = n - HALFGRID; r = half_m2 + half_n*half_n; if(r >= sqrR) continue;
			
      IMAGE[m][n] = 0.0;
			for(k=0; k<ANGLE; k++)
			{
				s = half_m*cos_[k] - half_n*sin_[k] + halfRat;
				j = int(s); c = s - j;
				if(j < RAYS-1) IMAGE[m][n] += RADON[k][j] + c*(RADON[k][j+1]-RADON[k][j]);
        else IMAGE[m][n] += RADON[k][j];
			}
    }
  }

  delete[] cos_; delete[] sin_;
}

void ReconAnalytical::FBP_ParER(double **Radon,
                                int ANGLE,
                                int RAYS,
                                double **IMAGE,
                                int SIZE,
                                double mu,
                                int method,
                                int filter)
{
  ///////////////////////////
  // filtering projections //
  ///////////////////////////
  int m, n, i, j, k;
	double HALFGRID = SIZE/2.-0.5;
	double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double bin_step = IMAGE_DIAMETER/RAYS;
	double attenu = mu/(pi*2);

	double *R = new double[RAYS];
	double *Filter = new double[RAYS*2+1];

  if(filter >= 2) ProjectionSimulator::FilterMedian(Radon, ANGLE, RAYS);

  if(method == 1) // for Tretiak-Metz formula
  {
	  for(i=0; i<=2*RAYS; i++) Filter[i] = 0.5*RampSL((i-RAYS)*bin_step, bin_step, attenu);
	  for(k=0; k<ANGLE; k++)
	  {
		  memcpy(R, Radon[k], sizeof(double)*RAYS);
		  for(m=0; m<RAYS; m++)
		  {
			  Radon[k][m] = 0.0;
			  for(n=0; n<RAYS; n++) Radon[k][m] += R[n]*Filter[m-n+RAYS];
		  } 
      if(filter >= 1) ProjectionSimulator::SavGol5(Radon[k], RAYS, R);
	  }
  }
  else if(method == 2) // for Novikov formula
  {
    for(i=0; i<=2*RAYS; i++) Filter[i] = HilbertAttCos((i-RAYS)*bin_step, bin_step, mu);
	  for(k=0; k<ANGLE; k++)
	  {
		  memcpy(R, Radon[k], sizeof(double)*RAYS);
		  for(m=0; m<RAYS; m++)
		  {
			  Radon[k][m] = 0.0;
			  for(n=0; n<RAYS; n++) Radon[k][m] += R[n]*Filter[m-n+RAYS];
		  } 
      for(m=1; m<RAYS-1; m++) R[m] = (Radon[k][m+1] - Radon[k][m-1])*0.5;
      R[0] = R[1]; R[RAYS-1] = R[RAYS-2];
      memcpy(Radon[k], R, sizeof(double)*RAYS);
      if(filter >= 1) ProjectionSimulator::SavGol5(Radon[k], RAYS, R);
	  }
    /*
    double **prj  = CreateMatrix<double>(ANGLE, RAYS);
    double *aCos = new double[RAYS];
    double *aSin = new double[RAYS];
    for(m=0; m<RAYS; m++)
    {
      aCos[m] = cos((m-RAYS/2.+0.5)*bin_step*mu);
      aSin[m] = sin((m-RAYS/2.+0.5)*bin_step*mu);
    }
	  for(k=0; k<ANGLE; k++)
	  {
		  for(m=0; m<RAYS; m++)
      {
        prj[k][m] = Radon[k][m]*aSin[m];
        Radon[k][m] = Radon[k][m]*aCos[m];
      }
    }
	  for(i=0; i<=2*RAYS; i++) Filter[i] = 0.25*HilbertKL(i-RAYS)*bin_step;

	  for(k=0; k<ANGLE; k++)
	  {
		  memcpy(R, Radon[k], sizeof(double)*RAYS);
		  for(m=0; m<RAYS; m++)
		  {
			  Radon[k][m] = 0.0;
			  for(n=0; n<RAYS; n++) Radon[k][m] += R[n]*Filter[m-n+RAYS];
		  }
      for(m=0; m<RAYS; m++) Radon[k][m] *= aCos[m];
      for(m=1; m<RAYS-1; m++) R[m] = (Radon[k][m+1] - Radon[k][m-1])*0.5;
      memcpy(Radon[k], R, sizeof(double)*RAYS);

		  memcpy(R, prj[k], sizeof(double)*RAYS);
		  for(m=0; m<RAYS; m++)
		  {
			  prj[k][m] = 0.0;
			  for(n=0; n<RAYS; n++) prj[k][m] += R[n]*Filter[m-n+RAYS];
		  }
      for(m=0; m<RAYS; m++) prj[k][m] *= aSin[m];
      for(m=1; m<RAYS-1; m++) R[m] = (prj[k][m+1] - prj[k][m-1])*0.5;
      memcpy(prj[k], R, sizeof(double)*RAYS);

      for(m=0; m<RAYS; m++) Radon[k][m] += prj[k][m];
      if(filter >= 1) ProjectionSimulator::SavGol5(Radon[k], RAYS, R);
	  }

    delete[] aCos;
    delete[] aSin;
    FreeMatrix(prj);*/
  }

	delete []R; delete []Filter;

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

	double r, s, c, s1, r1, mm, nn;
	double recon_step = IMAGE_DIAMETER/SIZE;
  double sqr = HALFGRID*HALFGRID;
	double angle_step = 2*pi/ANGLE;
	r1 = mu*recon_step;
  s = recon_step/bin_step;
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*angle_step);
    sin_[i] = sin(i*angle_step);
  }

  for(m=0; m<SIZE; m++)
  {
    mm = m-HALFGRID;
		for(n=0; n<SIZE; n++)
		{
      nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;
			
      IMAGE[m][n] = 0.0;
			for(k=0; k<ANGLE; k++)
			{
				s1 = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*s; j = int(s1); c = s1-j;
				if(j < RAYS-1)
					IMAGE[m][n] += exp(r1*(cos_[k]*nn + sin_[k]*mm))*(Radon[k][j]+c*(Radon[k][j+1]-Radon[k][j]));
				else
					IMAGE[m][n] += Radon[k][j];
			}
		}
  }
  
  delete[] cos_; delete[] sin_;
}

// implementation of FBP with classical filters based on Novikov's formula
void ReconAnalytical::FBP_ParYJ(double **attRadon,  // attenuated projections
                                int ANGLE,          // viewing angle numbers
                                int RAYS,           // projection bins in each view
                                double **IMAGE,     // image to be reconstructed
                                int SIZE,           // the size of reconstructed image
                                double **Atten,     // the attenuation map, has to have the same size with projection bins
                                int filter
                                )
{
  ////////////////////////////////////
  // preprocessing attenuation map  //
  ////////////////////////////////////

  if(filter >= 2) ProjectionSimulator::FilterMedian(attRadon, ANGLE, RAYS);

	int    i, k, j, m, n, jj;
	double ANGLE_STEP = 2*pi/ANGLE;
  double *hilbert  = new double[2*RAYS+1];
  float  ***divAtt = CreateVolume<float>(ANGLE, RAYS, RAYS);
  double **attPrj  = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjH = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjC = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjS = CreateMatrix<double>(ANGLE, RAYS);

  for(i=0; i<2*RAYS+1; i++) hilbert[i] = HilbertKL(i-RAYS);

  for(k=0; k<ANGLE; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, divAtt[k], RAYS, k*ANGLE_STEP);
		for(i=0; i<RAYS; i++)
    {
			for(j=RAYS-2; j>=0; j--) divAtt[k][i][j] += divAtt[k][i][j+1];
      attPrj[k][i] = divAtt[k][i][0]/2.;
		  for(j=0; j<RAYS; j++) 
      {
        divAtt[k][i][j] -= float(attPrj[k][i]);
        divAtt[k][i][j] = float(exp(divAtt[k][i][j]));
      }
    }

    for(i=0; i<RAYS; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[RAYS-j+i];
    }
  }

  ///////////////////////////
  // filtering projections //
  ///////////////////////////

	double *SLFilter = new double[RAYS];
  double *tempVec  = new double[RAYS];
  double **attPrjCH = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjSH = CreateMatrix<double>(ANGLE, RAYS);
  float  ***divAttTemp = CreateVolume<float>(ANGLE, RAYS, RAYS);
  double cs, sn;

	for(i=0; i<RAYS; i++) SLFilter[i] = RampSL(i);

	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = attRadon[k][i]*cos(attPrjH[k][i])*exp(attPrj[k][i]);
      attPrjS[k][i] = attRadon[k][i]*sin(attPrjH[k][i])*exp(attPrj[k][i]);
    }
  }

	for(k=0; k<ANGLE; k++)
  {
    memcpy(tempVec, attPrjC[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = 0.0; attPrjCH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjC[k][i] += tempVec[j]*SLFilter[abs(i-j)];
        attPrjCH[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }
    if(filter >= 1) ProjectionSimulator::SavGol5(attPrjC[k], RAYS, tempVec);
    
    memcpy(tempVec, attPrjS[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjS[k][i] = 0.0; attPrjSH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjS[k][i] += tempVec[j]*SLFilter[abs(i-j)];
        attPrjSH[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }
    if(filter >= 1) ProjectionSimulator::SavGol5(attPrjS[k], RAYS, tempVec);
  }

	for(k=0; k<ANGLE; k++)
  {
		for(i=0; i<RAYS; i++)
    {
      cs = cos(attPrjH[k][i]);
      sn = sin(attPrjH[k][i]);
			for(j=0; j<RAYS; j++) 
      {
        divAttTemp[k][i][j] = float(divAtt[k][i][j]*cs);
			  divAtt[k][i][j] = float(divAtt[k][i][j]*sn);
      }
    }
  }

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

	double r, s, t, cy, mm, nn;
	double HALFGRID = SIZE/2.-0.5;
	double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double RECON_STEP = IMAGE_DIAMETER/SIZE;
	double PROJ_STEP = IMAGE_DIAMETER/RAYS;
  double sqr = HALFGRID*HALFGRID;
  double ratio = RECON_STEP/PROJ_STEP;
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP);
    sin_[i] = sin(i*ANGLE_STEP);
  }

  memset(IMAGE[0], 0, sizeof(double)*SIZE*SIZE);
  for(m=0; m<SIZE; m++)
  {
    mm = m-HALFGRID;
		for(n=0; n<SIZE; n++)
		{
      nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;

			for(k=0; k<ANGLE; k++)
			{
				s = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*ratio; j = int(s); cy = s-j; jj = int(s+0.5);
        t = (sin_[k]*mm + cos_[k]*nn + HALFGRID)*ratio; i = int(t+0.5);
        if(jj>0 && jj<RAYS-1)
        {
          IMAGE[m][n] += divAttTemp[k][jj][i]*(attPrjC[k][j]+cy*(attPrjC[k][j+1]-attPrjC[k][j]));
          IMAGE[m][n] += divAtt[k][jj][i]*(attPrjS[k][j]+cy*(attPrjS[k][j+1]-attPrjS[k][j]));
          IMAGE[m][n] += (divAttTemp[k][jj+1][i]-divAttTemp[k][jj-1][i])*(attPrjCH[k][j]+cy*(attPrjCH[k][j+1]-attPrjCH[k][j]))/2.;
          IMAGE[m][n] += (divAtt[k][jj+1][i]-divAtt[k][jj-1][i])*(attPrjSH[k][j]+cy*(attPrjSH[k][j+1]-attPrjSH[k][j]))/2.;
        }
			}
		}
  }

  double coef = 0.5/ANGLE;
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  ////////////////////////////////
  // releasing allocated memory //
  ////////////////////////////////

  FreeVolume(divAtt);
  FreeVolume(divAttTemp);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(attPrjC);
  FreeMatrix(attPrjS);
  FreeMatrix(attPrjCH);
  FreeMatrix(attPrjSH);
  delete[] hilbert;
  delete[] tempVec;
  delete[] SLFilter;
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::FBP_ParYJ1(double **attRadon,  // attenuated projections
                                 int ANGLE,          // viewing angle numbers
                                 int RAYS,           // projection bins in each view
                                 double **IMAGE,     // image to be reconstructed
                                 int SIZE,           // the size of reconstructed image
                                 double **Atten,     // the attenuation map, has to have the same size with projection bins
                                 int filter
                                 )
{
  ////////////////////////////////////
  // preprocessing attenuation map  //
  ////////////////////////////////////

  if(filter >= 2) ProjectionSimulator::FilterMedian(attRadon, ANGLE, RAYS);

	int    i, k, j, m, n;
	double ANGLE_STEP = 2*pi/ANGLE;
  float  ***divAtt = CreateVolume<float>(ANGLE, RAYS, RAYS);
  float  ***attKer = CreateVolume<float>(ANGLE, RAYS, RAYS);
  double **attPrj  = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjH = CreateMatrix<double>(ANGLE, RAYS);
  double **attRadonH = CreateMatrix<double>(ANGLE, RAYS);
  double **attCos = CreateMatrix<double>(ANGLE, RAYS);
  double **attSin = CreateMatrix<double>(ANGLE, RAYS);
  double **attExp = CreateMatrix<double>(ANGLE, RAYS);
  double *hilbert  = new double[2*RAYS+1];
  double *tempBuf  = new double[RAYS];
  double *tempBuf1  = new double[RAYS];

  for(i=0; i<2*RAYS+1; i++) hilbert[i] = HilbertSTD(i-RAYS);

  for(k=0; k<ANGLE; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, divAtt[k], RAYS, k*ANGLE_STEP);
		for(i=0; i<RAYS; i++)
    {
			for(j=RAYS-2; j>=0; j--) divAtt[k][i][j] += divAtt[k][i][j+1];
      attPrj[k][i] = divAtt[k][i][0]/2.;
		  for(j=0; j<RAYS; j++) divAtt[k][i][j] = float(exp(divAtt[k][i][j]));
    }

    for(i=0; i<RAYS; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[RAYS+i-j];
    }
  }
  for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      attExp[k][i] = exp(attPrj[k][i]);
      attCos[k][i] = cos(attPrjH[k][i]);
      attSin[k][i] = sin(attPrjH[k][i]);
    }
  }

  for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
      for(j=0; j<RAYS; j++) 
        attKer[k][i][j] = float(hilbert[i-j+RAYS]*(attCos[k][i]*attCos[k][j]+attSin[k][i]*attSin[k][j])*attExp[k][j]/attExp[k][i]);
  }

  ///////////////////////////
  // filtering projections //
  ///////////////////////////

	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      attRadonH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) attRadonH[k][i] += attKer[k][i][j]*attRadon[k][j];
      for(j=0; j<RAYS; j++) divAtt[k][i][j] = float(divAtt[k][i][j]*attRadonH[k][i]);
    }
    for(j=0; j<RAYS; j++)
    {
      for(i=0; i<RAYS; i++) tempBuf[i] = divAtt[k][i][j];
      if(filter >= 1) ProjectionSimulator::SavGol5(tempBuf, RAYS, tempBuf1);
      for(i=1; i<RAYS-1; i++) divAtt[k][i][j] = float(0.5*(tempBuf[i+1]-tempBuf[i-1]));
      divAtt[k][0][j] = divAtt[k][1][j];
      divAtt[k][RAYS-1][j] = divAtt[k][RAYS-2][j];
    }
  }

  FreeMatrix(attCos);
  FreeMatrix(attSin);
  FreeMatrix(attExp);
  FreeVolume(attKer);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(attRadonH);
  delete[] hilbert;
  delete[] tempBuf;
  delete[] tempBuf1;

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

	double r, s, t, cy, mm, nn;
	double HALFGRID = SIZE/2.-0.5;
	double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double RECON_STEP = IMAGE_DIAMETER/SIZE;
	double PROJ_STEP = IMAGE_DIAMETER/RAYS;
  double sqr = HALFGRID*HALFGRID;
  double ratio = RECON_STEP/PROJ_STEP;
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP);
    sin_[i] = sin(i*ANGLE_STEP);
  }

  memset(IMAGE[0], 0, sizeof(double)*SIZE*SIZE);
  for(m=0; m<SIZE; m++)
  {
    mm = m-HALFGRID;
		for(n=0; n<SIZE; n++)
		{
      nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;

			for(k=0; k<ANGLE; k++)
			{
				s = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*ratio; j = int(s); cy = s-j;
        t = (sin_[k]*mm + cos_[k]*nn + HALFGRID)*ratio; i = int(t+0.5);
        if(j>0 && j<RAYS-1)
        {
          IMAGE[m][n] += divAtt[k][j][i];// + cy*(divAtt[k][j+1][i] - divAtt[k][j][i]);
        }
			}
		}
  }

  double coef = 0.5/ANGLE;
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  ////////////////////////////////
  // releasing allocated memory //
  ////////////////////////////////

  FreeVolume(divAtt);
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::FBP_ParKL(double **attRadon,  // attenuated projections
                                int ANGLE,          // viewing angle numbers
                                int RAYS,           // projection bins in each view
                                double **IMAGE,     // image to be reconstructed
                                int SIZE,           // the size of reconstructed image
                                double **Atten,     // the attenuation map, has to have the same size with projection bins
                                int filter
                                )
{
  ////////////////////////////////////
  // preprocessing attenuation map  //
  ////////////////////////////////////

  if(filter >= 2) ProjectionSimulator::FilterMedian(attRadon, ANGLE, RAYS);

	int    i, k, j, m, n, jj;
	double ANGLE_STEP = 2*pi/ANGLE;
  double *hilbert  = new double[2*RAYS+1];
  float  ***divAttCos = CreateVolume<float>(ANGLE, RAYS, RAYS);
  float  ***divAttSin = CreateVolume<float>(ANGLE, RAYS, RAYS);
  double **attPrj  = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjH = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjC = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjS = CreateMatrix<double>(ANGLE, RAYS);

  for(i=0; i<=2*RAYS; i++) hilbert[i] = HilbertKL(i-RAYS);

  // step P-I: sign-convolution of attenuation
	for(k=0; k<ANGLE; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, divAttSin[k], RAYS, k*ANGLE_STEP);
		for(i=0; i<RAYS; i++)
    {
			for(j=RAYS-2; j>=0; j--) divAttSin[k][i][j] += divAttSin[k][i][j+1];
      attPrj[k][i] = divAttSin[k][i][0]/2.;
		  for(j=0; j<RAYS; j++) 
      {
        divAttSin[k][i][j] -= float(attPrj[k][i]);
        divAttSin[k][i][j] = float(exp(divAttSin[k][i][j]));
      }
    }
    for(i=0; i<RAYS; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[RAYS+i-j];
    }
  }

  ///////////////////////////
  // filtering projections //
  ///////////////////////////

  double *tempVec  = new double[RAYS];
  double cs, sn;

	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = attRadon[k][i]*cos(attPrjH[k][i])*exp(attPrj[k][i]);
      attPrjS[k][i] = attRadon[k][i]*sin(attPrjH[k][i])*exp(attPrj[k][i]);
    }
  }

	for(k=0; k<ANGLE; k++)
  {
    memcpy(tempVec, attPrjC[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjC[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }

    memcpy(tempVec, attPrjS[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjS[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjS[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }

		for(i=0; i<RAYS; i++)
    {
      cs = cos(attPrjH[k][i]);
      sn = sin(attPrjH[k][i]);
			for(j=0; j<RAYS; j++) 
      {
        divAttCos[k][i][j] = float(divAttSin[k][i][j]*cs*attPrjC[k][i]);
			  divAttSin[k][i][j] = float(divAttSin[k][i][j]*sn*attPrjS[k][i]);
      }
    }
  }

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

  memset(IMAGE[0], 0, sizeof(double)*SIZE*SIZE);

	double r, s, t, cy, mm, nn;
	double HALFGRID = SIZE/2.-0.5;
	double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double RECON_STEP = IMAGE_DIAMETER/SIZE;
	double PROJ_STEP = IMAGE_DIAMETER/RAYS;
  double sqr = HALFGRID*HALFGRID;
  double ratio = RECON_STEP/PROJ_STEP;

  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];  
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP);
    sin_[i] = sin(i*ANGLE_STEP);
  }

  for(m=0; m<SIZE; m++)
  {
    mm = m-HALFGRID;
		for(n=0; n<SIZE; n++)
		{
      nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;

			for(k=0; k<ANGLE; k++)
			{
				s = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*ratio; j = int(s); cy = s-j; jj = int(s+0.5);
        t = (sin_[k]*mm + cos_[k]*nn + HALFGRID)*ratio; i = int(t+0.5);
        if(jj>0 && jj<RAYS-1)
        {
          IMAGE[m][n] += (divAttCos[k][jj+1][i]-divAttCos[k][jj-1][i])/2.;
          IMAGE[m][n] += (divAttSin[k][jj+1][i]-divAttSin[k][jj-1][i])/2.;
        }
			}
		}
  }

  double coef = 0.5/ANGLE;
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  ////////////////////////////////
  // releasing allocated memory //
  ////////////////////////////////

  FreeVolume(divAttSin);
  FreeVolume(divAttCos);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(attPrjC);
  FreeMatrix(attPrjS);
  delete[] hilbert;
  delete[] tempVec;
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::FBP_Par180( double **attRadon,  // attenuated projections
                                  int ANGLE,          // viewing angle numbers
                                  int RAYS,           // projection bins in each view
                                  double **IMAGE,     // image to be reconstructed
                                  double **Atten,     // the attenuation map, has to have the same size with projection bins
                                  int filter,         // filter to smooth the projections
                                  int iter            // number of iterations
                                  )
{
  if(filter >= 2) ProjectionSimulator::FilterMedian(attRadon, ANGLE, RAYS);

	double ANGLE_STEP = pi/ANGLE;
  int HalfAng = ANGLE;
  ANGLE = 2*ANGLE;

  ////////////////////////////////////
  // preprocessing attenuation map  //
  ////////////////////////////////////
  //
	int    i, k, j, m, n, jj;
  double *hilbert  = new double[2*RAYS+1];
  float  ***divAtt = CreateVolume<float>(ANGLE, RAYS, RAYS);
  float  ***divAttCos = CreateVolume<float>(ANGLE, RAYS, RAYS);
  float  ***divAttSin = CreateVolume<float>(ANGLE, RAYS, RAYS);
  double **attPrj  = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjH = CreateMatrix<double>(ANGLE, RAYS);

  for(i=0; i<2*RAYS+1; i++) hilbert[i] = HilbertKL(i-RAYS);
  double intensity1 = 0.0;
  for(k=0; k<HalfAng; k++) for(i=0; i<RAYS; i++) intensity1 += attRadon[k][i];

  for(k=0; k<ANGLE; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, divAtt[k], RAYS, k*ANGLE_STEP);
		for(i=0; i<RAYS; i++)
    {
			for(j=RAYS-2; j>=0; j--) divAtt[k][i][j] += divAtt[k][i][j+1];
      attPrj[k][i] = divAtt[k][i][0]/2.;
		  for(j=0; j<RAYS; j++) 
      {
        divAtt[k][i][j] = float(exp(divAtt[k][i][j]-attPrj[k][i]));
      }
    }
    for(i=0; i<RAYS; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[RAYS-j+i];
    }
  }

  ///////////////////////////
  // filtering projections //
  ///////////////////////////
  //
  double **attPrjC = CreateMatrix<double>(ANGLE, RAYS);
  double **attPrjS = CreateMatrix<double>(ANGLE, RAYS);
  double *tempVec  = new double[RAYS];
  double cs, sn;

  for(k=0; k<HalfAng; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = attRadon[k][i]*cos(attPrjH[k][i])*exp(attPrj[k][i]);
      attPrjS[k][i] = attRadon[k][i]*sin(attPrjH[k][i])*exp(attPrj[k][i]);
    }
  }

  for(k=0; k<HalfAng; k++)
  {
    memcpy(tempVec, attPrjC[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjC[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjC[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }
    memcpy(tempVec, attPrjS[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      attPrjS[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        attPrjS[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }
  }

	for(k=0; k<HalfAng; k++)
  {
		for(i=0; i<RAYS; i++)
    {
      cs = cos(attPrjH[k][i]);
      sn = sin(attPrjH[k][i]);
			for(j=0; j<RAYS; j++) 
      {
        divAttCos[k][i][j] = float(divAtt[k][i][j]*cs*attPrjC[k][i]);
			  divAttSin[k][i][j] = float(divAtt[k][i][j]*sn*attPrjS[k][i]);
      }
    }
  }

  ////////////////////////////////
  // backprojecting projections //
  ////////////////////////////////

	double r, s, t, cy, mm, nn;
	double HALFGRID = RAYS/2.-0.5;
	double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double RECON_STEP = IMAGE_DIAMETER/RAYS;
	double PROJ_STEP = IMAGE_DIAMETER/RAYS;
  double sqr = HALFGRID*HALFGRID;
  double ratio = RECON_STEP/PROJ_STEP;
  
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP);
    sin_[i] = sin(i*ANGLE_STEP);
  }

  memset(IMAGE[0], 0, sizeof(double)*RAYS*RAYS);
  for(m=0; m<RAYS; m++)
  {
    mm = m-HALFGRID;
		for(n=0; n<RAYS; n++)
		{
      nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;

			for(k=0; k<HalfAng; k++)
			{
				s = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*ratio; j = int(s); cy = s-j; jj = int(s+0.5);
        t = (sin_[k]*mm + cos_[k]*nn + HALFGRID)*ratio; i = int(t+0.5);
        if(jj>0 && jj<RAYS-1)
        {
          if(filter && jj >= 3 && jj < RAYS-3)
          {
            IMAGE[m][n] += 0.286*(divAttCos[k][jj+1][i]-divAttCos[k][jj-1][i]) + 0.1715*(divAttCos[k][jj+2][i]-divAttCos[k][jj-2][i]) 
                         - 0.043*(divAttCos[k][jj+3][i]-divAttCos[k][jj-3][i]);
            IMAGE[m][n] += 0.286*(divAttSin[k][jj+1][i]-divAttSin[k][jj-1][i]) + 0.1715*(divAttSin[k][jj+2][i]-divAttSin[k][jj-2][i])
                         - 0.043*(divAttSin[k][jj-3][i]-divAttSin[k][jj+3][i]); 
          }
          else
          {
            IMAGE[m][n] += (divAttCos[k][jj+1][i]-divAttCos[k][jj-1][i])/2.;
            IMAGE[m][n] += (divAttSin[k][jj+1][i]-divAttSin[k][jj-1][i])/2.;
          }
        }
			}
		}
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // start the iteration
  //
  double intensity2 = 0.0;
  double **recnImg = CreateMatrix<double>(RAYS, RAYS);
  double **tempImg = CreateMatrix<double>(RAYS, RAYS);
  memcpy(recnImg[0], IMAGE[0], sizeof(double)*RAYS*RAYS);
  int iterations = 0;
  while( iterations < iter )
  {
    // calculate the intensity balance factor
    intensity2 = 0.0;
    for(k=0; k<HalfAng; k++)
    {
      ProjectionSimulator::RotateCarte(IMAGE, tempImg, RAYS, k*ANGLE_STEP);
      for(i=0; i<RAYS; i++) 
      {
        attPrjS[k][i] = 0.0;
        for(j=0; j<RAYS; j++) 
        {
          attPrjS[k][i] += tempImg[i][j]/divAtt[k][i][j];
        }
        intensity2 += attPrjS[k][i] * exp(-attPrj[k][i]);
      }
    }
    intensity2 = intensity1/intensity2;

    // filter the missing projections
    for(k=HalfAng; k<ANGLE; k++)
    {
      ProjectionSimulator::RotateCarte(IMAGE, tempImg, RAYS, k*ANGLE_STEP);
      for(i=0; i<RAYS; i++) 
      {
        attPrjS[k][i] = 0.0;
        for(j=0; j<RAYS; j++) 
        {
          attPrjS[k][i] += tempImg[i][j]/divAtt[k][i][j];
        }
        attPrjC[k][i] = attPrjS[k][i]*cos(attPrjH[k][i])*intensity2;
        attPrjS[k][i] = attPrjS[k][i]*sin(attPrjH[k][i])*intensity2;
      }
      memcpy(tempVec, attPrjC[k], sizeof(double)*RAYS);
      for(i=0; i<RAYS; i++)
      {
        attPrjC[k][i] = 0.0;
        for(j=0; j<RAYS; j++) 
        {
          attPrjC[k][i] += tempVec[j]*hilbert[RAYS-j+i];
        }
      }
      memcpy(tempVec, attPrjS[k], sizeof(double)*RAYS);
      for(i=0; i<RAYS; i++)
      {
        attPrjS[k][i] = 0.0;
        for(j=0; j<RAYS; j++) 
        {
          attPrjS[k][i] += tempVec[j]*hilbert[RAYS-j+i];
        }
      }

		  for(i=0; i<RAYS; i++)
      {
        s = cos(attPrjH[k][i]);
        t = sin(attPrjH[k][i]);
			  for(j=0; j<RAYS; j++) 
        {
          divAttCos[k][i][j] = float(divAtt[k][i][j]*s*attPrjC[k][i]);
			    divAttSin[k][i][j] = float(divAtt[k][i][j]*t*attPrjS[k][i]);
        }
      }
    }

    // backproject the missing projections
    memset(IMAGE[0], 0, sizeof(double)*RAYS*RAYS);
    for(m=0; m<RAYS; m++)
    {
      mm = m-HALFGRID;
		  for(n=0; n<RAYS; n++)
		  {
        nn = n-HALFGRID; r = (mm*mm+nn*nn); if(r >= sqr) continue;

			  for(k=HalfAng; k<ANGLE; k++)
			  {
				  s = (cos_[k]*mm - sin_[k]*nn + HALFGRID)*ratio; j = int(s); cy = s-j; jj = int(s+0.5);
          t = (sin_[k]*mm + cos_[k]*nn + HALFGRID)*ratio; i = int(t+0.5);
          if(jj>0 && jj<RAYS-1)
          {
            if(filter && jj >= 3 && jj < RAYS-3)
            {
              IMAGE[m][n] += 0.286*(divAttCos[k][jj+1][i]-divAttCos[k][jj-1][i]) + 0.1715*(divAttCos[k][jj+2][i]-divAttCos[k][jj-2][i]) 
                           - 0.043*(divAttCos[k][jj+3][i]-divAttCos[k][jj-3][i]);
              IMAGE[m][n] += 0.286*(divAttSin[k][jj+1][i]-divAttSin[k][jj-1][i]) + 0.1715*(divAttSin[k][jj+2][i]-divAttSin[k][jj-2][i])
                           - 0.043*(divAttSin[k][jj-3][i]-divAttSin[k][jj+3][i]); 
            }
            else
            {
              IMAGE[m][n] += (divAttCos[k][jj+1][i]-divAttCos[k][jj-1][i])/2.;
              IMAGE[m][n] += (divAttSin[k][jj+1][i]-divAttSin[k][jj-1][i])/2.;
            }
          }
			  }
		  }
    }

    for(i=0; i<RAYS; i++) for(j=0; j<RAYS; j++) IMAGE[i][j] += recnImg[i][j];
    iterations ++;
  }
  
  double coef = 0.5/ANGLE;
  for(m=0; m<RAYS; m++) for(n=0; n<RAYS; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, RAYS, RAYS);

  ////////////////////////////////
  // releasing allocated memory //
  ////////////////////////////////

  FreeVolume(divAtt);
  FreeVolume(divAttSin);
  FreeVolume(divAttCos);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(attPrjC);
  FreeMatrix(attPrjS);
  FreeMatrix(recnImg);
  FreeMatrix(tempImg);
  delete[] hilbert;
  delete[] tempVec;
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::FBP_FanAng( double **PROJECTION,
                                  int PROJECT_RAY,
                                  int ROTAT_ANGLE,
                                  double **IMAGE,
                                  int SIZE,
                                  int filter  )
{
	int i, j, k, m, n;
	double half_ray = PROJECT_RAY/2. - 0.5;
	double PROJECT_RAY_STEP = RAY_VIEW/PROJECT_RAY;
  double DIST_SOURCE = RADIUS_FOCA*SIZE/(2*RADIUS_IMAGE); // in pixel

  // begin the convolution procedure.	
	double *R = new double[PROJECT_RAY];
	double *AA = new double[PROJECT_RAY];
	double *BB = new double[PROJECT_RAY];
	double *FILTER = new double[PROJECT_RAY];

	AA[0] = 1.0;
	for(i=1; i<PROJECT_RAY; i++)
  {
		AA[i] = i*PROJECT_RAY_STEP/sin(i*PROJECT_RAY_STEP);
    AA[i] *= AA[i];
  }
	
	for(i=0; i<PROJECT_RAY; i++)
  {
		BB[i] = cos((i-half_ray)*PROJECT_RAY_STEP);
		FILTER[i] = RampSL(i); //Shepp-Logan filter
  }
	
	for(k=0; k<ROTAT_ANGLE; k++)
  {
		memcpy(R, PROJECTION[k], sizeof(double)*PROJECT_RAY);
		for(i=0; i<PROJECT_RAY; i++)
		{
			PROJECTION[k][i] = 0.0;
			for(n=0; n<PROJECT_RAY; n++)
			{	
				m = abs(i-n);
				PROJECTION[k][i] += R[n]*BB[n]*AA[m]*FILTER[m];
			}
    }
  }
	delete []AA; delete []BB; delete []FILTER; delete []R;

  //begin the backprojection procedure
	double a, b, r, e, aa, theta, c1, c2, c3;
	double half_size = SIZE/2. - 0.5; 
	double image_step = RADIUS_FOCA*sin(RAY_VIEW/2)/half_size;
  double sqr = half_size*half_size*image_step*image_step;
	double c = RADIUS_FOCA*RADIUS_FOCA;
	double ROTAT_ANGLE_STEP = (pi+pi)/ROTAT_ANGLE;
  double *cos_ = new double[ROTAT_ANGLE];
  double *sin_ = new double[ROTAT_ANGLE];
  for(i=0; i<ROTAT_ANGLE; i++)
  {
    cos_[i] = cos(i*ROTAT_ANGLE_STEP);
    sin_[i] = sin(i*ROTAT_ANGLE_STEP);
  }
	
  double beta, val;
  int betaIndex = 0;
	for(i=0; i<SIZE; i++)
  {
    a = image_step*(i-half_size);
		for(j=0; j<SIZE; j++)
		{
			b = -image_step*(j-half_size); 
      r = (a*a+b*b);
			if(r < sqr)
			{
				IMAGE[i][j] = 0.0;
				for(m=0; m<ROTAT_ANGLE; m++)
				{
          c2 = (sin_[m]*a - cos_[m]*b);
          c1 = c + r + 2*RADIUS_FOCA*c2;
					theta = atan((cos_[m]*a+sin_[m]*b)/(RADIUS_FOCA+c2));
					e = theta/PROJECT_RAY_STEP + half_ray; n = int(e); aa = e-n;
          //IMAGE[i][j] += (PROJECTION[m][n]+aa*(PROJECTION[m][n+1]-PROJECTION[m][n]))/c1;

          // this segment is for partial scan
          beta = m*ROTAT_ANGLE_STEP + 2*theta + pi; 
          betaIndex = int(beta/ROTAT_ANGLE_STEP + 0.5); 
          betaIndex %= ROTAT_ANGLE;
          val = PROJECTION[betaIndex][PROJECT_RAY-1-int(e+0.5)];
          c3 = 2*RADIUS_FOCA*c2;
          if( (val < -1.e-8 || val > 1.e-8) && (PROJECTION[m][n] < -1.e-8 || PROJECTION[m][n] > 1.e-8) )
            IMAGE[i][j] += 0.5*(PROJECTION[m][n]+aa*(PROJECTION[m][n+1]-PROJECTION[m][n]))/c1;
          else
            IMAGE[i][j] += (PROJECTION[m][n]+aa*(PROJECTION[m][n+1]-PROJECTION[m][n]))/c1*(0.5+0.5*sqrt((c1)/(c+r-c3)));
				}
			}
			else
				IMAGE[i][j] = 0.0;
		}
  }
  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  double coef = PROJECT_RAY*ROTAT_ANGLE_STEP/pi;
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  delete[] cos_; delete[] sin_;
}

void ReconAnalytical::FBP_FanHTT( double **PROJECTION,
                                  int PROJECT_RAY,
                                  int ROTAT_ANGLE,
                                  double **IMAGE,
                                  int SIZE,
                                  int partial)
{
	int i, j, k, m, n;
	double half_ray = PROJECT_RAY/2. - 0.5;
	double PROJECT_RAY_STEP = RAY_VIEW/PROJECT_RAY;
	double ROTAT_ANGLE_STEP = (pi+pi)/ROTAT_ANGLE;
	double ratio = PROJECT_RAY_STEP/ROTAT_ANGLE_STEP;

	double *R = new double[PROJECT_RAY];
	double *hilbertA = new double[2*PROJECT_RAY+1];
	for(i=0; i<=2*PROJECT_RAY; i++) hilbertA[i] = HilbertAng(i-PROJECT_RAY, PROJECT_RAY_STEP);

	// calculate the partial derivatives
	double** proj_f = CreateMatrix<double>(ROTAT_ANGLE, PROJECT_RAY);

	if( partial )
	{
		for(k=1; k<ROTAT_ANGLE-1; k++)
		{
			for(i=1; i<PROJECT_RAY-1; i++)
				proj_f[k][i] = (PROJECTION[k][i+1]-PROJECTION[k][i-1])-(PROJECTION[k+1][i]-PROJECTION[k-1][i])*ratio;
		}
		for(i=1; i<PROJECT_RAY-1; i++)
			proj_f[0][i] = (PROJECTION[1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;
		for(i=1; i<PROJECT_RAY-1; i++)
			proj_f[ROTAT_ANGLE-1][i] = (PROJECTION[ROTAT_ANGLE-1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;

		if(partial == 1) // over scan
			memset(proj_f[ROTAT_ANGLE*5/8], 0, ROTAT_ANGLE/4*PROJECT_RAY*sizeof(double));
		else if(partial == 2) // short scan
			memset(proj_f[ROTAT_ANGLE*7/12], 0, ROTAT_ANGLE/3*PROJECT_RAY*sizeof(double));
		else if(partial == 3) // ROI scan 1
		{
			memset(proj_f[ROTAT_ANGLE/12], 0, ROTAT_ANGLE/6*PROJECT_RAY*sizeof(double));
			memset(proj_f[ROTAT_ANGLE*5/12], 0, ROTAT_ANGLE/6*PROJECT_RAY*sizeof(double));
			memset(proj_f[ROTAT_ANGLE*3/4], 0, ROTAT_ANGLE/6*PROJECT_RAY*sizeof(double));
		}
		else if(partial == 4) // ROI scan 2
		{
			memset(proj_f[0], 0, ROTAT_ANGLE/12*PROJECT_RAY*sizeof(double));
			memset(proj_f[ROTAT_ANGLE*5/12], 0, ROTAT_ANGLE*7/12*PROJECT_RAY*sizeof(double));
		}
	}
	else
	{
		for(k=0; k<ROTAT_ANGLE; k++)
		{
			for(i=1; i<PROJECT_RAY-1; i++)
				proj_f[k][i] = (PROJECTION[k][i+1]-PROJECTION[k][i-1]);
			proj_f[k][0] = proj_f[k][1];
			proj_f[k][PROJECT_RAY-1] = proj_f[k][PROJECT_RAY-2];
		}
	}
	memcpy(PROJECTION[0], proj_f[0], sizeof(double)*ROTAT_ANGLE*PROJECT_RAY);

	// begin the angular Hilbert transform	
	for(k=0; k<ROTAT_ANGLE; k++)
	{
		memcpy(R, PROJECTION[k], sizeof(double)*PROJECT_RAY);
		for(i=0; i<PROJECT_RAY; i++)
		{
			PROJECTION[k][i] = 0.0;
			for(n=0; n<PROJECT_RAY; n++)
				PROJECTION[k][i] += R[n]*hilbertA[PROJECT_RAY+i-n];
		}
	}

	delete[] hilbertA;
	FreeMatrix(proj_f);
	delete []R;

	//begin the backprojection procedure
	double a, b, r, e, aa, sigma, c1;
	double half_size = SIZE/2. - 0.5; 
	double image_step = RADIUS_FOCA*sin(RAY_VIEW/2)/half_size;
	double sqr = half_size*half_size*image_step*image_step;
	double c = RADIUS_FOCA*RADIUS_FOCA;
	double *cos_ = new double[ROTAT_ANGLE];
	double *sin_ = new double[ROTAT_ANGLE];
	for(i=0; i<ROTAT_ANGLE; i++)
	{
		cos_[i] = cos(i*ROTAT_ANGLE_STEP);
		sin_[i] = sin(i*ROTAT_ANGLE_STEP);
	}

	double beta = 0.;
	double v1, v2;
	int betaIndex = 0;
	double eps = 1.e-12;
	double rcos, rsin, k2, delta;
	for(i=0; i<SIZE; i++)
	{
		a = (i-half_size)*image_step; 
		for(j=0; j<SIZE; j++)
		{
			b = -(j-half_size)*image_step; 
			r = (a*a+b*b);
			if(r < sqr)
			{
				IMAGE[i][j] = 0.0;
				for(m=0; m<ROTAT_ANGLE; m++)
				{
					rcos = (cos_[m]*a+sin_[m]*b);
					rsin = (sin_[m]*a - cos_[m]*b);
					c1 = RADIUS_FOCA + rsin;
					sigma = atan(rcos/c1);
					k2 = (c1+rsin)*RADIUS_FOCA + r;
					delta = 1./sqrt(k2);
					e = sigma/PROJECT_RAY_STEP + half_ray; n = int(e); aa = e-n;

					beta = m*ROTAT_ANGLE_STEP + sigma + sigma + pi; 
					betaIndex = int(beta/ROTAT_ANGLE_STEP + 0.5); 
					betaIndex %= ROTAT_ANGLE;
					k = PROJECT_RAY-1-n;
					v1 = PROJECTION[m][n]+aa*(PROJECTION[m][n+1]-PROJECTION[m][n]);
					v2 = PROJECTION[betaIndex][k];
					if( (v1 < -eps || v1 > eps) && (v2 < -eps || v2 > eps) )
						IMAGE[i][j] += 0.5*v1*delta;
					else
						IMAGE[i][j] += v1*delta;
				}
			}
			else
				IMAGE[i][j] = 0.0;
		}
	}
	delete[] cos_; delete[] sin_;

	ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
	double coef = PROJECT_RAY/(2.*ROTAT_ANGLE);
	for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;
}

void ReconAnalytical::FBP_FanAHT( double **PROJECTION,
                                  int PROJECT_RAY,
                                  int ROTAT_ANGLE,
                                  double **IMAGE,
                                  int SIZE,
								                  int partial)
{
	int i, j, k, m, n;
	double half_ray = PROJECT_RAY/2. - 0.5;
	double PROJECT_RAY_STEP = RAY_VIEW/PROJECT_RAY;
	double ROTAT_ANGLE_STEP = (pi+pi)/ROTAT_ANGLE;
	double ratio = PROJECT_RAY_STEP/ROTAT_ANGLE_STEP;

  // begin the angular Hilbert transform	
	double *R = new double[PROJECT_RAY];
	double  *hilbertA = new double[2*PROJECT_RAY+1];
	for(i=0; i<=2*PROJECT_RAY; i++) hilbertA[i] = HilbertAng(i-PROJECT_RAY, PROJECT_RAY_STEP);

	for(k=0; k<ROTAT_ANGLE; k++)
	{
		memcpy(R, PROJECTION[k], sizeof(double)*PROJECT_RAY);
		for(i=0; i<PROJECT_RAY; i++)
		{
			PROJECTION[k][i] = 0.0;
			for(n=0; n<PROJECT_RAY; n++)
				PROJECTION[k][i] += R[n]*hilbertA[PROJECT_RAY+i-n];
		}
	}
	delete[] hilbertA;

	// calculate the partial derivatives
	double** proj_f = CreateMatrix<double>(ROTAT_ANGLE, PROJECT_RAY);

	if(partial)
	{
		for(k=1; k<ROTAT_ANGLE-1; k++)
		{
			for(i=1; i<PROJECT_RAY-1; i++)
				proj_f[k][i] = (PROJECTION[k][i+1]-PROJECTION[k][i-1])-(PROJECTION[k+1][i]-PROJECTION[k-1][i])*ratio;
		}
		for(i=1; i<PROJECT_RAY-1; i++)
			proj_f[0][i] = (PROJECTION[1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;
		for(i=1; i<PROJECT_RAY-1; i++)
			proj_f[ROTAT_ANGLE-1][i] = (PROJECTION[ROTAT_ANGLE-1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;
	}
	else
	{
		for(k=0; k<ROTAT_ANGLE; k++)
		{
			for(i=1; i<PROJECT_RAY-1; i++)
				proj_f[k][i] = (PROJECTION[k][i+1]-PROJECTION[k][i-1]);
			proj_f[k][0] = proj_f[k][1];
			proj_f[k][PROJECT_RAY-1] = proj_f[k][PROJECT_RAY-2];
		}
	}

	// partial scan processing
	if(partial == 1) // over scan
		memset(proj_f[ROTAT_ANGLE*5/8], 0, ROTAT_ANGLE/4*PROJECT_RAY*sizeof(double));
	else if(partial == 2) // short scan
		memset(proj_f[ROTAT_ANGLE*7/12], 0, ROTAT_ANGLE/3*PROJECT_RAY*sizeof(double));
	else if(partial == 3) // ROI scan 1
	{
		memset(proj_f[ROTAT_ANGLE/9], 0, ROTAT_ANGLE/9*PROJECT_RAY*sizeof(double));
		memset(proj_f[ROTAT_ANGLE*4/9], 0, ROTAT_ANGLE/9*PROJECT_RAY*sizeof(double));
		memset(proj_f[ROTAT_ANGLE*7/9], 0, ROTAT_ANGLE/9*PROJECT_RAY*sizeof(double));
	}
	else if(partial == 4) // ROI scan 2
	{
		memset(proj_f[0], 0, ROTAT_ANGLE/12*PROJECT_RAY*sizeof(double));
		memset(proj_f[ROTAT_ANGLE*5/12], 0, ROTAT_ANGLE*7/12*PROJECT_RAY*sizeof(double));
	}

	// weighting filtered projections
	for(i=0; i<PROJECT_RAY; i++) R[i] = RADIUS_FOCA*cos((i-half_ray)*PROJECT_RAY_STEP);
	for(i=0; i<PROJECT_RAY; i++) for(k=0; k<ROTAT_ANGLE; k++) proj_f[k][i] /= R[i];

	memcpy(PROJECTION[0], proj_f[0], sizeof(double)*ROTAT_ANGLE*PROJECT_RAY);
	FreeMatrix(proj_f);
	delete []R;

	//begin the backprojection procedure
	double a, b, r, e, aa, bb, sigma;
	double half_size = SIZE/2. - 0.5; 
	double image_step = RADIUS_FOCA*sin(RAY_VIEW/2)/half_size;
	double sqr = half_size*half_size*image_step*image_step;
	double c = RADIUS_FOCA*RADIUS_FOCA;
	double *cos_ = new double[ROTAT_ANGLE];
	double *sin_ = new double[ROTAT_ANGLE];
	for(i=0; i<ROTAT_ANGLE; i++)
	{
		cos_[i] = cos(i*ROTAT_ANGLE_STEP);
		sin_[i] = sin(i*ROTAT_ANGLE_STEP);
	}

	int betaIndex = 0, kk;
	double eps = 1.e-12;
	double v1, v2, rcos, beta = 0.;
	for(i=0; i<SIZE; i++)
	{
		a = (i-half_size)*image_step; 
		for(j=0; j<SIZE; j++)
		{
			b = -(j-half_size)*image_step; 
			r = (a*a+b*b);
			if(r < sqr)
			{
				IMAGE[i][j] = 0.0;
				for(m=0; m<ROTAT_ANGLE; m++)
				{
					rcos = (cos_[m]*a+sin_[m]*b);
					sigma = asin(rcos/RADIUS_FOCA);
					e = sigma/PROJECT_RAY_STEP + half_ray; n = int(e); aa = e-n;
					beta = m*ROTAT_ANGLE_STEP - sigma + 2*pi;
					e = beta/ROTAT_ANGLE_STEP; k = int(e); bb = e-k; k %= ROTAT_ANGLE;

					beta = m*ROTAT_ANGLE_STEP + sigma + pi;
					betaIndex = int(beta/ROTAT_ANGLE_STEP + 0.5); 
					betaIndex %= ROTAT_ANGLE;
					kk = PROJECT_RAY-1-n;
					if(k < ROTAT_ANGLE-1)
						v1 = PROJECTION[k][n]+aa*(PROJECTION[k][n+1]-PROJECTION[k][n])+bb*(PROJECTION[k+1][n]-PROJECTION[k][n])+
							aa*bb*(PROJECTION[k+1][n+1]+PROJECTION[k][n]-PROJECTION[k+1][n]-PROJECTION[k][n+1]);
					else
						v1 = PROJECTION[k][n]+aa*(PROJECTION[k][n+1]-PROJECTION[k][n]);
					v2 = PROJECTION[betaIndex][kk];
					if( (v1 < -eps || v1 > eps) && (v2 < -eps || v2 > eps) )
						IMAGE[i][j] += 0.5*v1;
					else
						IMAGE[i][j] += v1;
				}
			}
			else
				IMAGE[i][j] = 0.0;
		}
	}
	delete[] cos_; delete[] sin_;

	ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
	double coef = PROJECT_RAY/2./ROTAT_ANGLE;
	for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;
}

void ReconAnalytical::FBP_FanSpc( double **PROJECTION,
				                          int    PROJECT_RAY,
				                          int    ROTAT_ANGLE,
				                          double **IMAGE,
				                          int    SIZE,
                                  int    filter)
{
  if(filter == 1) ProjectionSimulator::FilterMedian(PROJECTION, ROTAT_ANGLE, PROJECT_RAY);

	int i, j, k, m, n;
	double a, b, dd;
	double half_ray = PROJECT_RAY/2.-0.5;
	double PROJ_STEP = RAY_FOV/PROJECT_RAY;

  // begin the convolution procedure.	
	double *R, *AA, *FILTER;
	R = new double[PROJECT_RAY];
	AA = new double[PROJECT_RAY];
	FILTER = new double[PROJECT_RAY];
	
	a = RADIUS_FOCA+RADIUS_DETC; 
  dd = a;
  a *= a;
  dd *= a;
	for(i=0; i<PROJECT_RAY; i++)
	{
		b = (i-half_ray)*PROJ_STEP; b *= b;
		AA[i] = dd/sqrt(a+b);
		FILTER[i] = RampSL(i); //Shepp-Logan filter
	}
	
	for(k=0; k<ROTAT_ANGLE; k++)
  {
		memcpy(R, PROJECTION[k], sizeof(double)*PROJECT_RAY);
		for(i=0; i<PROJECT_RAY; i++)
		{
			PROJECTION[k][i] = 0.0;
			for(n=0; n<PROJECT_RAY; n++) PROJECTION[k][i] += R[n]*AA[n]*FILTER[abs(i-n)];
    }
    if(filter == 1) ProjectionSimulator::SavGol5(PROJECTION[k], PROJECT_RAY, R);
  }
  delete []AA; delete []FILTER; delete []R;

  // begin the backprojection procedure
	double c, r, d, e, aa, c1;
	double half_size = SIZE/2. - 0.5;	
  double sqr = half_size*half_size;
  double IMAGE_DIAMETER = RADIUS_IMAGE*2;
	double ANGLE_STEP = (pi+pi)/ROTAT_ANGLE;
	double IMAGE_STEP = IMAGE_DIAMETER/SIZE*(RADIUS_FOCA+RADIUS_DETC)/RADIUS_FOCA;

  double *cos_ = new double[ROTAT_ANGLE];
  double *sin_ = new double[ROTAT_ANGLE];
  for(i=0; i<ROTAT_ANGLE; i++)
  {
    cos_[i] = cos(i*ANGLE_STEP);
    sin_[i] = sin(i*ANGLE_STEP);
  }

	c = RADIUS_FOCA+RADIUS_DETC;
	for(i=0; i<SIZE; i++)
  {
		for(j=0; j<SIZE; j++)
		{
			a = i-half_size; b = j-half_size; r = (a*a+b*b);
      a *= IMAGE_STEP; b *= -IMAGE_STEP;
			if(r < sqr)
			{
				IMAGE[i][j] = 0.0;
				for(m=0; m<ROTAT_ANGLE; m++)
				{
					c1 = 1./(c+sin_[m]*a - cos_[m]*b);
          d = c*c1*(cos_[m]*a + sin_[m]*b);
          c1 *= c1;
					e = d/PROJ_STEP + half_ray; n = int(e); aa = e-n;
					IMAGE[i][j] += (PROJECTION[m][n]+aa*(PROJECTION[m][n+1]-PROJECTION[m][n]))*c1;
				}
			}
			else
				IMAGE[i][j] = 0.0;
		}
  }
  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  double coef = 1.732*0.25/ROTAT_ANGLE;
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  delete[] cos_; delete[] sin_;
}

void ReconAnalytical::FBP_FanAng( double **fanRadon,
                                  int ANGLE,
                                  int RAYS,
                                  double **IMAGE,
                                  int SIZE,
                                  double **Atten,
                                  int attSize,
                                  double **fanAtt,
                                  int filter
                                  )
{
  // ============================= //
  // preprocessing attenuation map //
  // ============================= //

	int     i, k, j, m, n, kk, ii, jj;
  int     attAngle = ANGLE;
  double  attAngleStep = 2*pi/attAngle;
  double  attRayStep = RAY_VIEW/RAYS;
  float   ***attSgn = CreateVolume<float>(attAngle, attSize, attSize);
	float ***attSgnT = CreateVolume<float>(attAngle, attSize, attSize);
  double  **attPrj = CreateMatrix<double>(attAngle, attSize);
  double  **attPrjH = CreateMatrix<double>(attAngle, attSize);
  double  **fanAttH = CreateMatrix<double>(ANGLE, RAYS);
  double  *hilbert = new double[2*attSize+1];
  double  *hilbertA = new double[2*RAYS+1];
  double cs, sn;
  
  for(i=0; i<=2*attSize; i++) hilbert[i] = HilbertKL(i-attSize);
  for(i=0; i<=2*RAYS; i++) hilbertA[i] = HilbertAng(i-RAYS, attRayStep);

  // step P-I: Sign convolution of attenuation map
	for(k=0; k<attAngle; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, attSgn[k], attSize, k*attAngleStep);
		for(i=0; i<attSize; i++)
    {
			for(j=attSize-2; j>=0; j--) attSgn[k][i][j] += attSgn[k][i][j+1];
      attPrj[k][i] = attSgn[k][i][0]/2.;
		  for(j=0; j<attSize; j++) attSgn[k][i][j] = float(exp(attSgn[k][i][j]-attPrj[k][i]));
    }
  }

  // step P-II: Hilbert transform of attenuation projections
	for(k=0; k<attAngle; k++)
  {
    for(i=0; i<attSize; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<attSize; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[attSize-j+i];
    }
  }

	for(k=0; k<ANGLE; k++)
  {
		for(i=0; i<attSize; i++)
    {
      cs = cos(attPrjH[k][i]); 
      sn = sin(attPrjH[k][i]);
			for(j=0; j<attSize; j++) 
      {
        attSgnT[k][i][j] = float(attSgn[k][i][j]*cs);
        attSgn[k][i][j] = float(attSgn[k][i][j]*sn);
      }
    }
  }

  // step P-III: angular Hilbert transform of attenuation projections
	for(k=0; k<ANGLE; k++)
  {
    for(j=0; j<RAYS; j++) fanAtt[k][j] /= 2;
    for(i=0; i<RAYS; i++)
    {
      //using interpolation to calculate the angular Hilbert transform
      //kk = int((k*attAngleStep + (i-RAYS/2)*attRayStep)/attAngleStep+0.5);
      //kk = (kk+ANGLE)%ANGLE;
      //ii = int(RADIUS_FOCA*sin((i-RAYS/2)*attRayStep)*RAYS/(RADIUS_DETC*2)+attSize/2+0.5);
      //fanAttH[k][i] = attPrjH[kk][ii]/(2*pi);
      //fanAtt[k][i] = attPrj[kk][ii];

      // using fan-beam projected attenuation to calculate the angular Hilbert transform
      fanAttH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) fanAttH[k][i] += fanAtt[k][j]*hilbertA[RAYS-j+i]/(2*pi);  // this is a mystery constant
    }
  }

  // ======================//
  // filtering projections //
  // ======================//

  // step P-IV: weighting the attenuated projections
  double **fanRadonC = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonS = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonCH = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonSH = CreateMatrix<double>(ANGLE, RAYS);
  double *tempVec = new double[RAYS];
	double *SLFilter = new double[RAYS];

	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      fanRadonC[k][i] = fanRadon[k][i]*cos(fanAttH[k][i])*exp(fanAtt[k][i])*RADIUS_FOCA*cos(attRayStep*(i-RAYS/2.+0.5));
      fanRadonS[k][i] = fanRadon[k][i]*sin(fanAttH[k][i])*exp(fanAtt[k][i])*RADIUS_FOCA*cos(attRayStep*(i-RAYS/2.+0.5));
    }
  }
  if(filter == 1)
  {
    ProjectionSimulator::FilterMedian(fanRadonC, ANGLE, RAYS);
    ProjectionSimulator::FilterMedian(fanRadonS, ANGLE, RAYS);
  }

  SLFilter[0] = RampSL(0);
  double attRayST; 
	for(i=1; i<RAYS; i++) 
  {
    attRayST = i*attRayStep;
    SLFilter[i] = RampSL(i)*attRayST*attRayST;
    attRayST = sin(attRayST);
    SLFilter[i] /= (attRayST*attRayST);
  }

	for(k=0; k<ANGLE; k++)
  {
    memcpy(tempVec, fanRadonC[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      fanRadonC[k][i] = 0.0; 
      fanRadonCH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        fanRadonC[k][i] += tempVec[j]*SLFilter[abs(i-j)];
        fanRadonCH[k][i] += tempVec[j]*hilbertA[RAYS-j+i];
      }
    }
    if(filter == 1)
      ProjectionSimulator::SavGol5(fanRadonC[k], RAYS, tempVec);

    memcpy(tempVec, fanRadonS[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      fanRadonS[k][i] = 0.0; fanRadonSH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        fanRadonS[k][i] += tempVec[j]*SLFilter[abs(i-j)];
        fanRadonSH[k][i] += tempVec[j]*hilbertA[RAYS-j+i];
      }
    }
    if(filter == 1)
      ProjectionSimulator::SavGol5(fanRadonS[k], RAYS, tempVec);
  }

  // ===========================//
  // backprojecting projections //
  // ===========================//
	double a, b, r, e, aa, sigma, theta, c1, c2, ss, tt, ct, st;
	double half_size = SIZE/2. - 0.5;
  double half_ray = RAYS/2. - 0.5;
  double sqr = half_size*half_size;
  double DIST_SOURCE = RADIUS_FOCA*SIZE/(2*RADIUS_IMAGE); // in pixel
  double dist2 = 2*DIST_SOURCE;
	double c = DIST_SOURCE*DIST_SOURCE;
	double image_step = DIST_SOURCE*sin(RAY_VIEW/2)/half_size;
  double *tempAttVec = new double[attSize];

  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*attAngleStep);
    sin_[i] = sin(i*attAngleStep);
  }

  memset(IMAGE[0], 0, sizeof(double)*SIZE*SIZE);
	for(i=0; i<SIZE; i++)
  {
		for(j=0; j<SIZE; j++)
		{
      a = i-half_size; b = j-half_size; r = (a*a+b*b); if( r >= sqr ) continue;
      a *= image_step; b *= -image_step;

      for(m=0; m<ANGLE; m++)
			{
        c2 = (sin_[m]*a - cos_[m]*b);
        c1 = c + r + dist2*c2; //c1 *= c1;
				sigma = atan((cos_[m]*a+sin_[m]*b)/(DIST_SOURCE+c2));
				e = sigma/attRayStep + half_ray; n = int(e); aa = e-n;

        theta = sigma + m*attAngleStep; ct = cos(theta); st = sin(theta);
        kk = int(theta/attAngleStep + 0.5); kk = (kk+ANGLE)%ANGLE;
				ss = (ct*a - st*b)/image_step + half_size; jj = int(ss+0.5); 
        tt = (st*a + ct*b)/image_step + half_size; ii = int(tt+0.5); 
        
        if( jj >= 1 && jj < attSize-1 )
        {
          IMAGE[i][j] += attSgnT[kk][jj][ii]*(fanRadonC[m][n]+aa*(fanRadonC[m][n+1]-fanRadonC[m][n]))/c1;
          IMAGE[i][j] += attSgn[kk][jj][ii]*(fanRadonS[m][n]+aa*(fanRadonS[m][n+1]-fanRadonS[m][n]))/c1;
          c1 = sqrt(c1)*SIZE*2;
          IMAGE[i][j] += (attSgnT[kk][jj+1][ii]-attSgnT[kk][jj-1][ii])*(fanRadonCH[m][n]+aa*(fanRadonCH[m][n+1]-fanRadonCH[m][n]))/c1;
          IMAGE[i][j] += (attSgn[kk][jj+1][ii]-attSgn[kk][jj-1][ii])*(fanRadonSH[m][n]+aa*(fanRadonSH[m][n+1]-fanRadonSH[m][n]))/c1;
        }
			}
		}
  }
  double coef = RAYS*RAYS/(4.*ANGLE);
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  // ===========================//
  // releasing allocated memory //
  // ===========================//
  FreeVolume(attSgn);
  FreeVolume(attSgnT);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(fanAttH);
  FreeMatrix(fanRadonC);
  FreeMatrix(fanRadonS);
  FreeMatrix(fanRadonCH);
  FreeMatrix(fanRadonSH);
  delete[] hilbert;
  delete[] hilbertA;
  delete[] SLFilter;
  delete[] tempVec;
  delete[] tempAttVec;
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::FBP_FanSpc( double **fanRadon,
                                  int ANGLE,
                                  int RAYS,
                                  double **IMAGE,
                                  int SIZE,
                                  double **Atten,
                                  int attSize,
                                  double **fanAtt,
                                  int filter
                                  )
{
  // ============================= //
  // preprocessing attenuation map //
  // ============================= //

	int     i, k, j, m, n, kk, ii, jj;
  double  attAngleStep = 2*pi/ANGLE;
  double  prjFanStep = RAY_FOV/RAYS*RADIUS_FOCA/(RADIUS_FOCA + RADIUS_DETC);
  double  prjParStep = 2*RADIUS_DETC/attSize;
  float   ***attSgn = CreateVolume<float>(ANGLE, attSize, attSize);
	float  ***attSgnT = CreateVolume<float>(ANGLE, attSize, attSize);
  double  **attPrj = CreateMatrix<double>(ANGLE, attSize);
  double  **attPrjH = CreateMatrix<double>(ANGLE, attSize);
  double  **fanAttH = CreateMatrix<double>(ANGLE, RAYS);
  double  *hilbert = new double[2*attSize+1];
  double  u, uu, cs, sn;
  double  D = RADIUS_FOCA, DD = D*D;
  
  for(i=0; i<=2*attSize; i++) hilbert[i] = HilbertKL(i-attSize);

  // step P-I: Sign convolution of attenuation map
	for(k=0; k<ANGLE; k++)
  {
    ProjectionSimulator::RotateCarte(Atten, attSgn[k], attSize, k*attAngleStep);
		for(i=0; i<attSize; i++)
    {
			for(j=attSize-2; j>=0; j--) attSgn[k][i][j] += attSgn[k][i][j+1];
      attPrj[k][i] = attSgn[k][i][0]/2.;
		  for(j=0; j<attSize; j++) attSgn[k][i][j] = float(exp(attSgn[k][i][j]-attPrj[k][i]));
    }

    for(i=0; i<attSize; i++)
    {
      attPrjH[k][i] = 0.0;
      for(j=0; j<attSize; j++) attPrjH[k][i] += attPrj[k][j]*hilbert[attSize-j+i];
    }

		for(i=0; i<attSize; i++)
    {
      cs = cos(attPrjH[k][i]); 
      sn = sin(attPrjH[k][i]);
			for(j=0; j<attSize; j++) 
      {
        attSgnT[k][i][j] = float(attSgn[k][i][j]*cs);
        attSgn[k][i][j] = float(attSgn[k][i][j]*sn);
      }
    }
  }

  // angular Hilbert transform of attenuation projections
	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      //using interpolation to calculate the angular Hilbert transform
      u = (i-RAYS/2.+0.5)*prjFanStep; uu = u*u;
      kk = int(k + atan(u/D)/attAngleStep + ANGLE + 0.5);
      kk = kk%ANGLE;
      ii = int(D*u/sqrt(DD+uu)/prjFanStep + attSize/2 + 0.5);
      fanAtt[k][i] = attPrj[kk][ii];
      fanAttH[k][i] = attPrjH[kk][ii]/(2*pi);
    }
  }

  // ======================//
  // filtering projections //
  // ======================//

  // step P-IV: weighting the attenuated projections
  double **fanRadonC = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonS = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonCH = CreateMatrix<double>(ANGLE, RAYS);
  double **fanRadonSH = CreateMatrix<double>(ANGLE, RAYS);
  double *tempVec = new double[RAYS];
	double *SLFilter = new double[RAYS];

	for(i=0; i<RAYS; i++) SLFilter[i] = RampSL(i);

	for(k=0; k<ANGLE; k++)
  {
    for(i=0; i<RAYS; i++)
    {
      u = (i-RAYS/2.+0.5)*prjFanStep; uu = u*u; u = D/sqrt(DD+uu);
      fanRadonC[k][i] = fanRadon[k][i]*cos(fanAttH[k][i])*exp(fanAtt[k][i])*u;
      fanRadonS[k][i] = fanRadon[k][i]*sin(fanAttH[k][i])*exp(fanAtt[k][i])*u;
    }
  }
  if(filter == 1)
  {
    ProjectionSimulator::FilterMedian(fanRadonC, ANGLE, RAYS);
    ProjectionSimulator::FilterMedian(fanRadonS, ANGLE, RAYS);
  }

	for(k=0; k<ANGLE; k++)
  {
    memcpy(tempVec, fanRadonC[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      fanRadonC[k][i] = 0.0; 
      for(j=0; j<RAYS; j++) fanRadonC[k][i] += tempVec[j]*SLFilter[abs(i-j)];
    }
    for(j=0; j<RAYS; j++) 
    {
      u = (j-RAYS/2.+0.5)*prjFanStep; uu = u*u;
      tempVec[j] *= D/sqrt(DD+uu);
    }
    for(i=0; i<RAYS; i++)
    {
      fanRadonCH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) fanRadonCH[k][i] += tempVec[j]*hilbert[RAYS-j+i];
    }
    if(filter == 1)
      ProjectionSimulator::SavGol5(fanRadonC[k], RAYS, tempVec);

    memcpy(tempVec, fanRadonS[k], sizeof(double)*RAYS);
    for(i=0; i<RAYS; i++)
    {
      fanRadonS[k][i] = 0.0;
      for(j=0; j<RAYS; j++) fanRadonS[k][i] += tempVec[j]*SLFilter[abs(i-j)];
    }
    for(j=0; j<RAYS; j++) 
    {
      u = (j-RAYS/2.+0.5)*prjFanStep; uu = u*u;
      tempVec[j] *= D/sqrt(DD+uu);
    }
    for(i=0; i<RAYS; i++)
    {
      fanRadonSH[k][i] = 0.0;
      for(j=0; j<RAYS; j++) 
      {
        fanRadonSH[k][i] += tempVec[j]*hilbert[RAYS-j+i];
      }
    }
    if(filter == 1)
      ProjectionSimulator::SavGol5(fanRadonS[k], RAYS, tempVec);
  }

  // ===========================//
  // backprojecting projections //
  // ===========================//
	double a, b, r, aa, sigma, theta, ss, tt, ct, st, w;
	double half_size = SIZE/2. - 0.5;
	double image_step = 2*RADIUS_IMAGE/SIZE;
  double half_ray = RAYS/2. - 0.5;
  double sqr = half_size*half_size;

  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(i=0; i<ANGLE; i++)
  {
    cos_[i] = cos(i*attAngleStep);
    sin_[i] = sin(i*attAngleStep);
  }

  memset(IMAGE[0], 0, sizeof(double)*SIZE*SIZE);
	for(i=0; i<SIZE; i++)
  {
		for(j=0; j<SIZE; j++)
		{
      a = i-half_size; b = j-half_size; 
      r = (a*a+b*b); if( r >= sqr ) continue;
      a *= image_step; b *= -image_step;

      for(m=0; m<ANGLE; m++)
			{
				w = D/(D+sin_[m]*a - cos_[m]*b);
        u = w*(cos_[m]*a + sin_[m]*b);
				sigma = atan(u/D);
        u = u/prjFanStep + half_ray; n = int(u); aa = u-n;
        
        theta = sigma + m*attAngleStep;
        kk = int(theta/attAngleStep + ANGLE + 0.5); kk = kk%ANGLE;

        ct = cos(theta); st = sin(theta);
				ss = (ct*a - st*b)/image_step + half_size; jj = int(ss+0.5); 
        tt = (st*a + ct*b)/image_step + half_size; ii = int(tt+0.5); 

        if( jj >= 1 && jj < attSize-1 )
        {
          //w *= w;
          IMAGE[i][j] += attSgnT[kk][jj][ii]*(fanRadonC[m][n]+aa*(fanRadonC[m][n+1]-fanRadonC[m][n]))*w;
          IMAGE[i][j] += attSgn[kk][jj][ii]*(fanRadonS[m][n]+aa*(fanRadonS[m][n+1]-fanRadonS[m][n]))*w;

          w = sqrt(w)/(2*SIZE);
          IMAGE[i][j] += (attSgnT[kk][jj+1][ii]-attSgnT[kk][jj-1][ii])*(fanRadonCH[m][n]+aa*(fanRadonCH[m][n+1]-fanRadonCH[m][n]))*w;
          IMAGE[i][j] += (attSgn[kk][jj+1][ii]-attSgn[kk][jj-1][ii])*(fanRadonSH[m][n]+aa*(fanRadonSH[m][n+1]-fanRadonSH[m][n]))*w;
        }
			}
		}
  }

  double coef = RAYS*RAYS/(4.*ANGLE);
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;

  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  // ===========================//
  // releasing allocated memory //
  // ===========================//
  FreeVolume(attSgn);
  FreeVolume(attSgnT);
  FreeMatrix(attPrj);
  FreeMatrix(attPrjH);
  FreeMatrix(fanAttH);
  FreeMatrix(fanRadonC);
  FreeMatrix(fanRadonS);
  FreeMatrix(fanRadonCH);
  FreeMatrix(fanRadonSH);
  delete[] hilbert;
  delete[] SLFilter;
  delete[] tempVec;
  delete[] cos_;
  delete[] sin_;
}

void ReconAnalytical::DBF_FanAng( double **PROJECTION,
                                  int PROJECT_RAY, 
                                  int ROTAT_ANGLE, 
                                  double **IMAGE, 
                                  int SIZE, 
                                  double mu)
{
	int i, j, k, m, n;
	double half_ray = PROJECT_RAY/2. - 0.5;
	double PROJECT_RAY_STEP = RAY_VIEW/PROJECT_RAY;
	double ROTAT_ANGLE_STEP = (pi+pi)/ROTAT_ANGLE;
	double ratio = PROJECT_RAY_STEP/ROTAT_ANGLE_STEP;

	// calculate the partial derivatives
	double** proj_f = CreateMatrix<double>(ROTAT_ANGLE, PROJECT_RAY);

	for(k=1; k<ROTAT_ANGLE-1; k++)
	{
		for(i=1; i<PROJECT_RAY-1; i++)
			proj_f[k][i] = (PROJECTION[k][i+1]-PROJECTION[k][i-1])-(PROJECTION[k+1][i]-PROJECTION[k-1][i])*ratio;
	}
	for(i=1; i<PROJECT_RAY-1; i++)
		proj_f[0][i] = (PROJECTION[1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;
	for(i=1; i<PROJECT_RAY-1; i++)
		proj_f[ROTAT_ANGLE-1][i] = (PROJECTION[ROTAT_ANGLE-1][i+1]-PROJECTION[ROTAT_ANGLE-1][i-1]) - (PROJECTION[1][i]-PROJECTION[ROTAT_ANGLE-1][i])*ratio;
	
	// weighting filtered projections
  double *R = new double[PROJECT_RAY];
	for(i=0; i<PROJECT_RAY; i++) R[i] = RADIUS_FOCA*cos((i-half_ray)*PROJECT_RAY_STEP);
	for(i=0; i<PROJECT_RAY; i++) for(k=0; k<ROTAT_ANGLE; k++) proj_f[k][i] /= R[i];
	memcpy(PROJECTION[0], proj_f[0], sizeof(double)*ROTAT_ANGLE*PROJECT_RAY);

  FreeMatrix(proj_f);
	delete []R;

  // backprojection process
	const double half_size = SIZE/2. - 0.5; 
	const double image_step = RADIUS_FOCA*sin(RAY_VIEW/2)/half_size;
	double *cos_ = new double[ROTAT_ANGLE];
	double *sin_ = new double[ROTAT_ANGLE];
	for(i=0; i<ROTAT_ANGLE; i++)
	{
		cos_[i] = cos(i*ROTAT_ANGLE_STEP);
		sin_[i] = sin(i*ROTAT_ANGLE_STEP);
	}

	double eps = 1.e-12;
	double x, y, rcos, rsin, sigma, beta;
  double sigmaW, betaW, temp;
	for(i=0; i<SIZE; i++)
	{
		x = (i-half_size)*image_step; 
		for(j=0; j<SIZE; j++)
		{
			y = -(j-half_size)*image_step; 

      IMAGE[i][j] = 0.0;
			for(m=0; m<ROTAT_ANGLE/2; m++)
			{
				rcos = (cos_[m]*x+sin_[m]*y);
        rsin = (sin_[m]*x-cos_[m]*y);

				sigma = asin(rcos/RADIUS_FOCA);
				temp = sigma/PROJECT_RAY_STEP + half_ray; if(temp <= 0.0 || temp >= (PROJECT_RAY - 0.5)) continue;
        n = int(temp); sigmaW = temp-n;

				beta = m*ROTAT_ANGLE_STEP - sigma; if(beta < 0.0)  beta += 2*pi;
				temp = beta/ROTAT_ANGLE_STEP; 
        k = int(temp); betaW = temp-k; k %= ROTAT_ANGLE;
        if(k == ROTAT_ANGLE-1 || n == PROJECT_RAY-1) 
        {
          IMAGE[i][j] += (exp(rsin*mu)*PROJECTION[k][n]);
          continue;
        }
        temp = sigmaW*(PROJECTION[k+1][n+1]+PROJECTION[k][n]-PROJECTION[k+1][n]-PROJECTION[k][n+1]);
        temp = betaW*(temp + PROJECTION[k+1][n]-PROJECTION[k][n]);
				temp += PROJECTION[k][n]+sigmaW*(PROJECTION[k][n+1]-PROJECTION[k][n]);

        IMAGE[i][j] += (exp(rsin*mu)*temp);
			}
		}
	}
	delete[] cos_; delete[] sin_;

	//ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
	//double coef = PROJECT_RAY/2./ROTAT_ANGLE;
	//for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;
}

void ReconAnalytical::CHD_FanAng( double **RADON,
                                  int RAY_ANGLE, 
                                  int ANGLE, 
                                  double **IMAGE, 
                                  int SIZE, 
                                  double mu)
{
	int j, k, m, n;
	double a, b, r, s, delta, theta, da, c, d;
	double RAY_STEP = RAY_VIEW/RAY_ANGLE;
	double ANGLE_STEP = 2*pi/ANGLE;

  double HALFRAY = RAY_ANGLE/2.-0.5;
	int half_ray = (RAY_ANGLE+1)/2;
	int half_angle = ANGLE/2;

	double **RADON1 = CreateMatrix<double>(RAY_ANGLE, ANGLE);
	double *R = new double[2*ANGLE];
	double *_COS = new double[half_angle];
	double *_SIN = new double[half_angle];
	
  // rebinning the projection using Fourier series shift
	for(k=0; k<RAY_ANGLE; k++)
	{
		for(j=0;j<ANGLE;j++)
		{
			R[2*j] = RADON[k][j]; R[2*j+1] = 0.0;
		}
		FFT(R, ANGLE, 1);
		
		theta = (k - HALFRAY)*RAY_STEP;
		for(j=0; j<half_angle; j++)
		{
			_COS[j] = cos(j*theta);
			_SIN[j] = sin(j*theta);
		}
		for(j=0; j<half_angle; j++) RADON[k][j] = R[2*j]*_COS[j]-R[2*j+1]*_SIN[j];
		for(j=0; j<half_angle; j++) RADON1[k][j] = R[2*j]*_SIN[j]+R[2*j+1]*_COS[j];
		for(j=half_angle; j<ANGLE; j++) RADON[k][j] = R[2*j]*_COS[ANGLE-j-1]+R[2*j+1]*_SIN[ANGLE-j-1];
		for(j=half_angle; j<ANGLE; j++) RADON1[k][j] = -R[2*j]*_SIN[ANGLE-j-1]+R[2*j+1]*_COS[ANGLE-j-1];
	}
	delete []_COS; delete []_SIN;

  // uneven sampling integral
	double **RADON2 = CreateMatrix<double>(half_ray, ANGLE*2);
  memset(RADON2[0], 0, sizeof(double)*half_ray*ANGLE*2);
  if( mu < eps ) mu = eps;
  double atten = pi/mu;
  double band = (RADIUS_FOCA+RADIUS_IMAGE)*RAY_STEP;
  double atten2;
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(k=0; k<ANGLE; k++)
  {
    cos_[k] = cos(k*ANGLE_STEP);
    sin_[k] = sin(k*ANGLE_STEP);
  }

	for(m=half_ray; m<RAY_ANGLE; m++)
	{	
		r = RADIUS_FOCA*sin((m-HALFRAY)*RAY_STEP);
    atten2 = mu*r;
		for(j=0; j<RAY_ANGLE; j++)
		{	
			s = RADIUS_FOCA*sin((j-HALFRAY)*RAY_STEP);
			delta = (sin((j-HALFRAY+0.5)*RAY_STEP)-sin((j-HALFRAY-0.5)*RAY_STEP)); 
			for(k=0; k<ANGLE; k++)
			{
				R[2*k] = (RampSL(r*cos_[k]-s, band) - RampSL(r*cos_[k]-s, atten))*exp(atten2*sin_[k]);
				R[2*k+1] = 0.0;
			}
			
			FFT(R, ANGLE, -1);
			for(k=0; k<ANGLE; k++)
			{
				RADON2[m-half_ray][2*k] += (R[2*k]*RADON[j][k] - R[2*k+1]*RADON1[j][k])*delta;
				RADON2[m-half_ray][2*k+1] += (R[2*k]*RADON1[j][k] + R[2*k+1]*RADON[j][k])*delta;
			}
		}
	}

  for(j=half_ray; j<RAY_ANGLE; j++)
	{
		FFT(RADON2[j-half_ray], ANGLE, -1);
		for(k=0; k<ANGLE; k++)
		{
			RADON[j][k] = RADON2[j-half_ray][2*k];
		}
	}
	delete []R;
  delete[] cos_;
  delete[] sin_;
  FreeMatrix(RADON1);
	FreeMatrix(RADON2);

  // coordinate transform
	double half_size = SIZE/2.-0.5;
  double ratio = 2.*RADIUS_IMAGE/(RADIUS_FOCA*SIZE);
	for(j=0; j<SIZE; j++)
	{
		for(k=0; k<SIZE; k++)
		{
			a = j - half_size;
			b = half_size - k;
			r = sqrt(a*a+b*b);
			if(r > half_size)
      {
        IMAGE[j][k] = 0.0;
        continue;
      }
      if(b >= 0.0)
				theta = acos(a/r);
			else 
				theta = pi+pi-acos(a/r);
			s = theta/ANGLE_STEP; n = int(s); d = s-n;
			da = asin(r*ratio); s = da/RAY_STEP; s += half_ray; m = int(s); c = s-m;
			
			if(m < RAY_ANGLE-1)
			{
				if(n < ANGLE-1)
					IMAGE[j][k]=
						RADON[m][n] + c*(RADON[m+1][n]-RADON[m][n]) +
						d*(RADON[m][n+1]-RADON[m][n]+c*(RADON[m+1][n+1]-RADON[m][n+1]+RADON[m][n]-RADON[m+1][n]));
				else
				{
					IMAGE[j][k]=
						RADON[m][n]+c*(RADON[m+1][n]-RADON[m][n]) + 
						d*(RADON[m][0]-RADON[m][n]+c*(RADON[m+1][0]-RADON[m][0]+RADON[m][n]-RADON[m+1][n]));
				}
			}
			else if(m==0)
			{
				if(n < ANGLE-1)
					IMAGE[j][k]= RADON[m][n] + d*(RADON[m][n+1]-RADON[m][n]);
				else
				{
					IMAGE[j][k]= RADON[m][n] + d*(RADON[m][0]-RADON[m][n]);
				}	
			}
		}
	}
  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);
}

void ReconAnalytical::CHD_FanSpc(	double **PROJECTION,
                                  int    RAYS,
                                  int    ANGLE,
                                  double **IMAGE,
                                  int    SIZE,
                                  double attenuat
                                )
{
	int    j, k, m, n;
	double a, b, r, s, delta, theta, c, d, step, tilt;
	double ANGLE_STEP = 2*pi/ANGLE;
	double RAY_STEP = RAY_FOV/RAYS;
  double HALFRAYGRID = RAYS/2. - 0.5;
	int half_angle = ANGLE/2;
	int half_ray = (RAYS+1)/2;
  
  if( attenuat < eps ) attenuat = eps;
	
	double *R = new double[2*ANGLE];
	double **PROJECTION1 = CreateMatrix<double>(RAYS, ANGLE);
	double *AA = new double[half_angle];
	double *BB = new double[half_angle];
	
	for(k=0; k<RAYS; k++)
	{
		for(j=0; j<ANGLE; j++)
		{
			R[2*j] = PROJECTION[k][j]; R[2*j+1] = 0.0;
		}
		FFT(R, ANGLE, 1);

		step = (k-HALFRAYGRID)*RAY_STEP;
		tilt = step/(RADIUS_FOCA+RADIUS_IMAGE);
		theta = atan(tilt);
		for(j=0; j<half_angle; j++)
		{
			AA[j] = sin(j*theta);
			BB[j] = cos(j*theta);
		}

  // rebinning projection by fourier shift
		for(j=0; j<half_angle; j++) PROJECTION[k][j] = R[j+j]*BB[j]-R[j+j+1]*AA[j];
		for(j=0; j<half_angle; j++) PROJECTION1[k][j] = R[j+j]*AA[j]+R[j+j+1]*BB[j];
		for(j=half_angle; j<ANGLE; j++) PROJECTION[k][j] = R[j+j]*BB[ANGLE-j-1]+R[j+j+1]*AA[ANGLE-j-1];
		for(j=half_angle; j<ANGLE; j++) PROJECTION1[k][j] = -R[j+j]*AA[ANGLE-j-1]+R[j+j+1]*BB[ANGLE-j-1];
	}
	delete[] AA; delete[] BB;

//================= integral calculation=================
	double **PROJECTION2 = CreateMatrix<double>(half_ray, 2*ANGLE);
  memset(PROJECTION2[0], 0, sizeof(double)*RAYS*ANGLE);

	double delta1, delta2, step1, step2;
	double ds = (RADIUS_FOCA+RADIUS_IMAGE)*(RADIUS_FOCA+RADIUS_IMAGE);
	double coeff = RADIUS_FOCA/(RADIUS_FOCA+RADIUS_IMAGE);
	double band = coeff*RAY_STEP;
  double atten = pi/attenuat;
  double atten2;
  double *cos_ = new double[ANGLE];
  double *sin_ = new double[ANGLE];
  for(k=0; k<ANGLE; k++)
  {
    cos_[k] = cos(k*ANGLE_STEP);
    sin_[k] = sin(k*ANGLE_STEP);
  }

	for(m=0; m<half_ray; m++)
	{	
		step = (m+0.5)*RAY_STEP;
		tilt = step/sqrt(step*step+ds);
		r = RADIUS_FOCA*tilt;
    atten2 = attenuat*r;
											
		for(j=0; j<RAYS; j++)
		{	
			step = (j-HALFRAYGRID)*RAY_STEP;
			tilt = step/sqrt(step*step+ds);
			s = RADIUS_FOCA*tilt;

			step1 = (j-HALFRAYGRID-0.5)*RAY_STEP;
			tilt = step1/sqrt(step1*step1+ds);
			delta1 = RADIUS_FOCA*tilt;

			step2 = (j-HALFRAYGRID+0.5)*RAY_STEP;
			tilt = step2/sqrt(step2*step2+ds);
			delta2 = RADIUS_FOCA*tilt;
			
			delta = delta2-delta1;

			for(k=0; k<ANGLE; k++)
			{
				R[k+k] = (RampSL(r*cos_[k]-s, band)-RampRect(r*cos_[k]-s, atten))*exp(atten2*sin_[k]);
				R[k+k+1] = 0.0;
			}
			
			FFT(R, ANGLE, -1);
			for(k=0; k<ANGLE; k++)
			{
				PROJECTION2[m][k+k] += (R[k+k]*PROJECTION[j][k]-R[k+k+1]*PROJECTION1[j][k])*delta;
				PROJECTION2[m][k+k+1] += (R[k+k]*PROJECTION1[j][k]+R[k+k+1]*PROJECTION[j][k])*delta;
			}
		}
	}

	for(j=0;j<half_ray;j++)
	{		
		FFT(PROJECTION2[j], ANGLE, -1);
		for(k=0;k<ANGLE;k++)
		{
			PROJECTION[j][k] = PROJECTION2[j][k+k]/(ANGLE);
		}
	}

  FreeMatrix(PROJECTION1);
  FreeMatrix(PROJECTION2);
	delete[] R; delete[] cos_, delete[] sin_;

///////// coordinate transform from polar to Cartesian
	double half_size = SIZE/2. - 0.5;
	double RECON_IMAGE_STEP = RAY_FOV*RADIUS_FOCA/sqrt(ds+RAY_FOV*RAY_FOV/4.0)/RAYS;
	for(j=0; j<SIZE; j++)
	{
		for(k=0; k<SIZE; k++)
		{
			a = j - half_size;
			b = half_size - k;
			r = sqrt(a*a+b*b);
			if(r > half_size)
      {
				IMAGE[j][k]=0.0;
        continue;
      }

      if(b >= 0.0)
				theta = acos(a/r);
			else 
				theta = pi+pi-acos(a/r);
      s = theta/(ANGLE_STEP);
			n = int(s); d = s - n;

			r = r*RECON_IMAGE_STEP;
			s = r*(RADIUS_FOCA+RADIUS_IMAGE)/sqrt(RADIUS_FOCA*RADIUS_FOCA-r*r);
			s = s/RAY_STEP; m = int(s); c = s-m;
			
			if(m < half_ray-1)
			{
				if(n<ANGLE-1)
					IMAGE[j][k]=
						PROJECTION[m][n]+c*(PROJECTION[m+1][n]-PROJECTION[m][n])+d*(PROJECTION[m][n+1]-PROJECTION[m][n]
						+c*(PROJECTION[m+1][n+1]-PROJECTION[m][n+1]+PROJECTION[m][n]-PROJECTION[m+1][n]));
				else
				{
					IMAGE[j][k]=
						PROJECTION[m][n]+c*(PROJECTION[m+1][n]-PROJECTION[m][n])+d*(PROJECTION[m][0]-PROJECTION[m][n]
						+c*(PROJECTION[m+1][0]-PROJECTION[m][0]+PROJECTION[m][n]-PROJECTION[m+1][n]));
				}
			}
			else
			{
				m = half_ray-1;
				if((n < ANGLE-1))
					IMAGE[j][k]= PROJECTION[m][n]+d*(PROJECTION[m][n+1]-PROJECTION[m][n]);
				else
				{
					IMAGE[j][k]= PROJECTION[m][n]+d*(PROJECTION[m][0]-PROJECTION[m][n]);
				}	
			}
		}
	}
  ProjectionSimulator::Smooth5(IMAGE, SIZE, SIZE);

  double coef = band*RAY_FOV/(ANGLE*8);
  for(m=0; m<SIZE; m++) for(n=0; n<SIZE; n++) IMAGE[m][n] *= coef;
}

//////////////////////////////////////////
// Filters used in image reconstruction //
//////////////////////////////////////////

double ReconAnalytical::HilbertSTD(int n)
{
  if(n == 0) return 0.;
  else return (1./(n*pi));
}

double ReconAnalytical::HilbertKL(int n)
{
  if(n >= 2) return ((n*log(1.-1./(n*n))+log((n+1.)/(n-1)))/pi);
  else if(n <= -2) return ((n*log(1.-1./(n*n))-log((n-1.)/(n+1)))/pi);
  else return (2*n*ln2/pi);
}

double ReconAnalytical::HilbertAttCos(double s, double step, double mu)
{
  if( fabs(s) < eps*step ) return 0.0;  
  return (cos(mu*s)/(pi*s));
}

double ReconAnalytical::HilbertRect(double s, double d)
{
  if( fabs(s) < eps*d ) return 0.0;
  else { double bb = s*pi; return (1.0 - cos(bb/d))/bb; }
}

double ReconAnalytical::HilbertCos(double s, double d)
{
  if( fabs(s) < eps*d ) return 0.0;  
  double bb = s*pi/d; 
  return ((1.0+sin(bb))/(pi*(2*s+d))+(1.0-sin(bb))/(pi*(2*s-d)));
}

double ReconAnalytical::HilbertAng(int n, double step)
{
  if( n == 0 ) return 0.0;
  return (step/(pi*sin(n*step)));
}

double ReconAnalytical::HilbertKer(double s, double step)
{
  if( s < step/4 && s > -step/4) return 0.0;
  return 1.0/s;
}

double ReconAnalytical::HilbertHyperCos(double s, double step, double mu)
{
  if( fabs(s) < eps*step ) return 0.0;
  else { double bb = exp(mu*s); bb = 0.5*(bb+1./bb); return (bb/s); }
}

double ReconAnalytical::HyperCosDiff(double s, double step, double mu)
{
  if( s < step/4 && s > -step/4) return 0.0;
  return ((exp(s*mu)+exp(-s*mu))/2.-1.0)/s;
}

double ReconAnalytical::HilbertHyperCosInv(double s, double step, double mu, double low, double up, double s1)
{
  if( fabs(s) < eps*step ) return 0.0;
  else 
  { 
    double bb = exp(mu*s); 
    bb = 0.5*s*(bb+1./bb);
    double aa = s*s;
    s1 *= -(aa*mu*mu);
    double sqUp = up*up;
    double sqLow = low*low;
    if(aa < sqLow) aa = sqrt(sqUp-aa)-sqrt(sqLow-aa);
    else if(aa < sqUp) aa = sqrt(sqUp-aa);
    else aa = 0.0;
    return (step*aa*(1+s1)/bb); 
  }
}

double ReconAnalytical::RampRect(double s, double d)
{
  if( fabs(s) < eps*d ) return (pi/(2*d*d));
	double filter, a= pi*s/d;
  filter = pi*sin(a)/(a*d*d) + (cos(a)-1.0)/(pi*s*s);
	return filter;
}

double ReconAnalytical::RampCos(double s, double d)
{
  double aa, bb, cc, dd;
  aa = s*pi/d;
  bb = d + 2*s;
  cc = d - 2*s;
  dd = sin(aa);
  if(bb > 0. && bb < eps*d) bb = eps*d;
  if(bb < 0. && bb > -eps*d) bb = -eps*d;
  if(cc > 0. && cc < eps*d) cc = eps*d;
  if(cc < 0. && cc > -eps*d) cc = -eps*d;

  double filter;
  filter = cos(aa)*(1./(d*bb)+1./(d*cc));
  filter -= (2*(1+dd)/(pi*bb*bb) +2*(1-dd)/(pi*cc*cc));
  return filter;
}

double ReconAnalytical::RampSL(double s, double d)
{
  double filter, bb, cc, a = sin(pi*s/d);
  bb = d + 2*s;
  cc = d - 2*s;
  if(bb > 0. && bb < eps*d) bb = eps*d;
  if(bb < 0. && bb > -eps*d) bb = -eps*d;
  if(cc > 0. && cc < eps*d) cc = eps*d;
  if(cc < 0. && cc > -eps*d) cc = -eps*d;

  filter = (1+a)/bb + (1-a)/cc;
  return (2*filter/(pi*d));
}

double ReconAnalytical::RampSL( double s, double d, double atten)
{
	double a, b, a1, b1;
	a = d+s+s;
	b = d-s-s;
	a1 = pi*atten*a;
	b1 = pi*atten*b;
	return (cos(a1)/a+cos(b1)/b)*(d*d/pi);
}

double ReconAnalytical::RampRect(int n)
{
  if( n == 0 ) return pi/2;
  else return (cos(n*pi)-1)/(pi*n*n);
}

double ReconAnalytical::RampCos(int n)
{
  double v = cos(n*pi)/(1-4*n*n);
  v -= 1./((1-2*n)*(1-2*n)*pi);
  v -= 1./((1+2*n)*(1+2*n)*pi);  
  return (2*v);
}

double ReconAnalytical::RampSL(int n)
{
  return (1.0/((0.25 - n*n)*pi));
}

double ReconAnalytical::Sinc(double s)
{
  s *= pi;
  if(s > eps || s < -eps) return sin(s)/s;
  else return 1.0;
}

void ReconAnalytical::FFT(double data[] , int nn, int isign)
{
  int ii, jj, n, mmax, m, j, i, istep;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr, tempi;

  n = 2*nn; // the size of array which is pointed by pointer data is 2*nn .
  j = 0;

// here is the part for re-ordering data
  for (ii=0; ii<nn; ii++)
  {
    i = 2*ii;
    if (j > i) //change two complexs
    {
      tempr = data[j];
      tempi = data[j+1];
      data[j] = data[i];
      data[j+1] = data[i+1];
      data[i] = tempr;
      data[i+1] = tempi;
    }
    m = n/2;
    while( m >= 2  && j > m-1 )
    {
      j = j-m;
      m = m/2;
    }
    j = j+m;
  }

// here is the part for Danielson-Lanczos
  mmax = 2;
  while( n > mmax )  // log(n) times repeatetion
  {
    istep = 2*mmax;
    theta = 6.28318530717959/(isign*mmax);
    wpr = cos(theta);
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (ii=1; ii<=mmax/2; ii++)
    {
      m = 2*ii - 1;
      for(jj=0; jj<=(n-m)/istep; jj++)
      {
        i =  m + jj*istep;
        j = i + mmax;
        tempr = wr*(data[j-1])-wi*(data[j]);
        tempi = wr*(data[j])+ wi*(data[j-1]);
        data[j-1] = data[i-1] - tempr;
        data[j] = data[i] - tempi;
        data[i-1] = data[i-1] + tempr;
        data[i] = data[i] + tempi;
      }
      wtemp = wr;
      wr = wr*wpr - wi*wpi;
      wi = wi*wpr + wtemp*wpi;
    }
    mmax = istep ;
  }
}
