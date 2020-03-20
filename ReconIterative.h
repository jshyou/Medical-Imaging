// ReconIterative.h: interface for the ReconIterative class.
//
/**
 * Copyright (c) 2003 - 2005  Cubic Imaging LLC.
 *
 * Permission to copy, modify and distribute the software and its documentation 
 * for noncommercial use is hereby granted without fee, provided that the above 
 * copyright notice appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation. Cubic Imaging LLC 
 * makes no representations about the suitability of this software for any purpose. 
 * It is provided "as is" without express or implied warranty.
 * 
 * Permission to copy, modify, and sell the software and its documentation in commercial use 
 * for profits has to be approved by written notice from Cubic Imaging LLC.
 * 
 * \author Jason You (jyou@cubic-imaging.com)
 * \date created at 12/01/2003
 *
 * ML-EM algorithm is derived using the statistic theory, and provides a general framework to handle
 * some difficult reconstruction problems when considering the noise filtering, attenuation compensation, 
 * PSF compensation and compalicated scanning geometries. The major drawback of ML-EM algorithm 
 * is the cost due to huge amount of computations compared with the FBP algorithm. The implementation 
 * here is based on rotation transformation to realize the projection and backprojection.
 * The implementation can be configured to run for both the EM and OS-EM algorithm. Currently, 
 * the parallel- and fan-beam data acquisition geometry are included. The attenuation is also considered.
 *
 * The ART algorithm also can be implemented using rotation transform, and will be added in the future.
 *
 **/

#ifndef _CUBICIMAGING_IMAGERECON_ITERATIVE_H_
#define _CUBICIMAGING_IMAGERECON_ITERATIVE_H_

//! class for iterative reconstruction algorithms
class ReconIterative  
{
public:
  ReconIterative();
  ~ReconIterative();

public:
  //! OS-EM algorithm for the parallel-beam 180 projection data
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of views for scanning the object image
    \param bins   - number of bins at each view, also the dimension of reconstructed image
    \param image  - square matrix to hold the reconstructed image
    \param step   - the interval length to regroup the index in OS-EM
    \param iter   - number of iterations to be performed

    The function implements the standard OS-EM algorithm.
  */
  void EM_Parall( double **Radon,     // projections
                  int views,          // viewing angle numbers in 180 degree
                  int bins,           // projection bins in each view
                  double **image,     // image to be reconstructed
                  int step,           // the interval length to regroup the index
                  int iter = 1        // number of iterations to be performed
                  );

  //! OS-EM algorithm for the attenuated parallel-beam 360 projection data
  /*!
    \param Radon  - 2D array to hold the attenuated projection data
    \param views  - number of views for scanning the object image
    \param bins   - number of bins at each view, also the dimension of the reconstructed image and attenuation map
    \param image  - square matrix of the reconstructed image
    \param atten  - square matrix of the attenuation map
    \param step   - the interval length to regroup the index in OS-EM
    \param iter   - number of iterations to be performed

    The function implements the OS-EM algorithm with attenuation compensation.
  */
  void EM_Parall( double **Radon,     // attenuated projections
                  int views,          // viewing angle numbers in 360 degree
                  int bins,           // projection bins in each view
                  double **image,     // image to be reconstructed
                  double **atten,     // attenuation map matrix
                  int step,           // the interval length to regroup the index
                  int iter = 1,       // number of iterations to be performed
                  bool half = false   // indication of the half scanning
                  );

  //! OS-EM algorithm for fan-beam projections with equiangular sampling
  /*!
    \param Radon  - 2D array for projections
    \param views  - number of views for scanning the object image
    \param rays   - number of bins at each fan view
    \param image  - square matrix of the reconstructed image
    \param dim    - dimension of the reconstructed image
    \param step   - the interval length to regroup the index in OS-EM
    \param iter   - number of iterations to be performed

    The function implements the OS-EM algorithm for fan-beam geometry with equiangular sampling.
  */
  void EM_FanAng( double **Radon,     // projections
                  int views,          // viewing angle numbers in 360 degree 
                  int rays,           // projection rays at each fan
                  double **image,     // image to be reconstructed 
                  int dim,            // dimension of the reconstructed image
                  int step,           // the interval length to regroup the index 
                  int iter = 1        // number of iterations to be performed
                  );

  //! OS-EM algorithm for attenuated fan-beam projections with equiangular sampling
  /*!
    \param Radon  - 2D array for projections
    \param views  - number of views for scanning the object image
    \param rays   - number of bins at each fan view
    \param image  - square matrix of the reconstructed image
    \param atten  - square matrix of the attenuation map
    \param dim    - dimension of the reconstructed image and attenuation map
    \param step   - the interval length to regroup the index in OS-EM
    \param iter   - number of iterations to be performed

    The function implements the OS-EM algorithm for attenuated fan-beam projections with euiangular sampling.
  */
  void EM_FanAng( double **Radon,     // attenuated projections 
                  int views,          // viewing angle numbers in 360 degree
                  int rays,           // projection rays at each fan 
                  double **image,     // image to be reconstructed 
                  int dim,            // dimension of the reconstructed image 
                  double **atten,     // attenuation map matrix 
                  int step,           // the interval length to regroup the index 
                  int iter = 1        // number of iterations to be performed
                  );

private:
  //! projector used in OS-EM for parallel-beam projections
  void project(double **image, int bins, double **Radon, int views, int start, int step);

  //! projector used in OS-EM for attenuated parallel-beam projections
  void project(double **image, double **atten, int bins, double **Radon, int views, int start, int step, double **att, double **img);

  //! backprojector for interleaved parallel-beam projections
  void backprj(double **Radon, int views, int start, int step, double **image, int bins);

  //! single view backprojector for parallel-beam projections
  void backprj(double  *Radon, int view, int bins, double **image);

  //! single view projector for fan-beam projections with equiangular sampling
  void project_fanA(double *Radon, int RAY, double **IMAGE, int SIZE, double pos, double **pImage);

  //! single view projector for attenuated fan-beam projections with equiangular sampling
  void project_fanA(double *Radon, int RAY, double **atten, double **IMAGE, int SIZE, double pos, double **pImage, double **pAtten);

  //! single view backprojector for fan-beam projections with equiangular sampling
  void backprj_fanA(double *Radon, int RAY, double **IMAGE, int SIZE, double pos, double **pImage);

private: // attributes
  int     angleNum;         
  double  angle_step;
  double  *AA, *BB;         
  void    Init(int Views);  
};

#endif
