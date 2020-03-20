// ProjectionSimulator.h: interface for the ProjectionSimulator class.
//
/**
 * Copyright (c) 2003 - 2005  Cubic Imaging LLC
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
 * Shepp-Logan phantom has been widely used for evaluating image reconstruction algorithms.
 * This class provides functions for generating the Shepp-Logan phantom and various formats of 
 * projection data. The current version supports parallel- and fan-beam data acquisition geometries.
 * The single circular orbit cone-beam geometry and the spiral-fan-beam geometries is under development
 * and will be added soon.
 *
 * A lot of utility functions are included for some basic image processing.
 *
 * More details reagrding the mathematical theory can be found in J. You's review paper
 * "Image reconstructin from projections in tomography - math and implementation in C++".
 **/

#ifndef _CUBICIMAGING_PROJECTION_SIMULATOR_H_
#define _CUBICIMAGING_PROJECTION_SIMULATOR_H_

//! class for synthesizing the Shepp-Logan phantom and projection data
class ProjectionSimulator  
{
public: // constructor/destructor
  //! default constructor
	ProjectionSimulator();

  //! default destructor
	~ProjectionSimulator();

public: // object/attenuation phantom
  //! synthesizes a discrete 2D Shepp-Logan phantom
  /*!
    \param slImg  - square matrix to hold the Shepp-Logan phantom
    \param dim    - dimension of the square matrix

    The function generates a discrete Shepp-Logan phantom as a square matrix. 
    The memory for slImg has to be allocated in advance. 
    I recommend use CreateMatrix(...) to allocate the memory.
  */
  void GenerateSLPhantom(double ** slImg, int dim);

  //! synthesizes a discrete 3D Shepp-Logan phantom
  /*!
    \param slVol  - volumetric data constructed by a sequence of square matrices
    \param dim    - dimension of the matrix
    \param slice  - number of matrices

    The function generates a 3D volumetric data containing a set of square matrices. To avoid the memory
    management, I recommend use the CreateVolume(...) to allocate the buffer.
  */
  void GenerateSLPhantom(double ***slVol, int dim, int slice);

  //! synthesizes a discrete 2D chest-like attenuation map
  /*!
    \param attMap - square matrix to hold the attneuation map
    \param dim    - dimension of the square matrix

    The function generates a discrete chest-like attenuation map.
  */
  void GenerateAttenMap(double ** attMap, int dim, bool unif = false);

public: // analytical methods to generate projection data
  //! synthesizes the Radon transform of the Shepp-Logan phantom using analytical formula
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of angles for scanning the object in 180
    \param bins   - number of bins used in each scanning view

    The function analytically synthesizes the Radon transform of the Shepp-Logan phantom. 
    The scanning angle range is [0, 180]. No attenuation is included.
  */
  void ProjectorAPar(double **Radon, int views, int bins);

  //! synthesizes the exponential Radon transform of the Shepp-Logan phantom using analytical formula
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of angles for scanning the object
    \param bins   - number of bins used in each scanning view
    \param mu     - the exponential factor

    The function analytically synthesizes the exponential Radon transform of the Shepp-Logan phantom. 
    The scanning angle range is [0, 360].
  */
  void ProjectorAPar(double **Radon, int views, int bins, double mu);

  //! synthesizes fan-beam projection data using analytical formula with equiangular sampling
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of angles for scanning the object
    \param rays   - number of rays used in each scanning view
    \param mu     - the exponential factor

    The function analytically synthesizes the exponential Radon transform of the Shepp-Logan phantom 
    for fan-beam geometry with equiangular sampling. The scanning angle range is [0, 360]. When mu is zero, 
    the generated projections reduce to the standard fan-beam projections.
  */
  void ProjectorAFanA(double **Radon, int views, int rays, double mu = 0.);

  //! synthesizes fan-beam projection data using analytical formula with equally spaced sampling
  /*!
    \param Radon  - 2D array to hold the projections
    \param views  - number of angles for scanning the object
    \param rays   - number of rays used in each scanning view
    \param mu     - the exponential factor

    The function analytically synthesizes the exponential Radon transform of the Shepp-Logan phantom 
    for fan-beam geometry with equally spaced sampling. The scanning angle range is [0, 360].
  */
  void ProjectorAFanS(double **Radon, int views, int rays, double mu = 0.);

public: // numerically generate projections using coordinate transformation
  //! synthesizes the Radon transform using coordinate transformation
  /*!
    \param Image  - square matrix of the discrete object image
    \param dim    - dimension of the square matrix
    \param Radon  - 2D array to hold the projections
    \param views  - number of angles for scanning the object image

    The function numerically synthesizes Radon transform of any object "Image". The scanning angle
    range is [0, 180]. No attenuation is included. As default, the bins in each view is "dim".
  */
  void ProjectorDPar(double **Image, int dim, double **Radon, int views);

  //! synthesizes the exponential Radon transform using coordinate transformation
  /*!
    \param Image  - square matrix of the discrete object image
    \param dim    - dimension of the square matrix
    \param Radon  - 2D array to hold the projections
    \param views  - number of angles for scanning the object
    \param mu     - exponential factor

    The function numerically synthesizes the exponential Radon transform of any object "Image" 
    with exponential factor "mu". The scanning angle range is [0, 360].
  */
  void ProjectorDPar(double **Image, int dim, double **Radon, int view, double mu);

  //! synthesizes the attenuated Radon transform using coordinate transformation
  /*!
    \param Image  - square matrix of the discrete object image
    \param dim    - dimension of the object image
    \param Radon  - 2D array to hold the attenuated Radon transform
    \param views  - number of angles for scanning the object
    \param atten  - square matrix of attenuation map, it has the same dimension with the object image

    The function numerically synthesizes the attenuated Radon transform of any object "Image". 
    The scanning angle range is [0, 360].
  */
  void ProjectorDPar(double **Image, int dim, double **Radon, int views, double **atten);

  //! synthesizes the attenuated fan-beam projection data with equiangular sampling
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of angles for scanning the object
    \param rays   - number of rays used at each fan
    \param Image  - square matrix of the discrete object image
    \param dim    - dimension of the square matrix
    \param atten  - square matrix of attenuation map

    The function numerically synthesizes the attenuated fan-beam projection data of any object "Image" 
    for fan-beam geometry with equiangular sampling. The scanning angle range is [0, 360].
  */
  void ProjectorDFanA(double **Radon, int views, int rays, double **Image, int dim, double **atten = 0);

  //! synthesizes the attenuated fan-beam projection data with equally spaced sampling
  /*!
    \param Radon  - 2D array to hold the projection data
    \param views  - number of angles for scanning the object
    \param rays   - number of rays used at each fan
    \param Image  - square matrix of the discrete object image
    \param dim    - dimension of the square matrix
    \param atten  - square matrix of attenuation map

    The function numerically synthesizes the attenuated fan-beam projection data of any object "Image" 
    for fan-beam geometry with equally spaced sampling. The scanning angle range is [0, 360].
  */
  void ProjectorDFanS(double **Radon, int views, int rays, double **Image, int dim, double **atten = 0);

  //! synthesize the attenuated cone-beam projection data using coordinate transformation for equally spaced 2D sampling
  void ProjectorDConeS(double ***Radon, int angle, int theta, int phi, double ***Image, int dimZ, int dimXY, double ***atten = 0);

public: // utility functions
  //! 3-point smoothing with kernel(0.25, 0.5, 0.25)
  static void Smooth3(double *data, int len, double *temp);

  //! 5-point smoothing with kernel(0.125, 0.5, 0.125)
  static void Smooth5(double **IMAGE, int row, int col);

  //! 9-point smoothing with kernel(0.0625, 0.25, 0.125, 0.0625)
  static void Smooth9(double **IMAGE, int row, int col);

  //! 5-point with 3-order polynomial filter
  static void Smooth53(double *data, int len, double *temp);

  //! 5-point Savitzky-Golay filter (-0.086, 0.343,  0.486)
  static void SavGol5(double *data, int len, double *temp);

  //! 9-point Savitzky-Golay filter (0.035, -0.128,  0.070, 0.315, 0.417)
  static void SavGol9(double *data, int len, double *temp);

  //! 3-point one-dimensional median filter
  static void FilterMedian(double *arr, int len);

  //! 5-point two-dimensional two-dimensional median filter
  static void FilterMedian(double **IMAGE, int row, int col);

  //! coordinate transformation between (x, y) and (s, t)
  static void RotateCarte(double **image_in, double **image_out, int dim, double rotationRadian);

  //! coordinate transformation between (x, y) and (s, t)
  static void RotateCarte(double **image_in, float **image_out, int dim, double rotationRadian);

  //! coordinate transformation from (x, y) to (phi, theta) for equiangular raysum in fan-beam geometry
  static void CarteToPolarA(double **carImage, int dim, double **polImage, int angle, double focalPos);

  //! coordinate transformation from (phi, theta) to (x, y) for equiangular raysum in fan-beam geometry
  static void PolarAToCarte(double **polImage, int angle, double focalPos, double **carImage, int dim);

  //! coordinate transformation from (x, y) to (phi, u) for equally spaced raysum in fan-beam geometry
  static void CarteToPolarS(double **carImage, int dim, double **polImage, int angle, double focalPos);

  //! coordinate transformation from (phi, u) to (x, y) for equally spaced raysum in fan-beam geometry
  static void PolarSToCarte(double **polImage, int angle, double focalPos, double **carImage, int dim);

  //! coordinate transformation from Cartesian to Spheric
  static void CarteToConeP(double ***carImage, int dimZ, int dimXY, double ***coneImage, int theta, int phi, double focalPos);

  //! coordinate transformation from Spheric to Cartesian
  static void ConePToCarte(double ***coneImage, int theta, int phi, double focalPos, double ***carImage, int dimZ, int dimXY);

private: // utility functions for manipulating the Shepp-Logan phantom
  double ellip(double x, double y, double x0, double y0, double phi, double a0, double b0);
  double ellip(double x, double y, double z, double x0, double y0, double z0, double theta, double phi, double a0, double b0, double c0);
  double ellipProj(double x0, double y0, double phi, double a, double b, double l, double theta, double mu);
};

///////////////////////////////////////////////
//! memory allocation utility functions
//

//! allocate memory for matrix data
template<class TT>
inline TT ** CreateMatrix(int row, int col)
{
  TT **matrix;
	matrix = new TT *[row];
  matrix[0] = new TT[row*col];
  memset(matrix[0], 0, sizeof(TT)*row*col);
	for(int i=1; i<row; i++) matrix[i] = matrix[0] + col*i;
  return matrix;
}

//! free the matrix memory
template<class TT>
inline void FreeMatrix(TT ** matrix)
{
  delete[] matrix[0]; delete[] matrix;
}

//! allocate memory for volumetric data
template<class TT>
inline TT *** CreateVolume(int slice, int row, int col)
{
  TT *** volume;
	volume = new TT **[slice];
  volume[0] = new TT*[slice*row];
	for(int i=1; i<slice; i++) volume[i] = volume[0] + row*i;
  volume[0][0] = new TT[slice*row*col];
  memset(volume[0][0], 0, sizeof(TT)*slice*row*col);
  for(int j=0; j<slice; j++) for(int i=0; i<row; i++) volume[j][i] = volume[0][0] + j*row*col + i*col;
  return volume;
}

//! free the volumetric data
template<class TT>
inline void FreeVolume(TT *** volume)
{
  delete[] volume[0][0]; delete[] volume[0]; delete[] volume;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// geometric parameters for parallel-, fan- and cone-beam scanning geometries
// the graphic description can be found in the support documentation
//

//! distance between the rotation center and the data detector plane
const double RADIUS_DETC = 1.0;

//! distance between the rotation center and the fan focal point
const double RADIUS_FOCA = 2.0; //for normal simulation

//! distance between the rotation center and the boundary of object image to be reconstructed
const double RADIUS_IMAGE = 1.0;

//! fan subtending angle covered by the detector, (60 degree from -30 to 30 in this implementation)
const double RAY_VIEW = (3.1415926535897932385/3.0);

//! the diameter of the fan field of view on the detector plane for eqaully-spaced sampling
const double RAY_FOV = (2.0*1.7320508075689);
///////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
