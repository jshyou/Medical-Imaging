// ReconAnalytical.h: interface for the ReconAnalytical class.
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
 * The Filtered Backprojection (FBP) algorithm is the most important method for 
 * image reconstruction and has been used in commercial products for years. 
 * It is available to three typical data acquisition geometries: parallel-, fan- and cone-beam.
 * The ReconAnalytical class collects commonly used flavors of FBP algorithms.  
 * The FBP algorithm for attenuated fan-beam projections and partial scan
 * was developed by J. You recently.
 *
 * The noise strategy developed by You is applied whenever possible. The key ideas if the 
 * combination of median filter and Savitzky-Golay filter. It's simple but quite effective.
 *
 * Circular Harmonic Decomposition (CHD) reconstruction method was developed by J. You 
 * to deal with the rotating geometries such as variable-focal-fan-beam geometry. 
 * That method is originated from Cormack's Nobel-honored work.
 * CHD is able to handle the uniform attenuation in an analytical way without introducing 
 * extra intermediate processing error such as interpolation or nonuniformity. 
 * This algorithm was recommended by Natterer to handle the standard fan-beam data to 
 * avoid the reconstruction nonuniformity in different reconstruction region. The new angular 
 * Hilbert transform based FBP algorithm by You seems to be the option without any restriction 
 * on the sampling grid.
 *
 * The cone-beam image reconstruction is currently under development.
 **/

#ifndef _CUBICIMAGING_IMAGERECON_ANALYTICAL_H_
#define _CUBICIMAGING_IMAGERECON_ANALYTICAL_H_

//! class for various flavors of FBP algorithms and CHD reconstruction algorithms
class ReconAnalytical
{
public:
	ReconAnalytical();
	~ReconAnalytical();

public:
  //! classical parallel-beam FBP algorithm for 180-degree projections
  /*!
    \param Radon  - 2D array to hold projections
    \param views  - number of views for scanning the object image
    \param bins   - number of scanning bins at each view
    \param image  - square matrix for the reconstructed image
    \param dim    - dimension of the reconstructed image
    \param filter - 0 for no noise treatment, and 1 for Savitzky filter only and 2 for median + Savitzky-Golay

    The function implements the classical parallel-beam FBP algorithm in Shepp-Logan's paper. 
    The Poisson noise treatment strategy in You's paper is included.
  */
  void FBP_ParSL( double **Radon,
                  int views, 
                  int bins, 
                  double **image,
                  int dim,
                  int filter = 1
                  );

  //! FBP algorithm for exponential Radon transform with 360-degree projections
  /*!
    \param expRadon - 2D array to hold the exponentially attenuated projections
    \param views    - number of views for scanning the object image
    \param bins     - number of scanning bins at each view
    \param image    - square matrix of the reconstructed image
    \param dim      - dimension of the reconstructed image
    \param mu       - exponential factor in the exponential Radon transform
    \param method   - 1 for Tretiak-Metz formula, and 2 for Novikov formula
    \param filter   - 0 for no noise treatment, and 1 for Savitzky filter only and 2 for median + Savitzky-Golay

    The function implements the Tretiak-Metz and Novikov inversion formulae 
    in the fashion of FBP procedure. With the advent of Novikov's invresion formula, 
    the inverse exponential Radon transform becomes less important as used to be. 
    Also, the Tretiak-Metz formula is noticed to be less stable than the Novikov formula due to 
    the nature of exterior reconstruction.
  */
  void FBP_ParER( double **expRadon,
                  int views,
                  int bins,
                  double **image,
                  int dim,
                  double mu,
                  int method,
                  int filter = 1
                  );

  //! FBP algorithm for attenuated parallel-beam 360-degree projections using Novikov's formula
  /*!
    \param Radon  - 2D array to hold the attenuated projections
    \param views  - number of views for scanning the object image
    \param bins   - number of scanning bins at each view
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param atten  - square matrix of attenuation map
    \param filter - 0 for no noise treatment, and 1 for Savitzky filter only and 2 for median + Savitzky-Golay

    The function implements Novikov's inversion formula using a two-term version of the FBP algorithm 
    based on You's recent work. The Poisson noise treatment strategy is included.
  */
  void FBP_ParYJ( double **attRadon,  // 2D arrar to hold the attenuated projections
                  int views,          // viewing angle numbers in 360 degree
                  int bins,           // number of projection bins in each view
                  double **image,     // square matrix of the reconstructed image
                  int dim,            // the dimension of the reconstructed image
                  double **atten,     // the attenuation map, has to have the same size with projection bins
                  int filter = 0      // filter to smooth the projections
                  );
  void FBP_ParYJ1(double **attRadon,  // 2D arrar to hold the attenuated projections
                  int views,          // viewing angle numbers in 360 degree
                  int bins,           // number of projection bins in each view
                  double **image,     // square matrix of the reconstructed image
                  int dim,            // the dimension of the reconstructed image
                  double **atten,     // the attenuation map, has to have the same size with projection bins
                  int filter = 0      // filter to smooth the projections
                  );

  //! FBP algorithm for attenuated parallel-beam 360-degree projections using Novikov's formula
  /*!
    \param Radon  - 2D arrar to hold the attenuated projections
    \param views  - number of views for scanning the object image
    \param bins   - number of scanning bins at each view
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param atten  - square matrix of attenuation map
    \param filter - 2 for median filter, Savitzky-Golay filter is doomed as default

    The function implements the FBP algorithm in Kunyansky's paper. Kunyansky's implementation seemed 
    to generate smooth reconstruction compared with You's implementation.
    The noise treatment strategy in You's paper is available.
  */
  void FBP_ParKL( double **attRadon,  // attenuated projections
                  int views,          // viewing angle numbers in limited angle
                  int bins,           // projection bins in each view
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double **atten,     // the attenuation map, has to have the same size with projection bins
                  int filter = 1      // filter to smooth the projections
                  );

  //! iterated FBP algorithm for attenuated parallel-beam contiguous 180-degree projections
  /*!
    \param Radon  - 2D array for the attenuated projections
    \param views  - number of views for scanning the object image in 180-degree
    \param bins   - number of scanning bins at each view
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param atten  - square matrix of attenuation map, has to have the same size with projection bins
    \param filter - 2 for median filter, Savitzky-Golay filter is doomed as default
    \param iter   - number of iterations to perform FBP

    The function iterates the Kunyansky's implementation for 180-degree attenuated projection data. 
    The algorithm was proven to be exponentially convergent recently by You.
  */
  void FBP_Par180(double **attRadon,  // attenuated projections
                  int views,          // viewing angle numbers in 360 degree
                  int bins,           // projection bins in each view
                  double **image,     // image to be reconstructed
                  double **atten,     // the attenuation map, has to have the same size with projection bins
                  int filter = 1,     // filter to smooth the projections
                  int iter = 3        // number of iterations
                  );

  //! classical FBP algorithm for fan-beam 360-degree projections with equiangular sampling
  /*!
    \param fanRadon - 2D array for fan-beam projections
    \param rays     - number of scanning bins at each fan view
    \param views    - number of views for scanning the object image
    \param image    - square matrix of reconstructed image
    \param dim      - dimension of reconstructed image
    \param filter   - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters

    The function implements the classical FBP algorithms by Herman and his colleagues 
    for fan-beam projections with equiangular sampling.
    The noise treatment strategy in You's paper is included.
  */

  void FBP_FanAng(double **fanRadon,  // fan-beam projections
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  int filter = 0      // filter to smooth the projections
                  );

  //! Hilbert transfomr based FBP algorithm for (partial) fan-beam projections with equiangular sampling
  /*!
    \param fanRadon - 2D array for fan-beam projections
    \param rays     - number of scanning bins at each fan view
    \param views    - number of views for scanning the object image
    \param image    - square matrix of reconstructed image
    \param dim      - dimension of reconstructed image
    \param filter   - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters

    The function implements the Hilbert transofrm based FBP algorithms with considering the partial fan-beam 
    projections with equiangular sampling. Due to Novikov and Katsevich's work, the Hilbert
    transform has become quite popular in image reconstruction recently. 
    The algorithm implemented here is originally derived in Noo et el's paper.
  */
  void FBP_FanHTT(double **fanRadon,  // fan-beam projections
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  int partial = 0
                  );

  //! Angular Hilbert transfomr based FBP algorithm for (partial) fan-beam projections with equiangular sampling
  /*!
    \param fanRadon - 2D array for fan-beam projections
    \param rays     - number of scanning bins at each fan view
    \param views    - number of views for scanning the object image
    \param image    - square matrix of reconstructed image
    \param dim      - dimension of reconstructed image
    \param filter   - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters

    The function implements the angular Hilbert transofrm based FBP algorithms with considering the partial fan-beam 
    projections with equiangular sampling. The implementation is done under parallel coordinate instead of 
    the conventional fan-beam coordinate. The algorithm would resolve the nonuniformity issue with the fan-beam 
    coordinate and was developed in You's recent work.
  */
  void FBP_FanAHT(double **fanRadon,  // fan-beam projections
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  int partial = 0
                  );

  //! classical FBP algorithm for fan-beam 360-degree projections with equally spaced sampling
  /*!
    \param Radon  - 2-dimensional array for projections
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param filter - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters
    
     The function implements the standard FBP algorithms for fan-beam projections with equally spaced sampling.
     The noise treatment strategy in You's paper is included.
  */
  void FBP_FanSpc(double **fanRadon,  // fan-beam projections
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  int filter = 0      // filter to smooth the projections
                  );

  //! FBP algorithm for attenuated fan-beam 360-degree projections with equiangular sampling
  /*!
    \param Radon  - 2-dimensional array for attenuated projections
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param atten  - square matrix of attenuation map
    \param attDim - dimension of the square matrix of attenuation map
    \param fanAtt - fan-beam projections of attenuation map
    \param filter - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters

    The function implements the generalized FBP algorithms for attenuated fan-beam projections with equiangular sampling. 
    The algorithm was developed by J. You recently.
  */
  void FBP_FanAng(double **fanRadon,  // attenuated fan-beam projections
                  int views,          // viewing angle numbers in 360 degree
                  int rays,           // projection rays at each fan
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double **atten,     // the attenuation map
                  int attDim,         // the size of attenuation map
                  double ** fanAtt,   // fan-beam-projected attenuation
                  int filter = 0      // filter to smooth the projections
                  );

  //! FBP algorithm for attenuated fan-beam 360-degree projections with equally spaced sampling
  /*!
    \param Radon  - 2-dimensional array for attenuated projections
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix of reconstructed image
    \param dim    - dimension of reconstructed image
    \param atten  - square matrix of attenuation map
    \param attDim - dimension of the square matrix of attenuation map
    \param fanAtt - fan-beam projections of attenuation map
    \param filter - 0 for no noise treatment, and 1 for noise treatment by median + Savitzky-Golay filters

    The function implements the generalized FBP algorithms for attenuated fan-beam projections with equally spaced sampling. 
    The algorithm was developed by J. You recently.
  */
  void FBP_FanSpc(double **fanRadon,  // attenuated fan-beam projections
                  int views,          // viewing angle numbers in 360 degree
                  int rays,           // projection rays at each fan
                  double **image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double **atten,     // the attenuation map
                  int attDim,         // the size of attenuation map
                  double ** fanAtt,   // fan-beam-projected attenuation
                  int filter = 0      // filter to smooth the projections
                  );

  //! DBF algorithm for fan-beam 360-degree projection data with equiangular sampling and uniform attenuation
  /*!
    \param Radon  - 2-dimensional array for exponential Radon transform
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix for the reconstructed image
    \param dim    - dimension of reconstructed image
    \param mu     - exponential factor

    The function implements the CHD algorithm for fan-beam 360-degree projections with equiangular sampling. 
    The projection should be the exponential Radon transform with fan-beam coordinate.
  */
  void DBF_FanAng(double **fanRadon,  // fan-beam exponential Radon transform
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **Image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double mu           // the exponential factor
                  );

  //! CHD algorithm for fan-beam 360-degree projection data with equiangular sampling and uniform attenuation
  /*!
    \param Radon  - 2-dimensional array for exponential Radon transform
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix for the reconstructed image
    \param dim    - dimension of reconstructed image
    \param mu     - exponential factor

    The function implements the CHD algorithm for fan-beam 360-degree projections with equiangular sampling. 
    The projection should be the exponential Radon transform with fan-beam coordinate.
  */
  void CHD_FanAng(double **fanRadon,  // fan-beam exponential Radon transform
                  int rays,           // projection rays at each fan
                  int views,          // viewing angle numbers in 360 degree
                  double **Image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double mu           // the exponential factor
                  );

  //! CHD algorithm for fan-beam 360-degree  projection data with equally spaced sampling and uniform attenuation
  /*!
    \param Radon  - 2-dimensional array for exponential Radon transform
    \param rays   - number of scanning bins at each fan view
    \param views  - number of views for scanning the object image
    \param image  - square matrix for the reconstructed image
    \param dim    - dimension of the reconstructed image
    \param mu     - exponential factor

    The function implements the CHD algorithm for fan-beam 360-degree projections with equally spaced sampling.
    The projection should be the exponential Radon transform with fan-beam coordinate.
  */
  void CHD_FanSpc(double **fanRadon,  // fan-beam exponential Radon transform
                  int rays,           // projection bins in each view
                  int views,          // viewing angle numbers in 360 degree
                  double **Image,     // image to be reconstructed
                  int dim,            // the dimension of reconstructed image
                  double mu           // the exponential factor
                  );

  //! attenuated Hilbert transform kernel
  double HilbertKer(double s, double step);
  double HilbertHyperCos(double s, double step, double mu);
  double HyperCosDiff(double s, double step, double mu);

  //! inverse attenuated Hilbert transform kernel
  double HilbertHyperCosInv(double s, double step, double mu, double low, double up, double s1);

private: 
  //! discrete version of the standard Hilbert transform
  double HilbertSTD(int n);

  //! Kunyansky's discrete Hilbert transform
  double HilbertKL(int n);

  //! discrete version of the Hilbert transform for inverse exponential Radon transform
  double HilbertAttCos(double s, double step, double mu);

  //! rectangle-weighted Hilbert transform
  double HilbertRect(double s, double step);

  //! cosine-weighted Hilbert transform
  double HilbertCos(double s, double d);

  //! discrete version of the angular Hilbert transform
  double HilbertAng(int n, double step);

  //! discrete version of rectangle-weighted ramp filter
  double RampRect(int n);

  //! rectangle-weighted ramp filter
  double RampRect(double s, double step);

  //! discrete version of cosine-weighted ramp filter
  double RampCos(int n);

  //! cosine-weighted ramp filter
  double RampCos(double s, double step);

  //! discrete version of Shepp-Logan filter
  double RampSL(int n);

  //! Shepp-Logan ramp filter
  double RampSL(double s, double step);

  //! Shepp-Logan ramp filetr including the uniform attenuator
  double RampSL(double s, double step, double attenuator);

  //! sinc function
  double Sinc(double s);

  //! fast Fourier transform
  void FFT(double data[], int nn, int isign);
};

#endif 
