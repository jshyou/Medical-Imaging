// dibimage.h

/**
 * Copyright (c) 2001 - 2005  Cubic Imaging LLC
 *
 * Permission to copy, modify and distribute the software and its documentation 
 * for noncommercial use is hereby granted without fee, provided that the above 
 * copyright notice appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation. Cubic Imaging LLC 
 * makes no representations about the suitability of this software for any purpose. 
 * It is provided "as is" without express or implied warranty.
 * 
 * Permission to copy, modify, and sell the software and its documentation 
 * in commercial use for profits has to be approved by written notice from 
 * Cubic Imaging LLC
 * 
 * \author Jason You (jyou@cubic-imaging.com)
 * \date created by 08/08/2001
 *
 * This file defines a DIBImage class for generic image objects and some utility 
 * structures. The standard memory management and some basic processing methods 
 * are included. Currently, the class fully supports 8-, 16-bit grayscale, and 24, 32-bit 
 * color images. The data buffer is allocated in a contiguous chunk and arranged in 
 * the format of Windows Screen Coordinate convention.
 *
 * The primary attributes contain: row, column, number of pixel bits and bits buffer. 
 *
 */

#ifndef _CUBICIMAGING_DEVICE_INDEPENDENT_BIT_IMAGE_H_
#define _CUBICIMAGING_DEVICE_INDEPENDENT_BIT_IMAGE_H_

//! Class for device-independent bitmap image objects
class DIBImage
{
public:
	//! -- default constructor
	DIBImage();

  //! -- copy constructor                         
	DIBImage(const DIBImage & dibImg);

  //! -- assignment operator
	DIBImage & operator= (const DIBImage & dibImg); 

  //! -- default destructor  
	~DIBImage();                        
  
public: // status of image object
  //! -- round length must be 4 in Windows environment
	enum {round = 4};

  //! -- check if the image object is empty
  bool IsEmpty();

  //! -- clear the image object
  void ClearImage();

public: // accessors to native image attributes
	//! -- return the height of image
	short row() const { return rows; } 

  //! -- return the width of image
	short col() const { return cols; }
  
  //! -- return the number of bits used for each pixel
	char pxlbit() const { return bitsPerPxl; }

  //! -- return the ownership of image bits buffer
  bool getOwnership() const { return ownership; }  

  //! -- return the address to image bits buffer
	unsigned char* buff() const { return buffer; }

public: // accessors to calculated image attributes
  //! -- return the length of round-padded row in bytes
	int GetRowLen() const;
  
  //! -- return the length of total allocated buffer in bytes
  int GetBuffLen() const; 

  //! -- return the buffer address of the row-th line, NULL if not available
  /*!
    \param row - image line number

    This function returns the buffer address of the row-th line of image data, 
    returns NULL if not available. The pointer needs to be casted to the
    correspnding type for 16-, 24-bit image.
  */
	unsigned char* getRow(int row) const;

  //! -- return the pointer to the pixel at location (row, col), NULL if not available
  /*!
    \param row - image row number
    \param col - image column number

    This function returns the pointer to the pixel at the location of(row, col)
    under the Windows Screen coordinate, returns NULL when not available. 
    The pointer needs to be casted to the correspnding type for 16-, 24-bit image.
  */
	unsigned char* getPxl(int row, int col) const;	

public: // manage the image size/data
	//! -- creates/resizes an image object according to predefined size
  /*!
    \param row - image height
    \param col - image width
    \param bits - number of bits for each pixel

    This function creates a new image object or resizes an old image object
    by the new dimension (row, col). The pixel values are initialized to 0.
    Return false if the image object creation failed.
  */
	bool CreateImage(int row, int col, int bits = 8);

  //! -- creates/resizes an image object by copying other image data
  /*!
    \param data - image bit buffer
    \param row - image height
    \param col - image width
    \param rowLen - length of the allocated buffer for each row
    \param bits - number of bits for each pixel

    This function creates a new image object or resizes an old image object
    through copying other image data or any matrix of bits. Return false if the copy failed.
  */
	bool CopyImage(const char* data, int row, int col, int rowLen, int bits = 8);

  //! -- initialize image object through attaching other image buffer;
  /*!
    \param data - image bit buffer
    \param row - image height
    \param col - image width
    \param bits - number of bits for each pixel

    This function creates a new image object or replaces an old image object through 
    attaching external image bits buffer. The data buffer "data" has to follow the Windows
    BMP data structure.
  */
  void AttachImage(const char* data, int row, int col, int bits = 8);
  
  //! -- modify part of the image 
  /*!
    \param row - starting line number
    \param col - starting column number
    \param subImg - the source image object

    This function replaces the region of interest starting from the location (row, col) 
    in the overlapping region of this image object and the external image object with subImg.
  */
	void SetSubImage(int row, int col, const DIBImage & subImg);

  //! -- copy part of image to other image object
  /*!
    \param row - starting row number
    \param col - starting column number
    \param subImg - the destination image object

    This function copies the region of interest starting from the location (row, col) 
    in the overlapping region of this image object and the external image object to the 
    external image object subImg.
  */
	void GetSubImage(int row, int col, DIBImage & subImg);

  //! -- copy image to a contiguous buffer "da", da must be allocated outside, any padded bits will be removed!!!
  /*!
    \param data - the destination bits buffer

    This function copies the current image data to the outside contiguous buffer.
    The padded bytes will be removed during the copy.
  */
  void CopyImageTo(void* data); 
                              
  //! -- read bmp image file, 8- and 24-bits only
  bool ReadBMPFile(const char* name);

  //! -- write bmp image file, 8- and 24-bits only
  bool WriteBMPFile(const char* name);

  //! -- read the raw pixel data, the raw pixel data structure has to follow the Windows BMP format
  bool ReadRawPixelData(int row, int col, int pxlBits, const char* fileName, int filePosShift = 0);

public:  // utility statistic functions  
  //! -- get the min/max value for 8-bit image
  void GetMinMax(unsigned char & min, unsigned char & max);

  //! -- get the min/max value for 16-bit image
  void GetMinMax(unsigned short & min, unsigned short & max);

  //! -- get the min/max value for 24-bit image
  void GetMinMax(unsigned char * min, unsigned char * max);

  //! -- get median
  int  GetMedian( );

  //! -- get the mean/deviation for 8-, 16-bit grayscale image
  void GetStat(double & mean, double & deviation);

  //! -- get the mean/deviation for 24-bit color image
  void GetStat(double mean[], double deviation[]);

  //! -- get the histogram buffer pointer, user must be responsible for the memory release
  int* GetHistogram();

  //! -- get the histogram buffer pointer for ROI, including the boundary. User must be responsible for the memory release
  int* GetHistogram(int top, int left, int bottom, int right);

public:  // utility color/coodinate processing functions
  //! -- reverse the image color
  void InvertColor();                     

  //! -- quantize the image into 1 or 2-bit image, only available for 8-bit image
  void QuantizeColor(int level = 2);
  
  //! -- horizontal flip
  void FlipHorizontal();
  
  //! -- vertical flip
  void FlipVertical();    
  
  //! -- detect the bounding box, only available for 8-bit image
  void GetBoundingBox(int & top, int & left, int & bottom, int & right);

  //! -- shift image to different position                
  void Shift(int rowShift, int colShift); 
  
  //! -- rotate image around the center              
  void Rotate(const double degree, bool enlargeSize = false); 

protected: // native image attributes
  //! -- image height
  short rows;

  //! -- image width
  short cols;

  //! -- the number of bits used for representing each pixel
  char  bitsPerPxl;

  //! -- ownership of the image buffer
  bool ownership;

	//! -- address to image bits buffer
  unsigned char* buffer;

private:
  // -- helper functions for adaptive thresholding
  void quantize2();
  void quantize3();
  void quantize4();
};

// the following data structures are for the color modeling

// color modeling
struct RGBColorStr;
struct RGBAlphaStr;
struct YIQColorStr;
struct HSIColorStr;

//! RGB color model strucure, following the Windows convention
struct RGBColorStr
{
  //! -- the blue component
  unsigned char B;
  //! -- the green component
  unsigned char G;
  //! -- the red component
  unsigned char R;

  //! -- default constructor
  RGBColorStr(unsigned char b=0, unsigned char g=0, unsigned char r=0) : B(b), G(g), R(r) { }

  //! -- convert RGB color model to YIQ color model
  void ToYIQ( YIQColorStr* yiq );

  //! -- convert RGB color model to HSI color model
  void ToHSI( HSIColorStr* hsi );
};

//! RGB color model strucure, following the Windows convention
struct RGBAlphaStr
{
  //! -- the blue component
  unsigned char B;
  //! -- the green component
  unsigned char G;
  //! -- the red component
  unsigned char R;
  //! -- the alpha value
  unsigned char A;

  //! -- default constructor
  RGBAlphaStr(unsigned char b=0, unsigned char g=0, unsigned char r=0, unsigned char a=0) 
    : B(b), G(g), R(r), A(a) { }
};

//! Luminance-Chrominance color model strucure
struct YIQColorStr
{
  //! -- the luminance
  double Y;
  //! -- the chrominiance
  double I;
  //! -- the chrominiance
  double Q;

  //! -- default constructor
  YIQColorStr(double y = 0., double i = 0., double q = 0.) : Y(y), I(i), Q(q) { }

  //! -- convert YIQ color model to RGB color model
  void ToRGB(RGBColorStr* rgb);
};

//! Hue-Saturation color model strucure
struct HSIColorStr
{
  //! -- the hue
  double H;
  //! -- the saturation
  double S;
  //! -- the intensity
  double I;

  //! -- default constructor
  HSIColorStr(double h = 0., double s = 0., double i = 0.) : H(h), S(s), I(i) { }

  //! -- convert HSI color model to RGB color model
  void ToRGB(RGBColorStr* rgb);
};

//! BMP file header structure
struct BMPFileHeader 
{
  short bfType; 
  short bfSize[2]; 
  short bfReserved1; 
  short bfReserved2; 
  short bfOffBits[2];
};

//! BMP info header structure
struct BMPInfoHeader
{
  int   biSize; 
  int   biWidth; 
  int   biHeight; 
  short biPlanes; 
  short biBitCount; 
  int   biCompression; 
  int   biSizeImage; 
  int   biXPelsPerMeter; 
  int   biYPelsPerMeter; 
  int   biClrUsed; 
  int   biClrImportant; 
}; 

#endif //_DEVICE_INDEPENDENT_BIT_IMAGE_H_
