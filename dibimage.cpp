// dib_image.cpp
 
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include "dibimage.h"
#include "utilFunc.h"

DIBImage::DIBImage()
  : rows(0), cols(0), bitsPerPxl(8), buffer(NULL), ownership(true)
{
}

DIBImage::DIBImage(const DIBImage & dibImg) 
  : rows(0), cols(0), bitsPerPxl(8), buffer(NULL), ownership(true)
{ 
	*this = dibImg; 
}

DIBImage::~DIBImage()
{
  if(ownership) free(buffer);
}

DIBImage & DIBImage::operator = (const DIBImage & dibImg)
{
  this->CreateImage(dibImg.row(), dibImg.col(), dibImg.pxlbit());
	int len = this->GetBuffLen();
	memcpy(buffer, dibImg.buff(), len);

	return *this;
}

bool DIBImage::IsEmpty()
{
  return !( rows>0 && cols>0 && buffer!=NULL && (bitsPerPxl%8==0) );
}

void DIBImage::ClearImage()
{
  if(ownership) free(buffer);

  buffer = NULL;
  rows = 0;
  cols = 0;
  bitsPerPxl = 8;
  ownership = true;
}

unsigned char* DIBImage::getRow(int r) const
{
	if( r<0 || r>=rows || buffer==NULL ) return NULL;
	return (buffer+r*GetRowLen());
}

unsigned char* DIBImage::getPxl(int r, int c) const
{
	if( r<0 || r>=rows || c<0 || c>=cols || buffer==NULL ) return NULL;
	return (buffer+r*GetRowLen()+c*bitsPerPxl/8);
}

int DIBImage::GetRowLen() const
{ 
  return (((cols*bitsPerPxl/8+round-1)/round)*round); 
}

int DIBImage::GetBuffLen() const
{ 
  return (row()*GetRowLen()); 
}  

bool DIBImage::CreateImage(int r, int c, int bits)
{
  if( r == 0 || c == 0 )
  {
    ClearImage();
    return false;
  }

  int rowLen = ((c*(bits/8)+round-1)/round)*round;
	int len = r*rowLen;

  if( !ownership ) // check if the current image is attached
  {
    buffer = (unsigned char*)calloc(1, len);
  }
	else if( GetBuffLen() != len )
	{
		if(buffer == NULL) buffer = (unsigned char*)calloc(1, len);
    else buffer = (unsigned char*)realloc(buffer, len);
	}

  if( buffer == NULL )
  {
    ClearImage();
    return false;
  }

	rows = r;
	cols = c;
	bitsPerPxl = bits;
  ownership = true;

  return true;
}

bool DIBImage::CopyImage(const char* da, int r, int c, int rowLen, int bits)
{
	if( !CreateImage(r, c, bits) ) return false;

	for(int i=0; i<r; i++) memcpy(getRow(i), da+i*rowLen, rowLen);

  return true;
}

void DIBImage::AttachImage(const char* da, int r, int c, int bits)
{
  ClearImage(); // must clear any existing image data

	rows = r;
	cols = c;
	bitsPerPxl = bits;
  buffer = (unsigned char*)da;
  ownership = false;
}

void DIBImage::CopyImageTo( void* da)
{
	int row_len = cols*bitsPerPxl/8;
	for(int i=0; i<rows; i++) memcpy((char*)da+i*row_len, getRow(i), row_len);
}

void DIBImage::SetSubImage(int r, int c, const DIBImage & subImg)
{
	if( r<0 || c<0 || r>=rows || c>=cols || buffer == NULL ) return;
	if(bitsPerPxl != subImg.pxlbit()) return;

	int row = subImg.row();
	int col = subImg.col();
	int rw = (rows-r)>row? row : (rows-r);
	int cl = (cols-c)>col? col : (cols-c);

  int bytes = bitsPerPxl/8;

	unsigned char* in;
	unsigned char* out;
	for(int i=0; i<rw; i++)
	{
		in = subImg.getRow(i);
		out = getRow(i+r);
		memcpy(out+c*bytes, in, cl*bytes);
	}
}

void DIBImage::GetSubImage(int r, int c, DIBImage & subImg)
{
	if( r<0 || c<0 || r>=rows || c>=cols || buffer == NULL ) return;

	if(bitsPerPxl != subImg.pxlbit()) return;

	int row = subImg.row();
	int col = subImg.col();
	int rw = (rows-r)>row? row : (rows-r);
	int cl = (cols-c)>col? col : (cols-c);

  int bytes = bitsPerPxl/8;

	unsigned char* in;
	unsigned char* out;
	for(int i=0; i<rw; i++)
	{
		in = getRow(i+r);
		out = subImg.getRow(i);
		memcpy(out, in+c*bytes, cl*bytes);
	}
}

bool DIBImage::ReadBMPFile(const char* name)
{
  FILE* fp;
  fp = fopen(name, "rb");
  if( fp == NULL ) return false;
  size_t rowLen, len = 0;
  int i;

  BMPFileHeader bmfh;
  BMPInfoHeader bmih;
  len = sizeof(BMPFileHeader);
  if(fread(&bmfh, 1, len, fp) != len) goto ErrorRet;
  if(bmfh.bfType != (('M' << 8) | 'B')) goto ErrorRet;
  len = sizeof(BMPInfoHeader);
  if(fread(&bmih, 1, len, fp) != len) goto ErrorRet;
  if(bmih.biBitCount != 8 && bmih.biBitCount != 24) goto ErrorRet;
  CreateImage(abs(bmih.biHeight), bmih.biWidth, bmih.biBitCount);
  if(buff() == NULL) goto ErrorRet;

  if(bitsPerPxl == 8)
  {
    if(fseek(fp, 1024, SEEK_CUR)) goto ErrorRet;
  }

  rowLen = GetRowLen();
  for(i=0; i<rows; i++)
  {
    if(fread(getRow(rows-1-i), 1, rowLen, fp) != rowLen) goto ErrorRet;
  }

  fclose(fp);
  return true;
ErrorRet:
  fclose(fp);
  return false;
}

bool DIBImage::WriteBMPFile(const char* name)
{
  FILE* fp;
  fp = fopen(name, "wb");
  if( fp == NULL ) return false;

  int i, dwBufferSize = GetBuffLen();
  size_t len, rowLen;
  RGBAlphaStr *cp = NULL;

// write file header
  BMPFileHeader bmfh;
	bmfh.bfType = ('M' << 8) | 'B';
	*((int*)bmfh.bfSize) = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader) + dwBufferSize;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	*((int*)bmfh.bfOffBits) = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader);
  if( bitsPerPxl == 8 )
  {
		*((int*)bmfh.bfSize) = *((int*)bmfh.bfSize) + 256 * sizeof(RGBAlphaStr);
		*((int*)bmfh.bfOffBits) = *((int*)bmfh.bfOffBits) + 256 * sizeof(RGBAlphaStr);
    cp = new RGBAlphaStr[256];
    for(i=0; i<256; i++)
    {
      cp[i].B = cp[i].G = cp[i].R = i;
    }
  }
  len = sizeof(BMPFileHeader);
  if(fwrite(&bmfh, 1, len, fp) != len) goto ErrorRet;

// write bmp header
  BMPInfoHeader bmih;
	bmih.biSize = sizeof(BMPInfoHeader);
	bmih.biWidth = col();
	bmih.biHeight = row();
	bmih.biPlanes = 1;
	bmih.biBitCount = bitsPerPxl;
	bmih.biCompression = 0;
	bmih.biSizeImage = 0;
	bmih.biXPelsPerMeter = 0;
	bmih.biYPelsPerMeter = 0;
	bmih.biClrUsed =0;
	bmih.biClrImportant = 0;
  len = sizeof(BMPInfoHeader);
  if(fwrite(&bmih, 1, len, fp) != len) goto ErrorRet;

// write color palette
  if( bitsPerPxl == 8 )
  {
    if(fwrite(cp, 1, 1024, fp) != 1024) 
    {
      delete[] cp;
      goto ErrorRet;
    }
    delete[] cp;
  }

// write pixle data
  rowLen = GetRowLen();
  for(i=0; i<rows; i++)
  {
    if(fwrite(getRow(rows-1-i), 1, rowLen, fp) != rowLen) goto ErrorRet;
  }

  fclose(fp);
  return true;
ErrorRet:
  fclose(fp);
  return false;
}

bool DIBImage::ReadRawPixelData(int row, int col, int pxlBits, const char* fileName, int filePosShift)
{
  FILE* fp;
  fp = fopen(fileName, "rb");
  if( fp == NULL ) return false;

  fseek(fp, filePosShift, SEEK_SET);
  DIBImage dib;
  dib.CreateImage(row, col, pxlBits);

  int rowLen = dib.GetRowLen();
  for(int i=0; i<row; i++)
  {
    fread(dib.getRow(row-1-i), rowLen, 1, fp);
  }

  fclose(fp);
  return true;
}

int* DIBImage::GetHistogram()
{
  if( buffer == NULL ) return NULL;
  int len;
  if ( bitsPerPxl == 8 ) len = 0xFF+1;
  else if( bitsPerPxl == 16 ) len = 0xFFFF+1;
  else if ( bitsPerPxl == 24 ) len = (0xFF+1)*3;
  else return NULL;

  int* histo = (int*) calloc(len, sizeof(int));
  int i, j;

  if( bitsPerPxl == 8 )
  {
    unsigned char* pt;
    for(i=0; i<rows; i++)
    {
      pt = getRow(i);
      for(j=0; j<cols; j++)
      {
        histo[pt[j]]++;
      }
    }
  }
  else if( bitsPerPxl == 16 )
  {
    unsigned short* pt;
    for(i=0; i<rows; i++)
    {
      pt = (unsigned short*)getRow(i);
      for(j=0; j<cols; j++)
      {
        histo[pt[j]]++;
      }
    }
  }
  else if ( bitsPerPxl == 24 )
  {
    RGBColorStr* pt;
    for(i=0; i<rows; i++)
    {
      pt = (RGBColorStr*)getRow(i);
      for(j=0; j<cols; j++)
      {
        histo[pt[j].B]++;
        histo[pt[j].G+256]++;
        histo[pt[j].R+512]++;
      }
    }
  }

  return histo;
}

int* DIBImage::GetHistogram(int top, int left, int bottom, int right)
{
  if( IsEmpty() ) return NULL;
  if( top < 0 || top >= rows || left < 0 || left >= cols || bottom >= rows || right >= cols ) return NULL;

  int len;
  if ( bitsPerPxl == 8 ) len = 0xFF+1;
  else if( bitsPerPxl == 16 ) len = 0xFFFF+1;
  else if ( bitsPerPxl == 24 ) len = (0xFF+1)*3;
  else return NULL;

  int* histo = (int*) calloc(len, sizeof(int));
  int i, j;

  if( bitsPerPxl == 8 )
  {
    unsigned char* pt;
    for(i=top; i<=bottom; i++)
    {
      pt = getRow(i);
      for(j=left; j<=right; j++)
      {
        histo[pt[j]]++;
      }
    }
  }
  else if( bitsPerPxl == 16 )
  {
    unsigned short* pt;
    for(i=top; i<=bottom; i++)
    {
      pt = (unsigned short*)getRow(i);
      for(j=left; j<=right; j++)
      {
        histo[pt[j]]++;
      }
    }
  }
  else if ( bitsPerPxl == 24 )
  {
    RGBColorStr* pt;
    for(i=top; i<=bottom; i++)
    {
      pt = (RGBColorStr*)getRow(i);
      for(j=left; j<=right; j++)
      {
        histo[pt[j].B]++;
        histo[pt[j].G+256]++;
        histo[pt[j].R+512]++;
      }
    }
  }

  return histo;
}

void DIBImage::GetMinMax(unsigned char & min, unsigned char & max)
{
  if( IsEmpty() || bitsPerPxl != 8 ) return;
  min = max = buffer[0];
  unsigned char* pt;
  int i, j;

  for(i=0; i<rows; i++)
  {
    pt = getRow(i);
    j = 0;
    while( j < cols )
    {
      if( *pt < min ) min = *pt;
      if( *pt > max ) max = *pt;

      j ++;
      pt ++;
    }
  }
}

void DIBImage::GetMinMax(unsigned short & min, unsigned short & max)
{
  if(IsEmpty() || bitsPerPxl != 16) return;
  min = max = buffer[0];
  unsigned short* pt;
  int i, j;

  for(i=0; i<rows; i++)
  {
    pt = (unsigned short*)getRow(i);
    j = 0;
    while( j < cols )
    {
      if( *pt < min ) min = *pt;
      if( *pt > max ) max = *pt;

      j ++;
      pt ++;
    }
  }
}

void DIBImage::GetMinMax(unsigned char * min, unsigned char * max)
{
  if( IsEmpty() || bitsPerPxl != 24 ) return;
  min[0] = max[0] = buffer[0];
  min[1] = max[1] = buffer[1];
  min[2] = max[2] = buffer[2];
  unsigned char* pt;
  int i, j;

  for(i=0; i<rows; i++)
  {
    pt = getRow(i);
    j = 0;
    while( j < cols )
    {
      if( *pt < *min ) *min = *pt;
      if( *pt > *max ) *max = *pt;
      pt ++;

      if( *pt < *(min+1) ) *(min+1) = *pt;
      if( *pt > *(max+1) ) *(max+1) = *pt;
      pt ++;

      if( *pt < *(min+2) ) *(min+2) = *pt;
      if( *pt > *(max+2) ) *(max+2) = *pt;
      pt ++;

      j ++;
    }
  }
}

int DIBImage::GetMedian()
{
  if( IsEmpty() ) return -1;
  if( bitsPerPxl != 8 && bitsPerPxl != 16) return -1;

  int* hist = GetHistogram();
  int half = row()*col()/2;
  int sum = 0;
  int dim = 1 << bitsPerPxl;

  for(int i = 0; i<dim; i++)
  {
    sum += hist[i];

    if( sum >= half )
    {
      half = i;
      break;
    }
  }

  free(hist);

  return half;
}

void DIBImage::GetStat(double & mean, double & deviation)
{
  mean = 0.0;
  deviation = 0.0;

  if( IsEmpty() ) return;
  if( bitsPerPxl != 16 && bitsPerPxl != 8 ) return;

  double mean_ = 0;
  double deviation_ = 0;
  int i, j;

  if( bitsPerPxl == 8 )
  {
    unsigned char* pt;
    for(i=0; i<rows; i++)
    {
      pt = getRow(i);
      j = 0;
      while( j < cols )
      {
        mean_ += *pt;
        deviation_ += (*pt)*(*pt);

        pt ++;
        j ++;
      }
    }
  }
  else if( bitsPerPxl == 16 )
  {
    unsigned short* pt;
    for(i=0; i<rows; i++)
    {
      pt = (unsigned short*)getRow(i);
      j = 0;
      while( j < cols )
      {
        mean_ += *pt;
        deviation_ += (*pt)*(*pt);

        pt ++;
        j ++;
      }
    }
  }

  if( rows*cols == 0 ) return;
  mean = mean_/double(rows*cols);
  deviation = deviation_/double(rows*cols);
  deviation -= mean*mean;
  if( deviation < 0. ) deviation = 0.;
  else deviation = sqrt(deviation);
}

void DIBImage::GetStat(double mean[], double deviation[])
{
  mean[0] = 0.0;
  mean[1] = 0.0;
  mean[2] = 0.0;
  deviation[0] = 0.0;
  deviation[1] = 0.0;
  deviation[2] = 0.0;

  if( IsEmpty() || bitsPerPxl != 24 ) return;
  if( rows*cols == 0 ) return;

  unsigned char* pt;
  int i, j;
  int *mean_ = (int*) calloc(3, sizeof(int));
  int *dev_ = (int*) calloc(3, sizeof(int));

  for(i=0; i<rows; i++)
  {
    pt = getRow(i);
    j = 0; 
    while( j < cols )
    {
      *mean_ += *pt;
      *dev_ += (*pt)*(*pt);
      pt ++;
      mean_ ++;
      dev_ ++;

      *mean_ += *pt;
      *dev_ += (*pt)*(*pt);
      pt ++;
      mean_ ++;
      dev_ ++;

      *mean_ += *pt;
      *dev_ += (*pt)*(*pt);
      pt ++;

      j ++;
      mean_ -= 2;
      dev_ -= 2;
    }
  }

  double sz = rows*cols;
  mean[0] = mean_[0]/sz;
  mean[1] = mean_[1]/sz;
  mean[2] = mean_[2]/sz;
  deviation[0] = dev_[0]/sz;
  deviation[1] = dev_[1]/sz;
  deviation[2] = dev_[2]/sz;
  deviation[0] -= mean[0]*mean[0];
  deviation[1] -= mean[1]*mean[1];
  deviation[2] -= mean[2]*mean[2];

  if( deviation[0] < 0. ) deviation[0] = 0.;
  else deviation[0] = sqrt(deviation[0]);
  if( deviation[1] < 0. ) deviation[1] = 0.;
  else deviation[1] = sqrt(deviation[1]);
  if( deviation[2] < 0. ) deviation[2] = 0.;
  else deviation[2] = sqrt(deviation[2]);

  free(mean_);
  free(dev_);
}

void DIBImage::InvertColor()
{
  if( IsEmpty() ) return;

  int i, j;

  if( bitsPerPxl == 8 )
  {
    unsigned char* pt;
    for(i=0; i<rows; i++)
    {
      pt = getRow(i);
      for(j=0; j<cols; j++)
      {
        pt[j] ^= 0xFF;
      }
    }
  }
  if( bitsPerPxl == 16 )
  {
    unsigned short* pt;
    for(i=0; i<rows; i++)
    {
      pt = (unsigned short*)getRow(i);
      for(j=0; j<cols; j++)
      {
        pt[j] ^= 0xFFFF;
      }
    }
  }
  else if( bitsPerPxl == 24 )
  {
    RGBColorStr* rgb;
    for(i=0; i<rows; i++)
    {
      rgb = (RGBColorStr*)getRow(i);
      for(j=0; j<cols; j++)
      {
        rgb[j].B ^= 0xFF;
        rgb[j].G ^= 0xFF;
        rgb[j].R ^= 0xFF;
      }
    }
  }
}

void DIBImage::QuantizeColor(int level)
{
  if( IsEmpty() || bitsPerPxl != 8 ) return;

  if( level == 2 )
  {
    quantize2();
    return;
  }
  else if( level == 3 )
  {
    quantize3();
    return;
  }
  else if( level == 4 )
  {
    quantize4();
    return;
  }
  else
  {
    return;
  }
}

void DIBImage::FlipHorizontal()
{
  if( IsEmpty() ) return;

  int i, j;

  if( bitsPerPxl == 8 )
  {
    unsigned char* buff;
    unsigned char temp;
    for(i=0; i<rows; i++)
    {
      buff = getRow(i);
      for(j=0; j<cols/2; j++)
      {
        temp = buff[j];
        buff[j] = buff[cols-1-j];
        buff[cols-1-j] = temp;
      }
    }
  }
  else if( bitsPerPxl == 24 )
  {
    RGBColorStr* rgb;
    RGBColorStr temp;
    for(i=0; i<rows; i++)
    {
      rgb = (RGBColorStr*)getRow(i);
      for(j=0; j<cols/2; j++)
      {
        temp.B = rgb[j].B;
        rgb[j].B = rgb[cols-1-j].B;
        rgb[cols-1-j].B = temp.B;

        temp.G = rgb[j].G;
        rgb[j].G = rgb[cols-1-j].G;
        rgb[cols-1-j].G = temp.G;

        temp.R = rgb[j].R;
        rgb[j].R = rgb[cols-1-j].R;
        rgb[cols-1-j].R = temp.R;
      }
    }
  }
}

void DIBImage::FlipVertical()
{
  if( IsEmpty() ) return;

  int rLen = GetRowLen();
  char* temp = (char*) calloc(1, rLen);
  unsigned char *buff1, *buff2;

  for(int i=0; i<rows/2; i++)
  {
    buff1 = getRow(i);
    buff2 = getRow(rows-1-i);

    memcpy(temp, buff1, rLen);
    memcpy(buff1, buff2, rLen);
    memcpy(buff2, temp, rLen);
  }

  free(temp);
}

void DIBImage::GetBoundingBox(int & top, int & left, int & bottom, int & right)
{
  int i, j;
  unsigned char * bf;

  top = 0;
  for(i=0; i<rows; i++)
  {
    bf = getRow(i);
    for(j=0; j<cols; j++)
    {
      if( bf[j] ) { top = i; left = j; i = rows; break; }
    }
  }

  bottom = rows-1;
  for(i=rows-1; i>=0; i--)
  {
    bf = getRow(i);
    for(j=0; j<cols; j++)
    {
      if( bf[j] ) { bottom = i; right = j; i = -1; break; }
    }
  }

  if( left > right )
  {
    i = left;
    left = right;
    right = i;
  }

  for(i=top; i<=bottom; i++)
  {
    bf = getRow(i);
    for(j=0; j<=left; j++)
    {
      if( bf[j] ) { left = Min(j, left); break; }
    }
  }

  for(i=top; i<=bottom; i++)
  {
    bf = getRow(i);
    for(j=cols-1; j>=right; j--)
    {
      if( bf[j] ) { right = Max(j, right); break; }
    }
  }
}

void DIBImage::Shift(int rowShift, int colShift)
{
  if( IsEmpty() ) return;

  DIBImage temp(*this);
  memset(buffer, 0, GetBuffLen());
  int i;
  unsigned char *buff1, *buff2;
  if( rowShift >= rows-1 || colShift >= cols-1 ) return;

  int desRStart = Max(0, rowShift);
  int desCStart = Max(0, colShift);
  int srcRStart = Max(0, -rowShift);
  int srcCStart = Max(0, -colShift);
  int rLen = rows - Abs(rowShift);
  int cLen = cols - Abs(colShift);
  int bytes = bitsPerPxl/8;


  for(i=desRStart; i<desRStart+rLen; i++)
  {
    buff1 = getRow(i);
    buff2 = temp.getRow(i-desRStart+srcRStart);

    memcpy(buff1+desCStart*bytes, buff2+srcCStart*bytes, cLen*bytes);
  }
}

void DIBImage::Rotate(const double degree, bool enlargeSize)
{
  if( IsEmpty() ) return;

  // Check if rotation is necessary
  double arc = degree*3.1415926535897932/180;
  if( fabs(sin(arc))*::Max(rows, cols) < 1 ) return ;

  const int SCALE = (1<<8);
  const int SCALE2 = (1<<16);
  DIBImage in(*this);
  int rowShift = 0, colShift = 0;
  int newRow = row(), newCol = col();
  int oldRow = row(), oldCol = col();

  if( enlargeSize )
  {
    newCol = Nearest(Max(0.0, Max(Max(oldRow*sin(arc), oldCol*cos(arc)+oldRow*sin(arc)), oldCol*cos(arc))) - 
                     Min(0.0, Min(Min(oldRow*sin(arc), oldCol*cos(arc)+oldRow*sin(arc)), oldCol*cos(arc))) );
    newRow = Nearest(Max(0.0, Max(Max(-oldCol*sin(arc), -oldCol*sin(arc)+oldRow*cos(arc)), oldRow*cos(arc))) - 
                     Min(0.0, Min(Min(-oldCol*sin(arc), -oldCol*sin(arc)+oldRow*cos(arc)), oldRow*cos(arc))) );
    rowShift = (newRow-oldRow)/2;
    colShift = (newCol-oldCol)/2;

    this->CreateImage(newRow, newCol, this->pxlbit());
  }
  else
  {
    this->CreateImage(oldRow, oldCol, this->pxlbit());
  }

  int i, j;
	short half_row = newRow/2;
	short half_col = newCol/2;
	short xx, yy;
	short cos_ = Nearest(cos(-arc)*SCALE);
	short sin_ = Nearest(sin(-arc)*SCALE);
	short alpha_x, alpha_y;
	int aa, bb;
  int x0, y0, x1, y1;

  if( in.pxlbit() == 8 )
  {
    unsigned char** rot = (unsigned char**) malloc(newRow*sizeof(unsigned char*));
    unsigned char** img = (unsigned char**) malloc(oldRow*sizeof(unsigned char*));
    for(i=0; i<newRow; i++)
    {
      rot[i] = getRow(i);
    }
    for(i=0; i<oldRow; i++)
    {
      img[i] = in.getRow(i);
    }

	  for( i = 0; i < newRow; i++)
    {
		  for( j = 0; j < newCol; j++)
		  {
			  x1 = j - half_col;
			  y1 = i - half_row;
			  x0 = x1*cos_ - y1*sin_;
			  y0 = x1*sin_ + y1*cos_;
			  aa = x0 + half_col*SCALE;
        bb = y0 + half_row*SCALE;
			  yy = bb/SCALE - rowShift;
        if( yy >= oldRow-1 || yy < 0 ) continue;
			  alpha_y = bb%SCALE;
			  xx = aa/SCALE - colShift;
        if( xx >= oldCol-1 || xx < 0) continue;
			  alpha_x = aa%SCALE;

				rot[i][j] = (img[yy][xx]*SCALE2 + alpha_x*(img[yy][xx+1]-img[yy][xx])*SCALE
								+alpha_y*((img[yy+1][xx]-img[yy][xx])*SCALE + alpha_x
								*(img[yy+1][xx+1]+img[yy][xx]-img[yy][xx+1]-img[yy+1][xx])))/SCALE2;
      }
    }

    free(rot);
    free(img);
  }
  else if( in.pxlbit() == 16 )
  {
    unsigned short** rot = (unsigned short**)malloc(newRow*sizeof(unsigned short*));
    unsigned short** img = (unsigned short**)malloc(oldRow*sizeof(unsigned short*));
    for(i=0; i<newRow; i++)
    {
      rot[i] = (unsigned short*)getRow(i);
    }
    for(i=0; i<oldRow; i++)
    {
      img[i] = (unsigned short*)in.getRow(i);
    }

	  for( i = 0; i < newRow; i++)
    {
		  for( j = 0; j < newCol; j++)
		  {
			  x1 = j - half_col;
			  y1 = i - half_row;
			  x0 = x1*cos_ - y1*sin_;
			  y0 = x1*sin_ + y1*cos_;
			  aa = x0 + half_col*SCALE;
        bb = y0 + half_row*SCALE;
			  yy = bb/SCALE - rowShift;
			  alpha_y = bb%SCALE;
        if( yy >= oldRow-1 || yy < 0 ) continue;
			  xx = aa/SCALE - colShift;;
			  alpha_x = aa%SCALE;
        if( xx >= oldCol-1 || xx < 0) continue;

				rot[i][j] = (img[yy][xx]*SCALE2 + alpha_x*(img[yy][xx+1]-img[yy][xx])*SCALE
								+alpha_y*((img[yy+1][xx]-img[yy][xx])*SCALE + alpha_x
								*(img[yy+1][xx+1]+img[yy][xx]-img[yy][xx+1]-img[yy+1][xx])))/SCALE2;
      }
    }

    free(rot);
    free(img);
  }
  else if( in.pxlbit() == 24 )
  {
    RGBColorStr** rot = (RGBColorStr**)malloc(newRow*sizeof(RGBColorStr*));
    RGBColorStr** img = (RGBColorStr**)malloc(oldRow*sizeof(RGBColorStr*));
    for(i=0; i<newRow; i++)
    {
      rot[i] = (RGBColorStr*)getRow(i);
    }
    for(i=0; i<oldRow; i++)
    {
      img[i] = (RGBColorStr*)in.getRow(i);
    }

	  for( i = 0; i < newRow; i++)
    {
		  for( j = 0; j < newCol; j++)
		  {
			  x1 = j - half_col;
			  y1 = i - half_row;
			  x0 = x1*cos_ - y1*sin_;
			  y0 = x1*sin_ + y1*cos_;
			  aa = x0 + half_col*SCALE;
        bb = y0 + half_row*SCALE;
			  yy = bb/SCALE - rowShift;
			  alpha_y = bb%SCALE;
        if( yy >= oldRow-1 || yy < 0 ) continue;
			  xx = aa/SCALE - colShift;;
			  alpha_x = aa%SCALE;
        if( xx >= oldCol-1 || xx < 0) continue;

				rot[i][j].B = (img[yy][xx].B*SCALE2 + alpha_x*(img[yy][xx+1].B-img[yy][xx].B)*SCALE
								+alpha_y*((img[yy+1][xx].B-img[yy][xx].B)*SCALE + alpha_x
								*(img[yy+1][xx+1].B+img[yy][xx].B-img[yy][xx+1].B-img[yy+1][xx].B)))/SCALE2;

				rot[i][j].G = (img[yy][xx].G*SCALE2 + alpha_x*(img[yy][xx+1].G-img[yy][xx].G)*SCALE
								+alpha_y*((img[yy+1][xx].G-img[yy][xx].G)*SCALE + alpha_x
								*(img[yy+1][xx+1].G+img[yy][xx].G-img[yy][xx+1].G-img[yy+1][xx].G)))/SCALE2;

				rot[i][j].R = (img[yy][xx].R*SCALE2 + alpha_x*(img[yy][xx+1].R-img[yy][xx].R)*SCALE
								+alpha_y*((img[yy+1][xx].R-img[yy][xx].R)*SCALE + alpha_x
								*(img[yy+1][xx+1].R+img[yy][xx].R-img[yy][xx+1].R-img[yy+1][xx].R)))/SCALE2;
      }
    }

    free(rot);
    free(img);
  }
}

void DIBImage::quantize2()
{
  int *histog = GetHistogram();
  if( histog == NULL ) return;

  double mean, dev;
  GetStat(mean, dev);
  double mean1, mean2;
  int index1, index2, i, index;

  unsigned char threshold1 = unsigned char(mean+0.5);
  unsigned char threshold2;
  index = 0;

  while(true)
  {
    index ++;
    mean1 = 0.;
    mean2 = 0.;
    index1 = 0;
    index2 = 0;

    for(i=0; i<threshold1; i++)
    {
      mean1 += histog[i]*i;
      index1 += histog[i];
    }
    for(i=threshold1; i<255; i++)
    {
      mean2 += histog[i]*i;
      index2 += histog[i];
    }

    if( index1 > 0 )
      mean1 /= index1;
    else
    {
    }
    if( index2 > 0 )
      mean2 /= index2;
    else
    {
    }

    threshold2 = unsigned char((mean1+mean2)/2+0.5);

    if( threshold1 != threshold2 )
    {
      threshold1 = threshold2;
      continue;
    }
    else
    {
      break;
    }

    if( index > 32 ) break;
  }

  unsigned char* buff;
  for(i=0; i<rows; i++)
  {
    buff = getRow(i);
    for(index=0; index<cols; index++)
    {
      if( buff[index] > threshold1 ) buff[index] = 255;
      else buff[index] = 0;
    }
  }

  free(histog);
}

void DIBImage::quantize3()
{
  int *histog = GetHistogram();
  if( histog == NULL ) return;

  double mean1, mean2, mean3;
  int i, index;
  int sum1, sum2, sum3;
  unsigned char threshold1, threshold2;
  unsigned char threshold1_, threshold2_;
  
  index = 0;
  sum1 = 0;
  sum2 = 0;
  while(index < 256)
  {
    sum1 += histog[index];
    sum2 += histog[index];
    if( sum1 >= rows*cols*0.333 && sum1 < rows*cols )
    {
      threshold1 = index;
      sum1 = rows*cols;
    }
    if( sum2 >= rows*cols*0.666 && sum2 < rows*cols )
    {
      threshold2 = index;
      break;
    }
    index++;
  }

  index = 0;
  while(true)
  {
    index ++;
    mean1 = 0.;
    mean2 = 0.;
    mean3 = 0.;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;

    for(i=0; i<threshold1; i++)
    {
      mean1 += histog[i]*i;
      sum1 += histog[i];
    }
    for(i=threshold1; i<threshold2; i++)
    {
      mean2 += histog[i]*i;
      sum2 += histog[i];
    }
    for(i=threshold2; i<255; i++)
    {
      mean3 += histog[i]*i;
      sum3 += histog[i];
    }

    if( sum1 > 0 )
      mean1 /= sum1;
    if( sum2 > 0 )
      mean2 /= sum2;
    if( sum3 > 0 )
      mean3 /= sum3;

    threshold1_ = unsigned char((mean1+mean2)/2+0.5);
    threshold2_ = unsigned char((mean2+mean3)/2+0.5);

    if( threshold1 != threshold1_ || threshold2 != threshold2_)
    {
      threshold1 = threshold1_;
      threshold2 = threshold2_;
      continue;
    }
    else
    {
      break;
    }

    if( index > 32 ) break;
  }

  unsigned char* buff;
  for(i=0; i<rows; i++)
  {
    buff = getRow(i);
    for(index=0; index<cols; index++)
    {
      if( buff[index] >= threshold2 ) buff[index] = 255;
      else if( buff[index] >= threshold1 && buff[index] < threshold2 ) buff[index] = 128;
      else buff[index] = 0;
    }
  }

  free(histog);
}

void DIBImage::quantize4()
{
  int *histog = GetHistogram();
  if( histog == NULL ) return;

  double mean1, mean2, mean3, mean4;
  int i, index;
  int sum1, sum2, sum3, sum4;
  unsigned char threshold1, threshold2, threshold3;
  unsigned char threshold1_, threshold2_, threshold3_;
  
  index = 0;
  sum1 = 0;
  sum2 = 0;
  sum3 = 0;
  while(index < 256)
  {
    sum1 += histog[index];
    sum2 += histog[index];
    sum3 += histog[index];
    if( sum1 >= rows*cols*0.25 && sum1 < rows*cols )
    {
      threshold1 = index;
      sum1 = rows*cols;
    }
    if( sum2 >= rows*cols*0.50 && sum2 < rows*cols )
    {
      threshold2 = index;
      sum2 = rows*cols;
    }
    if( sum3 >= rows*cols*0.750 && sum3 < rows*cols )
    {
      threshold3 = index;
      break;
    }
    index++;
  }

  index = 0;
  while(true)
  {
    index ++;
    mean1 = 0.;
    mean2 = 0.;
    mean3 = 0.;
    mean4 = 0.;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;

    for(i=0; i<threshold1; i++)
    {
      mean1 += histog[i]*i;
      sum1 += histog[i];
    }
    for(i=threshold1; i<threshold2; i++)
    {
      mean2 += histog[i]*i;
      sum2 += histog[i];
    }
    for(i=threshold2; i<threshold3; i++)
    {
      mean3 += histog[i]*i;
      sum3 += histog[i];
    }
    for(i=threshold3; i<255; i++)
    {
      mean4 += histog[i]*i;
      sum4 += histog[i];
    }

    if( sum1 > 0 )
      mean1 /= sum1;
    if( sum2 > 0 )
      mean2 /= sum2;
    if( sum3 > 0 )
      mean3 /= sum3;
    if( sum4 > 0 )
      mean4 /= sum4;

    threshold1_ = unsigned char((mean1+mean2)/2+0.5);
    threshold2_ = unsigned char((mean2+mean3)/2+0.5);
    threshold3_ = unsigned char((mean3+mean4)/2+0.5);

    if( threshold1 != threshold1_ || threshold2 != threshold2_ || threshold3 != threshold3_ )
    {
      threshold1 = threshold1_;
      threshold2 = threshold2_;
      threshold3 = threshold3_;
      continue;
    }
    else
    {
      break;
    }

    if( index > 32 )
      break;
  }

  unsigned char* buff;
  for(i=0; i<rows; i++)
  {
    buff = getRow(i);
    for(index=0; index<cols; index++)
    {
      if( buff[index] >= threshold3 ) buff[index] = 255;
      else if( buff[index] >= threshold2 && buff[index] < threshold3 ) buff[index] = 172;
      else if( buff[index] >= threshold1 && buff[index] < threshold2 ) buff[index] = 64;
      else buff[index] = 0;
    }
  }

  free(histog);
}

// operations for color models
void RGBColorStr::ToYIQ( YIQColorStr* yiq )
{
  yiq->Y = 0.299*R+0.587*G+0.114*B;
  yiq->I = 0.596*R-0.275*G-0.321*B;
  yiq->Q = 0.212*R-0.523*G+0.311*B;
}

void RGBColorStr::ToHSI( HSIColorStr* hsi )
{
  static const double pi = 3.1415926535897932;
  hsi->I = (R+G+B)/3.;
  if( hsi->I > 0. )
  {
    unsigned char v = (R>G)? G : R;
    hsi->S = 1-((v>B)?B:v)/hsi->I;
  }
  else hsi->S = 0;

  if( R==G && R == B) hsi->H = 0.;
  else hsi->H = acos( (R-(G+B)*0.5)/sqrt((double)(R-G)*(R-G)+(R-B)*(G-B)));
  
  if( B > G ) hsi->H = 2*pi - hsi->H;
}

void YIQColorStr::ToRGB(RGBColorStr * rgb)
{
  rgb->R = Min(Max(int(Y+0.955688*I+0.619858*Q+0.5), 0), 255);
  rgb->G = Min(Max(int(Y-0.271582*I-0.646874*Q+0.5), 0), 255);
  rgb->B = Min(Max(int(Y-1.108177*I+1.705065*Q+0.5), 0), 255);
}

void HSIColorStr::ToRGB(RGBColorStr * rgb)
{
  static const double pi = 3.1415926535897932;
  double h;
  int temp;

  if( H > 0. && H <= pi*2/3 )
  {
    temp = int(I*(1+S*cos(H)/cos(pi/3-H))+0.5);
    rgb->R = Min(Max(temp, 0), 255);
    temp = int(I*(1-S)+0.5);
    rgb->B = Min(Max(temp, 0), 255);
    if( I > 0. )
      temp = int(3*I*(1-(rgb->R+rgb->B)/(3*I))+0.5);
    else
      temp = 0;
    rgb->G = Min(Max(temp, 0), 255);
  }
  else if( H > pi*2/3 && H <= pi*4/3 )
  {
    h = H-pi*2/3;
    temp = int(I*(1+S*cos(h)/cos(pi/3-h))+0.5);
    rgb->G = Min(Max(temp, 0), 255);
    temp = int(I*(1-S)+0.5);
    rgb->R = Min(Max(temp, 0), 255);
    if( I > 0. )
      temp = int(3*I*(1-(rgb->R+rgb->G)/(3*I))+0.5);
    else
      temp = 0;
    rgb->B = Min(Max(temp, 0), 255);
  }
  else
  {
    h = H-pi*4/3;
    temp = int(I*(1+S*cos(h)/cos(pi/3-h))+0.5);
    rgb->B = Min(Max(temp, 0), 255);
    temp = int(I*(1-S)+0.5);
    rgb->G = Min(Max(temp, 0), 255);
    if( I > 0. )
      temp = int(3*I*(1-(rgb->B+rgb->G)/(3*I))+0.5);
    else temp = 0;
    rgb->R = Min(Max(temp, 0), 255);
  }
}

