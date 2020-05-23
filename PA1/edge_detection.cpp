#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include "stdlib.h"

using namespace std;

#define pi 3.141592653589793



/* Create a one-dimensional Gaussian mask
 "H" of scale "s" (sigma) on "Hsize" pixels.

   The values of the mask are scaled so that the
 sum of mask values = 1.0
*/



void Gauss (float s, int Hsize, float * H)

{

	int     i;

	float  cst,tssq,x,sum;



	cst=1./(s*sqrt(2.0*pi));

	tssq=1./(2*s*s);


	for (i=0; i<Hsize; i++)
	{

		x=(float)(i-Hsize/2);

		H[i]=(cst*exp(-(x*x*tssq)));

	}


	sum=0.0;

	for (i=0;i<Hsize;i++)

		sum += H[i];

	for(i=0;i<Hsize;i++)

	H[i] /= sum;

}


void WriteImage(char fname[], int **fimage, int M, int N, int Q)
{
 int i, j;
 unsigned char *image;
 ofstream ofp;

 image = (unsigned char *) new unsigned char [M*N];

 // convert the integer values to unsigned char

 for(i=0; i<N; i++)
   for(j=0; j<M; j++)
     image[i*M+j]=(unsigned char)fimage[i][j];

 ofp.open(fname, ios::out);

 if (!ofp) {
   cout << "Can't open file: " << fname << endl;
   exit(1);
 }

 ofp << "P5" << endl;
 ofp << M << " " << N << endl;
 ofp << Q << endl;

 ofp.write( reinterpret_cast<char *>(image), (M*N)*sizeof(unsigned char));

 if (ofp.fail()) {
   cout << "Can't write image " << fname << endl;
   exit(0);
 }

 ofp.close();

}

void ReadImage(char fname[], int ***fimage, int& M, int& N, int& Q)
{
 int i, j;
 unsigned char *image;
 char header [100], *ptr;
 ifstream ifp;

 ifp.open(fname, ios::in);

 if (!ifp) {
   cout << "Can't read image: " << fname << endl;
   exit(1);
 }

 // read header

 ifp.getline(header,100,'\n');
 if ( (header[0]!=80) ||    /* 'P' */
      (header[1]!=53) ) {   /* '5' */
      cout << "Image " << fname << " is not PGM" << endl;
      exit(1);
 }

 ifp.getline(header,100,'\n');
 while(header[0]=='#')
   ifp.getline(header,100,'\n');

 M=strtol(header,&ptr,0);
 N=atoi(ptr);

 ifp.getline(header,100,'\n');

 Q=strtol(header,&ptr,0);

 image = (unsigned char *) new unsigned char [M*N];

 *fimage = new int* [N];
 for(i=0; i<N; i++)
   (*fimage)[i] = new int[M];

 ifp.read( reinterpret_cast<char *>(image), (M*N)*sizeof(unsigned char));

 if (ifp.fail()) {
   cout << "Image " << fname << " has wrong size" << endl;
   exit(1);
 }

 ifp.close();

 //
 // Convert the unsigned characters to integers
 //

 for(i=0; i<N; i++)
   for(j=0; j<M; j++)
     (*fimage)[i][j]=(int)image[i*M+j];
}

void mul(float* A, float* B, float* C, int I, int J){
  int i, j;
  for(i = 0; i < I; i++){
    for(j = 0; j < J; j++){
      C[i*I+j] = A[i]*B[j];
    }
  }
}

void gauss_extension(float s, int HSize, float* Hx, float* Hy, float* H){
  Gauss(s,HSize,Hx);//generates first gauss 1d kernel
  Gauss(s,HSize,Hy);//generates second gauss 1d kernel
  mul(Hx,Hy,H,HSize,HSize); //multiplies them together
}

void normalize(int** image, int N, int M){
  int i, j;
  int imax = -1000000;
  int imin = 10000000;
  double pix;
	int oldpix;
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
      imax = max(imax, image[i][j]);
      imin = min(imin, image[i][j]);
    }
  }
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
			oldpix = image[i][j];
      pix = (double)((oldpix-imin)*255)/(double)(imax-imin);
      image[i][j] = (int)pix;
    }
  }

}

void apply2dMask(float mask[], int** image, int** newimg, int msize, int N, int M){
  int i, j, k, l;
  int index = 0;
  float sum = 0;

  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
      for(k = (-msize/2); k <= msize/2; k++ ){
        for(l = (-msize/2); l <= msize/2; l++){
          if((i-k)<0 || (i-k)>=N || (j-l)<0 || (j-l)>=M){
            sum+=0;
          }
          else{
            sum += image[i-k][j-l] * mask[index];
          }
          index++;
        }
      }
      newimg[i][j] = sum;
      index = 0;
      sum = 0;
    }
  }
}

void magn(int** himage, int** vimage, int** mimage, int N, int M){
  int i, j;
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
      mimage[i][j] = sqrt((himage[i][j] * himage[i][j]) + (vimage[i][j] * vimage[i][j]));
    }
  }
}

void thresh(int** mimage, int** timage, int N, int M, int x){
  int i, j;
  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
      if(mimage[i][j] < x){
        timage[i][j] = 0;
      }
      else{
        timage[i][j] = 255;
      }
    }
  }
}

void gaus2d_wrapper(float s, int Hsize, float* H, float* Hx, float* Hy, int** image, int N, int M){
  gauss_extension(s, Hsize, Hx, Hy, H);
  apply2dMask(H, image, image, Hsize, N, M);
}



void sobel_wrapper(float hmask[], float vmask[], int** himage, int** vimage, int** nhimage, int** nvimage, int** mimage, int** timage, int msize, int N, int M, int x){
  apply2dMask(hmask, himage,nhimage, msize, N, M);
  apply2dMask(vmask, vimage,nvimage, msize, N, M);
	magn(nhimage, nvimage, mimage, N, M);
	normalize(nhimage,N,M);
	normalize(nvimage,N,M);
	normalize(mimage,N,M);
  thresh(mimage, timage, N, M, x);
  //normalize(timage,N,M);
}

int main(int argc, char** argv)
{
		cout << "Pick a threshold size:";
    int x;//(int)(argv[1][0] - '0');
		cin >> x;
    float vsobel[9] = {-1,-2,-1,
                      0,0,0,
                      1,2,1};
    float hsobel[9] = {-1,0,1
                     ,-2,0,2,
                     -1,0,2};
    int ** hinput_lenna, ** vinput_lenna,**nhinput_lenna, **nvinput_lenna, ** minput_lenna, ** tinput_lenna; //input and output images
    int ** hinput_sf, ** vinput_sf,**nhinput_sf, **nvinput_sf, ** minput_sf, ** tinput_sf;
    int lenna_x_size, lenna_y_size, Q_lenna, sf_x_size, sf_y_size, Q_sf;
    char lenna_name[20] = "lenna.pgm";
    char sf_name[20] = "sf.pgm";
    char lenna_outfile_x[20] = "lennaout_x.pgm";
    char lenna_outfile_y[20] = "lennaout_y.pgm";
    char lenna_outfile_magn[20] = "lennaout_magn.pgm";
    char lenna_outfile_thresh[20] = "lennaout_thresh.pgm";
    char sf_outfile_x[20] = "sfout_x.pgm";
    char sf_outfile_y[20] = "sfout_y.pgm";
    char sf_outfile_magn[20] = "sfout_magn.pgm";
    char sf_outfile_thresh[20] = "sfout_thresh.pgm";
    float* mask1x;
    float* mask1y;
    float* mask1_2d;
    mask1x = new float[5];
    mask1y = new float[5];
    mask1_2d = new float[25];
    ReadImage(lenna_name, &hinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    ReadImage(lenna_name, &vinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    ReadImage(lenna_name, &minput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    ReadImage(lenna_name, &tinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
		ReadImage(lenna_name, &nhinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
		ReadImage(lenna_name, &nvinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);

    ReadImage(sf_name, &hinput_sf, sf_x_size, sf_y_size, Q_sf);
    ReadImage(sf_name, &vinput_sf, sf_x_size, sf_y_size, Q_sf);
    ReadImage(sf_name, &minput_sf, sf_x_size, sf_y_size, Q_sf);
    ReadImage(sf_name, &tinput_sf, sf_x_size, sf_y_size, Q_sf);
		ReadImage(sf_name, &nhinput_sf, sf_x_size, sf_y_size, Q_sf);
		ReadImage(sf_name, &nvinput_sf, sf_x_size, sf_y_size, Q_sf);
   	gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, hinput_lenna, lenna_x_size, lenna_y_size);
    gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, vinput_lenna, lenna_x_size, lenna_y_size);
    gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, hinput_sf, sf_x_size, sf_y_size);
    gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, vinput_sf, sf_x_size, sf_y_size);

    sobel_wrapper(hsobel, vsobel, hinput_lenna, vinput_lenna, nhinput_lenna, nvinput_lenna, minput_lenna, tinput_lenna, 3, lenna_x_size, lenna_y_size, x);
    sobel_wrapper(hsobel, vsobel, hinput_sf, vinput_sf, nhinput_sf, nvinput_sf, minput_sf, tinput_sf, 3, sf_x_size, sf_y_size, x);


    WriteImage(lenna_outfile_x, nhinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    WriteImage(lenna_outfile_y, nvinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    WriteImage(lenna_outfile_magn, minput_lenna, lenna_x_size, lenna_y_size, Q_lenna);
    WriteImage(lenna_outfile_thresh, tinput_lenna, lenna_x_size, lenna_y_size, Q_lenna);

    WriteImage(sf_outfile_x, nhinput_sf, sf_x_size, sf_y_size, Q_sf);
    WriteImage(sf_outfile_y, nvinput_sf, sf_x_size, sf_y_size, Q_sf);
    WriteImage(sf_outfile_magn, minput_sf, sf_x_size, sf_y_size, Q_sf);
    WriteImage(sf_outfile_thresh, tinput_sf, sf_x_size, sf_y_size, Q_sf);

    //this is the end of c
    return 0;
}
