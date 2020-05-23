#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <cmath>
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

void gauss_extension(float s, int HSize, float* H){
  int i, j;
	double cst = 1./(2.0*pi*s*s);
	double e;
	for(i=0;i<HSize;i++){
		for(j=0;j<HSize;j++){
			float x=(float)(i-HSize/2);
			float y =(float)(j-HSize/2);
			e = exp(-(x*x+y*y)/(2*s*s));
			H[i*HSize+j] = e * cst;
		}
	}
	double sum = 0.0;
	for(i=0;i<HSize;i++){
		for(j=0;j<HSize;j++){
			sum += H[j+i*HSize];
		}
	}
	for(i=0;i<HSize;i++){
		for(j=0;j<HSize;j++){
			H[j+i*HSize] /= sum;
		}
	}
}

void apply2dMask(float* mask, int** image, int** img, int msize, int N, int M){
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
			//cout << "test" << endl;
      img[i][j] = sum;
      index = 0;
      sum = 0;
    }
  }
}

void gaus2d_wrapper(float s, int Hsize, float* H, int** image, int** img, int N, int M){
  H = new float[Hsize * Hsize];
	//	cout << "here" << endl;
  gauss_extension(s, Hsize, H);
	apply2dMask(H, image, img, Hsize, N, M);

}

void substract(int** a, int** b, int** c, int N, int M){
    int i, j;
    for(i = 0; i < N; i++){
      for(j = 0; j < M; j++){
        c[i][j] = b[i][j] - a[i][j];
      }
    }
}

void downsample(int** img, int** newimg, int s, int N, int M){
  int i,j, index_x, index_y;
  //newimg = new int*[N/s];
	//cout << s << " " << s << endl;
  //for(i = 0; i < N/s; i++){
  //  newimg[i] = new int[M/s];
  //}
  index_x = 0;
  index_y = 0;
  for(i = 0; i < N; i+=s){
    for(j = 0; j < M; j+=s){
      newimg[index_x][index_y] = img[i][j];
			index_y++;
    }
		index_y = 0;
		index_x++;
  }
	//cout << "here" << endl;
}

int main(int argc, char** argv)
{
    int i, j;
    int s, n;
    double k, t, sa, sigma;
    //n = (int)(argv[1][0] - '0');//levels
    //s = (int)(argv[2][0] - '0');//intermediate
    //sigma = (int)(argv[3][0] - '0');//scale
		cout << "Enter the amount of octaves:";
		cin >> n;
		cout << "Enter the amount of intermediate levels:";
		cin >> s;
		cout << "Enter the scale(sigma):";
		cin >> sigma;
		s++;
    k = 1;
    int ** input, ** output, ** difference; //input and output images
    int x_size, y_size, Q;
    char name[20] = "lenna.pgm";
		char outfile[30] = "dog_lenna_1.pgm";
    char out[30] = "lenna_1.pgm";
    float** mask1;
    mask1 = new float*[s];
    ReadImage(name, &input, x_size, y_size, Q);
		ReadImage(name, &output, x_size, y_size, Q);
    ReadImage(name, &difference, x_size, y_size, Q);
    WriteImage(out, input, x_size, y_size, Q);
		int x = x_size;
		int y = y_size;
		char ic[2];
		char jc[2];
		//cout << i << j << endl;
		ic[0] = 'a';
		ic[1] = '\0';
		jc[0] = 'b';
		jc[1] = '\0';
    //gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, input_1, x_size, y_size);
    for(i = 0; i < n; i++){
			//WriteImage(outfile, output, x, y, Q);
      for(j = 0; j <= s; j++){
        if(j != 0){
					sa = (double)(s+4)-j;
					//cout << sa << endl;
					t = 1.0/sa;
					k = pow(2.0, t);
				}
        //ReadImage(out, &input, x, y, Q);
        int temp = k * sigma;
				cout << k * sigma << " " << temp << endl;
        gaus2d_wrapper(k*sigma, 5*temp, mask1[j], input, output, x, y);

        strcpy(outfile, "a_lenna_");
        std::strncat(outfile, ic,1);
        std::strcat(outfile, "_");
        std::strncat(outfile, jc,1);
        std::strcat(outfile, ".pgm");
        WriteImage(outfile, output, x, y, Q);
        //k = pow(2, t);
        if(j > 0){
          substract(input, output, difference, x, y);
          strcpy(outfile, "dog_lenna_");
  				std::strncat(outfile, ic,1);
  				std::strcat(outfile, "_");
  				std::strncat(outfile, jc,1);
  				std::strcat(outfile, ".pgm");
          WriteImage(outfile, difference, x, y, Q);
        }
        strcpy(outfile, "a_lenna_");
        std::strncat(outfile, ic,1);
        std::strcat(outfile, "_");
        std::strncat(outfile, jc,1);
        std::strcat(outfile, ".pgm");
        ReadImage(outfile, &input, x, y, Q);
        //output = input;
				ic[0]++;
				jc[0]++;
      }
			k = 1;
      downsample(output, input, 2, x, y);
      x/=2;
      y/=2;
      //WriteImage(out, input, x, y, Q);
    }

    //THIS IS THE END OF B
    return 0;
}
