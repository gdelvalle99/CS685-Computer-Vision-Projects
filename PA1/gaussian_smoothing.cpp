#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>

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

void apply1dMask(float* mask, double* image, double* img,int msize, int isize){
  int i, j;
  int index = 0;
  double sum = 0;

  for(i = 0; i < isize; i++){
    for(j = (-msize/2); j <= msize/2; j++){
      if((i-j)<0 || (i-j)>=isize){
        sum += 0;
      }
      else{

        sum += (image[i-j]*mask[index]);
      }
      index++;
    }
    //cout << sum << endl;
		cout << sum << endl;
    img[i] = sum;
    sum = 0;
    index = 0;
  }
}

void apply1dhMask(float* mask, int** image, int msize, int N, int M){
  int i, j, k;
  int index = 0;
  float sum = 0;
/*  for(i = 0; i < msize; i++){
    sum += mask[i];
  }
  for(i = 0; i < msize; i++){
    mask[i] /= sum;
  }*/

  for(i = 0; i < N; i++){
    for(j = 0; j < M; j++){
      for(k = (-msize/2); k <= msize/2; k++){
        if((j-k)<0 || (j-k)>=M){
          sum += 0;
        }
        else{
          sum += (image[i][j-k]*mask[index]);
        }
        index++;
      }
      image[i][j] = sum;
      sum = 0;
      index = 0;
    }
  }
}

void apply1dvMask(float* mask, int** image, int msize, int N, int M){
  int i, j, k;
  int index = 0;
  float sum = 0;
/*  for(i = 0; i < msize; i++){
    sum += mask[i];
  }
  for(i = 0; i < msize; i++){
    mask[i] /= sum;
  }*/

  for(j = 0; j < M; j++){
    for(i = 0; i < N; i++){
      for(k = (-msize/2); k <= msize/2; k++){
        if((i-k)<0 || (i-k)>=N){
          sum += 0;
        }
        else{
          sum += (image[i-k][j]*mask[index]);
        }
        index++;
      }
      image[i][j] = sum;
      sum = 0;
      index = 0;
    }
  }
}

void apply2dMask(float* mask, int** image, int** newimg, int msize, int N, int M){
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

void gaus1d_wrapper(float s, int HSize, float* H, double* image, int isize, char file[]){
  double* img;
  img = new double[isize];
  Gauss(s, HSize, H);
  apply1dMask(H, image, img,HSize, isize);
	/*for(int j = 0; j < isize; j++){
    img[j] = image[j];
  }*/
  ofstream out;
  out.open(file, ios::out);
  for(int i = 0; i < isize; i++){
    out << img[i] << endl;
  }
}

void gaus2d_wrapper(float s, int Hsize, float* H, float* Hx, float* Hy, int** image, int** newimg,int N, int M){
  gauss_extension(s, Hsize, H);
  apply2dMask(H, image, newimg, Hsize, N, M);
}

void gaus2d1d_wrapper(float s,int HSize,float* H, int** image, int N, int M){
  Gauss(s, HSize, H);
  apply1dhMask(H, image, HSize, N, M);
  apply1dvMask(H, image, HSize, N, M);
}


int main(int argc, char** argv)
{
    //int ** input, ** output; //input and output images
    //int x_size, y_size, Q;
    //char name[20] = "lenna.pgm";
    //char outfile[20] = "lennaout.pgm";

    //ReadImage(name, &input, x_size, y_size, Q);

    //WriteImage(outfile, input, x_size, y_size, Q);
    //THIS IS A
    char file1[20] = "rect1.txt";
    char file5[20] = "rect5.txt";
    char file11[20] = "rect11.txt";
    float* mask1;
    float* mask5;
    float* mask11;
    mask1 = new float[5];
    mask5 = new float[25];
    mask11 = new float[55];
    double* rect1;
		double* rect5;
		double* rect11;
    double number;
    int index = 0;
    rect1 = new double[128];
		rect5 = new double[128];
		rect11 = new double[128];
    ifstream ifs("Rect_128.txt", ifstream::in);
    while(ifs >> number){
      rect1[index] = number;
			rect5[index] = number;
			rect11[index] = number;
			//cout << rect11[index] << endl;
      index++;
    }
    gaus1d_wrapper(1, 5, mask1, rect1, 128, file1);
    gaus1d_wrapper(5,25,mask5,rect5,128, file5);
    gaus1d_wrapper(11,55,mask11,rect11,128, file11);

    //THIS IS THE END OF A

    //THIS IS B
    int ** input_1, ** input_5, ** input_11, ** output1, ** output5, ** output11; //input and output images
    int x_size, y_size, Q;
    char name[20] = "lenna.pgm";
    char outfile1[20] = "lennaout1.pgm";
    char outfile5[20] = "lennaout5.pgm";
    char outfile11[20] = "lennaout11.pgm";
    float* mask1x;
    float* mask5x;
    float* mask11x;
    float* mask1y;
    float* mask5y;
    float* mask11y;
    float* mask1_2d;
    float* mask5_2d;
    float* mask11_2d;
    mask1x = new float[5];
    mask5x = new float[25];
    mask11x = new float[55];
    mask1y = new float[5];
    mask5y = new float[25];
    mask11y = new float[55];
    mask1_2d = new float[25];
    mask5_2d = new float[625];
    mask11_2d = new float[3025];
    ReadImage(name, &input_1, x_size, y_size, Q);
    ReadImage(name, &input_5, x_size, y_size, Q);
    ReadImage(name, &input_11, x_size, y_size, Q);
		ReadImage(name, &output1, x_size, y_size, Q);
		ReadImage(name, &output5, x_size, y_size, Q);
		ReadImage(name, &output11, x_size, y_size, Q);
    gaus2d_wrapper(1, 5, mask1_2d, mask1x, mask1y, input_1, output1, x_size, y_size);
    gaus2d_wrapper(5, 25, mask5_2d, mask5x, mask5y, input_5, output5, x_size, y_size);
    gaus2d_wrapper(11, 55, mask11_2d, mask11x, mask11y, input_11, output11, x_size, y_size);
    WriteImage(outfile1, output1, x_size, y_size, Q);
    WriteImage(outfile5, output5, x_size, y_size, Q);
    WriteImage(outfile11, output11, x_size, y_size, Q);

    //THIS IS THE END OF B

    //THIS IS C
		ReadImage(name, &input_1, x_size, y_size, Q);
    ReadImage(name, &input_5, x_size, y_size, Q);
    ReadImage(name, &input_11, x_size, y_size, Q);
    char outfile1c[20] = "lennaout1c.pgm";
    char outfile5c[20] = "lennaout5c.pgm";
    char outfile11c[20] = "lennaout11c.pgm";
    gaus2d1d_wrapper(1,5,mask1, input_1, x_size, y_size);
    gaus2d1d_wrapper(5,25,mask5, input_5, x_size, y_size);
    gaus2d1d_wrapper(11,55,mask11, input_11, x_size, y_size);
    WriteImage(outfile1c, input_1, x_size, y_size, Q);
    WriteImage(outfile5c, input_5, x_size, y_size, Q);
    WriteImage(outfile11c, input_11, x_size, y_size, Q);

    //this is the end of c
    return 0;
}
