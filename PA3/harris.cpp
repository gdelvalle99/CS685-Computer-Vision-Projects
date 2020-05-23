#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <opencv2/opencv.hpp>
using namespace std;

#define pi 3.141592653589793

void Gauss_deriv1 (float s, int Hsize, float * H)
{
  int     i;
  float  x,cst,tssq;

  cst = 1.0/(s*s*s*sqrt(2.0*pi)) ;
  tssq = -1.0/(2.0*s*s) ;

  for (i=0; i<Hsize; i++)
  {
    x = (float)(i-Hsize/2);
    H[i] = -x*cst*exp(x*x*tssq);
  }
}

void Gauss_derivy(float s, int Hsize, float* H){
  float x,y;
  for(int i = 0; i < Hsize; i++){
    for(int j = 0; j < Hsize; j++){
      x = (float)(i-Hsize/2);
      y = (float)(j-Hsize/2);
      H[i*Hsize+j] = -(y/(2*pi*s*s*s*s))*exp(-(((x*x)+(y*y))/(2*s*s)));
    }
  }
}

void Gauss_derivx(float s, int Hsize, float* H){
  float x,y;
  for(int i = 0; i < Hsize; i++){
    for(int j = 0; j < Hsize; j++){
      x = (float)(i-Hsize/2);
      y = (float)(j-Hsize/2);
      H[i*Hsize+j] = -(x/(2*pi*s*s*s*s))*exp(-(((x*x)+(y*y))/(2*s*s)));
    }
  }
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


int getSum(int** mat, int N, int M, int i, int j, cv::Mat src){
  int sum = 0;
  int max = 0;
  int ti, tj;
  for(int k = (-3/2); k <= 3/2; k++ ){
    for(int l = (-3/2); l <= 3/2; l++){
      if((i-k)<0 || (i-k)>=N || (j-l)<0 || (j-l)>=M){
        sum+=0;
      }
      else{
        sum += mat[i-k][j-l];
        //max = std::max(mat[i-k][j-l],max);
        if(max < mat[i-k][j-l]){
          max = mat[i-k][j-l];
          ti = i-k;
          tj = j-l;
        }
      }
    }
  }
  cv::Point center(tj,ti);
  cv::circle(src,center,1,cv::Scalar(0,125,0));
  return sum;
}
double getDeterminant(int a, int b, int c, int d){
  return (a * d) - (b * c);
}

double getTrace(int a, int b){
  return a + b;
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


int main(){
  char out4[20] = "out.pgm";
  char name[20] = "Test3.pgm";
  char out1[20] = "try7.pgm";
  char out2[20] = "try8.pgm";
  char out3[20] = "try9.pgm";
  cv::Mat src;
  src = cv::imread(name);
  int** fx,**out, **fx_new, **fy, **fy_new, **fxfy, **fxfy_new,**output, **fx_true, **fy_true, **fxfy_true;
  int x_size, y_size, Q;
  float* maskx = new float[9];
  float* masky = new float[9];
  float* mask_2d = new float[9];
  gauss_extension(1.5,3,mask_2d);
  //Gauss_derivx(1.05,3,mask);
  Gauss_derivx(1.05,3,maskx);
  Gauss_derivy(1.05,3,masky);
  for(int i = 0;  i < 3; i++){
    for(int j = 0; j < 3; j++){
      std::cout << maskx[i*3+j]<< " ";
    }
    std::cout << std::endl;
    //std::cout << mask[i] << std::endl;
  }
  for(int i = 0;  i < 3; i++){
    for(int j = 0; j < 3; j++){
      std::cout << masky[i*3+j]<< " ";
    }
    std::cout << std::endl;
    //std::cout << mask[i] << std::endl;
  }
  ReadImage(name,&out,x_size,y_size,Q);
  ReadImage(name,&fx,x_size,y_size,Q);
  ReadImage(name,&fy,x_size,y_size,Q);
  ReadImage(name,&fxfy,x_size,y_size,Q);
  ReadImage(name,&fx_new,x_size,y_size,Q);
  ReadImage(name,&fy_new,x_size,y_size,Q);
  ReadImage(name,&fxfy_new,x_size,y_size,Q);
  ReadImage(name,&fx_true,x_size,y_size,Q);
  ReadImage(name,&fy_true,x_size,y_size,Q);
  ReadImage(name,&fxfy_true,x_size,y_size,Q);
  std::cout << "here" << std::endl;
  //apply1dhMask(mask,fx,3,y_size,x_size);
  std::cout << "here" << std::endl;
  //apply1dhMask(mask,fx,5,x_size,y_size);
  //apply1dvMask(mask,fy,3,y_size,x_size);
  //apply1dvMask(mask,fy,5,x_size,y_size);
  //apply1dhMask(mask,fxfy,5,x_size,y_size);
  //apply1dvMask(mask,fxfy,5,x_size,y_size);
  apply2dMask(maskx,fx,fx_new,3,y_size,x_size);
  apply2dMask(masky,fy,fy_new,3,y_size,x_size);
  //apply2dMask(mask_2d,fxfy,fxfy_new,3,y_size,x_size);
  std::cout << "here" << std::endl;
//  WriteImage(out1,fx,x_size,y_size,Q);
//  WriteImage(out2,fy,x_size,y_size,Q);
//  WriteImage(out3,fxfy,x_size,y_size,Q);

  for(int i = 0; i <y_size; i++){
    for(int j = 0; j < x_size;j++){
      fx_new[i][j] = fx_new[i][j] * fx_new[i][j];
      fy_new[i][j] = fy_new[i][j] * fy_new[i][j];
      fxfy_new[i][j] = fx_new[i][j] * fy_new[i][j];
    }
  }

  apply2dMask(mask_2d,fx_new,fx_true,3,y_size,x_size);
  apply2dMask(mask_2d,fy_new,fy_true,3,y_size,x_size);
  apply2dMask(mask_2d,fxfy_new,fxfy_true,3,y_size,x_size);
  WriteImage(out1,fx_true,x_size,y_size,Q);
  WriteImage(out2,fy_true,x_size,y_size,Q);
  WriteImage(out3,fxfy_true,x_size,y_size,Q);
  //int** truex = new truex*[x_size*2];
  int fxd, fyd, fxfyd;
  double det;
  double corn_max = 0;
  double trace;
  double corner;
  double alpha = 0.06;
  double threshold = 0.01;
for(int i = 0; i < y_size; i++){
    for(int j = 0; j < x_size; j++){
      fxd = getSum(fx_true,y_size,x_size,i,j,src);
      fyd = getSum(fy_true,y_size,x_size,i,j,src);
      fxfyd = getSum(fxfy_true,y_size,x_size,i,j,src);
      det = getDeterminant(fxd,fxfyd,fxfyd,fyd);
      trace = getTrace(fxd,fyd);
      corner = det - alpha*trace*trace;
      corn_max = std::max(corn_max,corner);
    }
  }
std::cout << corn_max << std::endl;
for(int i = 0; i < y_size; i++){
  for(int j = 0; j < x_size; j++){
    //fxd = fx_new[i][j];//getSum(fx_new,x_size,y_size,i,j);
    //fyd = fy_new[i][j];//getSum(fy_new,x_size,y_size,i,j);
  //  fxfyd = fxfy_new[i][j];
    fxd = getSum(fx_true,y_size,x_size,i,j,src);
    fyd = getSum(fy_true,y_size,x_size,i,j,src);
    fxfyd = getSum(fxfy_true,y_size,x_size,i,j,src);
    det = getDeterminant(fxd,fxfyd,fxfyd,fyd);
    trace = getTrace(fxd,fyd);
    corner = det - alpha*trace*trace;
    if(corner > threshold*corn_max){
      //out[i][j] = 255;
      cv::Point center(j,i);
      cv::circle(src,center,5,(255));
    }
  }
}


  //cv::Mat M(x_size,y_size,CV_32SC3,cv::Scalar(0,0,255));
  WriteImage(out4,out,x_size,y_size,Q);
  cv::imshow("test",src);
  cv::imwrite("Harris3.jpg",src);
  cv::waitKey(0);
  return 0;
  //compute the horizontal and vertical derivatives of the image
  //computer the three images from the horizontal/vertical derivatives
  //convolve each of theese images with a larger Gaussian

  //computer scalar interest measure using one of the formulas discusses (r(aw)=det(aw)-alpha*trace^2(aw))
  //find local maxima above a certain threshold and report these as detected feature points
}
