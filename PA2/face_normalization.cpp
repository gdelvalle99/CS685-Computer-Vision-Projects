#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>

//using namespace std;

#define pi 3.141592653589793



/* Create a one-dimensional Gaussian mask
 "H" of scale "s" (sigma) on "Hsize" pixels.

   The values of the mask are scaled so that the
 sum of mask values = 1.0
*/
#define NR_END 1
#define FREE_ARG char*

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SQR(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define IMIN(a,b) ((a) < (b) ? (a) : (b))

/* Numerical Recipes standard error handler */
void nrerror(char error_text[])
{
 fprintf(stderr,"Numerical Recipes run-time error...\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"...now exiting to system...\n");
 exit(1);
}

/* allocate a float vector with subscript range v[nl..nh] */
float *vector(long nl,long nh)
{
 float *v;

 v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
 if (!v) nrerror("allocation failure in vector()");
 return v-nl+NR_END;
}

/* free a float vector allocated with vector() */
void free_vector(float* v,long nl,long nh)
{
 free((FREE_ARG) (v+nl-NR_END));
}

float pythag(float a, float b)
{
 float absa,absb;

 absa=fabs(a);
 absb=fabs(b);
 if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
 else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmp(float** a,int m,int n,float* w,float** v)
{
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}

void svbksb(float** u,float* w,float** v,int m,int n, float* b,float* x)
{
	int jj,j,i;
	float s,*tmp;

	tmp=vector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_vector(tmp,1,n);
}


void solve_system(int m, int n, float** a, float* x, float* b)

{
 int i,j;
 float w_min, w_max, **a_temp, *w, **v;

 w=(float *)malloc((n+1)*sizeof(float));

 v=(float **)malloc((n+1)*sizeof(float *));
 for(i=0; i<n+1; i++)
   v[i]=(float *)malloc((n+1)*sizeof(float));

 a_temp=(float **)malloc((m+1)*sizeof(float *));
 for(i=0; i<m+1; i++)
   a_temp[i]=(float *)malloc((n+1)*sizeof(float));

 /* copy a to a_temp because svdcmp changes a */

 for(i=0; i<m; i++)
   for(j=0; j<n; j++)
     a_temp[i+1][j+1]=a[i+1][j+1];

 svdcmp(a_temp,m,n,w,v);

 w_min=w_max=w[1];
 for(j=2; j<=n; j++) {
   if(w_min > w[j]) w_min=w[j];
   if(w_max < w[j]) w_max=w[j];
 }

 w_min=w_max*1.0e-6;
 for(j=1; j<=n; j++)
   if(w[j] < w_min) w[j]=0.0;

 svbksb(a_temp,w,v,m,n,b,x);

 for(i=0; i<m+1; i++)
   free(a_temp[i]);
 free(a_temp);

 free(w);

 for(i=0; i<n+1; i++)
   free(v[i]);
 free(v);
}

void WriteImage(char fname[], int **fimage, int M, int N, int Q)
{
 int i, j;
 unsigned char *image;
 std::ofstream ofp;

 image = (unsigned char *) new unsigned char [M*N];

 // convert the integer values to unsigned char

 for(i=0; i<N; i++)
   for(j=0; j<M; j++)
     image[i*M+j]=(unsigned char)fimage[i][j];

 ofp.open(fname, std::ios::out);

 if (!ofp) {
   std::cout << "Can't open file: " << fname << std::endl;
   exit(1);
 }

 ofp << "P5" << std::endl;
 ofp << M << " " << N << std::endl;
 ofp << Q << std::endl;

 ofp.write( reinterpret_cast<char *>(image), (M*N)*sizeof(unsigned char));

 if (ofp.fail()) {
   std::cout << "Can't write image " << fname << std::endl;
   exit(0);
 }

 ofp.close();

}

void ReadImage(char fname[], int ***fimage, int& M, int& N, int& Q)
{
 int i, j;
 unsigned char *image;
 char header [100], *ptr;
 std::ifstream ifp;

 ifp.open(fname, std::ios::in);

 if (!ifp) {
   std::cout << "Can't read image: " << fname << std::endl;
   exit(1);
 }

 // read header

 ifp.getline(header,100,'\n');
 if ( (header[0]!=80) ||    /* 'P' */
      (header[1]!=53) ) {   /* '5' */
      std::cout << "Image " << fname << " is not PGM" << std::endl;
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
   std::cout << "Image " << fname << " has wrong size" << std::endl;
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

void transform(float* a, float &cx, float &cy, int x, int y){
  cx = a[1]*x + a[2]*y + a[3];
  cy = a[4]*x + a[5]*y + a[6];
}

void invert(float k11, float k12, float k13,
            float k21, float k22, float k23,float* arr){
  arr[1] = -k22/((k12*k21)-(k11*k22));
  arr[2] = k12/((k12*k21)-(k11*k22));
  arr[3] = ((k13*k22)-(k12*k23))/((k12*k21)-(k11*k22));
  arr[4] = k21/((k12*k21)-(k11*k22));
  arr[5] = -k11/((k12*k21)-(k11*k22));
  arr[6] = -((k13*k21)-(k11*k23))/((k12*k21)-(k11*k22));
  //k11x1 + k12y1 + k13
  //k21x1 + k22y1 + k23
}

/*float convert(float OldMax, float OldMin, float NewMax, float NewMin, float OldValue){
  float OldRange = (OldMax - OldMin);
  float NewRange = (NewMax - NewMin);
//  std::cout<<(((OldValue - OldMin) * NewRange) / OldRange) + NewMin<<std::endl;
  return (((OldValue - OldMin) * NewRange) / OldRange) + NewMin;
}*/

int main(int argc, char** argv)
{
    int** input, **output;
    char name[20] = "20.pgm";
    char out[20] = "20_out.pgm";
    int x1, x2, x3, x4, y1, y2, y3, y4;
    //float px1 = 16, px2 = 16, px3 = 24, px4 = 32, py1 = 15, py2 = 25, py3 = 20, py4 = 20;
    float px1 = 20, px2 = 20, px3 = 30, px4 = 40, py1 = 9, py2 = 28, py3 = 19, py4 = 19;
    std::cout << "Input the coordinates for the left eye (Use this format: x y):";
    std::cin >> x1 >> y1;
    std::cout << "Input the coordinates for the right eye (Use this format: x y):";
    std::cin >> x2 >> y2;
    std::cout << "Input the coordinates for the nose (Use this format: x y):";
    std::cin >> x3 >> y3;
    std::cout << "Input the coordinates for the mouth center (Use this format: x y):";
    std::cin >> x4 >> y4;
    float P[9][7] ={ {0,0,0,0,0,0,0},
                     {0,x1, y1, 1,0,0,0},
                     {0,0,0,0,x1, y1, 1},
                     {0,x2, y2, 1,0,0,0},
                     {0,0,0,0,x2, y2, 1},
                     {0,x3, y3, 1,0,0,0},
                     {0,0,0,0,x3, y3, 1},
                     {0,x4, y4, 1,0,0,0},
                     {0,0,0,0,x4, y4, 1} };
    float* c = new float[3];
    float* cxy = new float[3];
    float pxy[9] = {0,px1,py1,px2,py2,px3,py3,px4,py4};
    float** Pt = new float*[9];
    float* pk = new float[7];
    float* ik = new float[7];
    float* ptxy = new float[8];
    for(int i = 0; i < 9; i++){
      Pt[i] = new float[7];
      for(int j = 1; j < 7; j++){
        Pt[i][j] = P[i][j];
      }
    }
    for(int i = 1; i < 9; i++){
      ptxy[i] = pxy[i];
    }
    solve_system(8,6,Pt,pk,ptxy);
    std::cout << pk[1] <<" " << pk[2] << " " << pk[3] << std::endl;
    std::cout << pk[4] <<" " << pk[5] << " " << pk[6] << std::endl;
    invert(pk[1],pk[2],pk[3],pk[4],pk[5],pk[6],ik);
    //std::cout << ik[1] <<" " << ik[2] << " " << ik[3] << std::endl;
    //std::cout << ik[4] <<" " << ik[5] << " " << ik[6] << std::endl;
    float xf, yf;
    int x_size, y_size, Q;
    ReadImage(name, &input, x_size, y_size, Q);
    ReadImage(name, &output, x_size, y_size, Q);
    float x, y;
    //std::cout << y_size << std::endl;
    for(int i = 0; i < 48; i++){
      for(int j = 0; j < 40;j++){
        //std::cout << i << " " << j << std::endl;
        transform(ik,xf,yf,i,j);
        //std::cout << xf << " " << yf << std::endl;
        x = abs(xf);
        y = abs(yf);
        output[i][j] = input[int(x)][int(y)];
        if((i == px1 && j == py1 ) || (i == px2 && j == py2 ) || (i == px3 && j == py3) || (i == px4 && j == py4 )){
          transform(pk,xf,yf,x,y);
          std::cout << "X pixel error for points "<< i << " " << j << ": " << abs(xf - px1) << " Y pixel error: " << abs(yf-py1) << std::endl;
        }

        xf = 0;
        yf = 0;
        //std::cout << "here" << std::endl;
      }
    }

    WriteImage(out, output, 40, 48, Q);
    return 0;
}
