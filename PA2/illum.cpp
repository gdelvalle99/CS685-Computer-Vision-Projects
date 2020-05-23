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

float transform(float x, float y, float* eq){
  return x*eq[1] + y*eq[2] + x*y*eq[3] + eq[4];
}

int main(int argc, char** argv){
  int** input, **output;
  char name[20] = "1_out.pgm";
  char out[20] = "1_light.pgm";
  char outfile[20] = "1_norm.pgm";
  int x_size, y_size, Q;
  ReadImage(name, &input, x_size, y_size, Q);
  ReadImage(name, &output, x_size, y_size, Q);
  std::cout << x_size << " " << y_size << std::endl;
  float c[x_size * y_size];
  for (int q = 0; q < y_size; q++){
    for (int t = 0; t < x_size; t++){
        c[q * x_size + t] = input[q][t];
    }
  }
  //std::cout << "here";
  float** nm = new float*[x_size * y_size + 1];
  float** temp = new float*[x_size * y_size];
  for(int i = 0; i < y_size; i++){
    for(int j = 0; j < x_size; j++){
      temp[i*x_size+j] = new float[4];
      //std::cout<<i*y_size+j<<std::endl;
      temp[i*x_size+j][0] = j;
      //std::cout << "here" << std::endl;
      temp[i*x_size+j][1] = i;
      temp[i*x_size+j][2] = i*j;
      //std::cout << "here" << std::endl;
      temp[i*x_size+j][3] = 1;
      //std::cout << i << " " << j << std::endl;
    }


  }

  //std::cout << "here" << std::endl;
  nm[0] = new float[5];
  for(int i = 1; i < (x_size * y_size) + 1 ; i++){
    nm[i] = new float[5];
    for(int j = 1; j < 5; j++){
      nm[i][j] = temp[i-1][j-1];
    }
  }
  //std::cout << "here" << std::endl;

  float* b = new float[x_size * y_size + 1];
  for(int q = 1; q < x_size*y_size+1; q++){
    b[q] = c[q-1];

  }
  float* k = new float[5];
  //std::cout << "here" << std::endl;
  solve_system(x_size*y_size,4,nm,k,b);
  std::cout << k[1] << " " << k[2] << " " << k[3] << " " << k[4] << std::endl;
  ReadImage(name, &input, x_size, y_size, Q);
  ReadImage(name, &output, x_size, y_size, Q);
  float mean = 0;
  for(int i = 0; i < y_size; i++){
    for(int j = 0; j < x_size; j++){
      output[i][j] = int(abs(transform(i,j,k)));
      mean += output[i][j];
    }
  }
  mean /= (x_size * y_size);
  WriteImage(out, output, x_size, y_size, Q);
  ReadImage(name, &input, x_size, y_size, Q);
  for(int i = 0; i < y_size; i++){
    for(int j = 0; j < x_size; j++){
      input[i][j] = int(abs(input[i][j] - transform(i,j,k)+mean));
    }
  }
  WriteImage(outfile, input, x_size, y_size, Q);
}
