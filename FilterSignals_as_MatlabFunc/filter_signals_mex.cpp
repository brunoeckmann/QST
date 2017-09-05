#include <matrix.h>
#include <mex.h>   
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include "boost/multi_array.hpp"

using boost::multi_array;
using boost::multi_array_ref;
using boost::const_multi_array_ref;
using boost::extents;
using boost::fortran_storage_order;

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
#define MWSIZE_MAX    281474976710655UL
#define MWINDEX_MAX   281474976710655UL
#define MWSINDEX_MAX  281474976710655L
#define MWSINDEX_MIN -281474976710655L
#else
#define MWSIZE_MAX    2147483647UL
#define MWINDEX_MAX   2147483647UL
#define MWSINDEX_MAX  2147483647L
#define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

#define MACH_EPS 1e-15
#define TINY 1.0e-20;
#define NR_END 1
#define FREE_ARG char*

const double ZERO=5e-3;

typedef boost::multi_array_types::extent_range range;

typedef boost::multi_array<double, 2> double2d; // for the MM matrix
typedef boost::multi_array<double, 3> double3d; // for signaling marginals
typedef boost::multi_array<double, 4> double4d; // for P(r,s,a,b)

void lubksb(double2d& a,int n,int *indx,double *b)
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

void ludcmp(double2d& a,int n,int *indx,double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  typedef boost::multi_array<double, 1> vector_type;
  typedef boost::multi_array_types::extent_range range;
  vector_type::extent_gen extents;
  vector_type vv(extents[range(1,n+1)]); // n+1 not included

  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0)
      throw std::runtime_error("Singular matrix in routine LUDCMP");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
}

void Remove_Signals(multi_array_ref<double, 4>& P_E_prev, multi_array_ref<double, 4>& P_E, int M, int W);
void Check_Signalling(double4d& P_E, int M, int W);

double Normalize(double4d& P_E_sign, int M, int W) {
  double deviation = 0.0;
  for (int n1=0;n1<M;n1++)
    for (int n2=0;n2<M;n2++) {
      double norm=0.;
      for (int w1=0;w1<W;w1++)
	for (int w2=0;w2<W;w2++)
	  norm+=P_E_sign[w1][w2][n1][n2];
      if (fabs(norm-1) > deviation)
	deviation = fabs(norm-1);
      for (int w1=0;w1<W;w1++)
	for (int w2=0;w2<W;w2++)
	  P_E_sign[w1][w2][n1][n2]/=norm;
    }
}

void LoadData(const char* filename, double4d& P, int M, int W) {
  FILE *pIn = fopen(filename, "r");
  if (pIn==NULL){
    printf("Sorry, I cannot open the file\n");
    exit(0);
  }
  double a;
  for (int n1=0;n1<M;n1++)
    for (int n2=0;n2<M;n2++)
      for (int w1=0;w1<W;w1++)
	for (int w2=0;w2<W;w2++) {
	  int ww=fscanf(pIn,"%lf",&a);
	  if (ww==EOF){
	    printf("Error, input file too short. Please, check variables M and W\n"); 
	    exit(0);
	  }
	  if (a<0) a=0;
	  P[w1][w2][n1][n2]=a;
	  P[w1][w2][n1][n2]=a;
	}
  int ww=fscanf(pIn,"%lf",&a);
  if (ww!=EOF){
    printf("Error, input file too long. Please, check variables M and W\n"); 
    exit(0);
  }
  fclose(pIn);
}

void Iterations(multi_array_ref<double, 4>& P_E_prev, multi_array_ref<double, 4>& P_E, int M, int W) {
  
  for (int ite=0;ite<10000;ite++) {
    Remove_Signals(P_E_prev, P_E, M, W);  // Questa funzione cambia P_E[]..[].
    
    double diff=0.,norm=0.;
    for (int n1=0;n1<M;n1++)
      for (int n2=0;n2<M;n2++){
	for (int w1=0;w1<W;w1++)
	  for (int w2=0;w2<W;w2++){
	    diff+=(P_E_prev[w1][w2][n1][n2]-P_E[w1][w2][n1][n2])
	      *(P_E_prev[w1][w2][n1][n2]-P_E[w1][w2][n1][n2]);
	    norm+=P_E[w1][w2][n1][n2]*P_E[w1][w2][n1][n2];
	  }
      }
    diff=sqrt(diff/norm);

    //	printf("%e\n",diff);
    if (diff<1e-14) break;

    for (int n1=0;n1<M;n1++)
      for (int n2=0;n2<M;n2++)
	for (int w1=0;w1<W;w1++)
	  for (int w2=0;w2<W;w2++)
	    P_E_prev[w1][w2][n1][n2]=P_E[w1][w2][n1][n2];

  }
  
}

void Remove_Signals(multi_array_ref<double, 4>& P_E_prev, multi_array_ref<double, 4>& P_E, int M, int W)
{
  double d[1],b[2*W*M*M+1];
  int indx[2*W*M*M+1];

  double4d P(boost::extents[W][W][M][M]);
  double3d P1(boost::extents[W][M][M]);
  double3d P2(boost::extents[W][M][M]);
  double2d::extent_gen double2d_extents;
  double2d MM(double2d_extents[range(1,2*W*M*M+1)][range(1,2*W*M*M+1)]);

  for (int m2=0;m2<M;m2++)
    for (int m1=0;m1<M;m1++){
      for (int w2=0;w2<W;w2++){
	P1[w2][m1][m2]=0.;
	P2[w2][m1][m2]=0.;
      }
      for (int w1=0;w1<W;w1++){
	for (int w2=0;w2<W;w2++){
	  P[w1][w2][m1][m2]=P_E_prev[w1][w2][m1][m2];
	  P1[w1][m1][m2]+=P[w1][w2][m1][m2];
	  P2[w2][m1][m2]+=P[w1][w2][m1][m2];
	}
      }
    }

  for (int k1=0;k1<2*W*M*M;k1++)
    for (int k2=0;k2<2*W*M*M;k2++)
      MM[k1+1][k2+1]=0.;

  for (int k1=0;k1<W*M*M;k1++){
    /**** definition of the matrix ****/
    int m2_1=k1%M;
    int m1_1=(k1/M)%M;
    int w_1=k1/M/M;

    MM[k1+1][k1+1]=P1[w_1][m1_1][m2_1]+ZERO;
    MM[k1+1+W*M*M][k1+1+W*M*M]=P2[w_1][m1_1][m2_1]+ZERO;


    for (int w=0;w<W;w++){
      int k2=(w*M+m1_1)*M+m2_1;
      MM[k1+1][k2+1+W*M*M]+=P[w_1][w][m1_1][m2_1];
      MM[k1+1+W*M*M][k2+1]+=P[w][w_1][m1_1][m2_1];
    }
    for (int w=0;w<W;w++)
      for (int m=0;m<M;m++){
	int k2=(w*M+m1_1)*M+m;
	MM[k1+1][k2+1+W*M*M]-=P[w_1][w][m1_1][m]/M;
	int k2_2=(w*M+m)*M+m2_1;
	MM[k1+1+W*M*M][k2_2+1]-=P[w][w_1][m][m2_1]/M;
      }
    /**********************************/
    b[k1+1]      =-P1[w_1][m1_1][m2_1];
    b[k1+1+W*M*M]=-P2[w_1][m1_1][m2_1];
    for (int m=0;m<M;m++){
      b[k1+1]      +=P1[w_1][m1_1][m]/M;
      b[k1+1+W*M*M]+=P2[w_1][m][m2_1]/M;
    }
  }

  ludcmp(MM,2*W*M*M,indx,d);
  lubksb(MM,2*W*M*M,indx,b);

  double3d eta1(boost::extents[W][M][M]);
  double3d eta2(boost::extents[W][M][M]);

  for (int m2=0;m2<M;m2++)
    for (int m1=0;m1<M;m1++)
      for (int w1=0;w1<W;w1++){
	int k1=(w1*M+m1)*M+m2;
	eta1[w1][m1][m2]=b[k1+1];
	eta2[w1][m1][m2]=b[k1+1+W*M*M];
#if 0
	for (int m=0;m<M;m++){
	  int k2=(w1*M+m1)*M+m;
	  int k3=(w1*M+m)*M+m2;
	  eta1[w1][m1][m2]-=b[k2+1]/M;
	  eta2[w1][m1][m2]-=b[k3+1]/M;
	}
#endif
	//		printf("%e  %e\n",eta1[w1][m1][m2],b[k1+1]);
	//		printf("%e  %e\n",eta2[w1][m1][m2],b[k1+1+W*M*M]);
      }


	
  for (int m2=0;m2<M;m2++)
    for (int m1=0;m1<M;m1++){
      double norm=0.;
      for (int w1=0;w1<W;w1++)
	for (int w2=0;w2<W;w2++){
	  int k1=(w1*M+m1)*M+m2;
	  int k2=(w2*M+m1)*M+m2;
	  P_E[w1][w2][m1][m2]=P_E_prev[w1][w2][m1][m2]*exp(b[k1+1]+b[k2+1+W*M*M]);
	  norm+=P_E[w1][w2][m1][m2];
	}
      for (int w1=0;w1<W;w1++)
	for (int w2=0;w2<W;w2++)
	  P_E[w1][w2][m1][m2]/=norm;
    }


  return;
}

void Check_Signalling(double4d& P_E, int M, int W)
{
  /****** Checking signalling ******/
  double3d P_1(boost::extents[W][M][M]);
  double3d P_2(boost::extents[W][M][M]);
  double2d P_1_mean(boost::extents[W][M]);
  double2d P_2_mean(boost::extents[W][M]);

  double diff1=0.;
  double diff2=0.;
  for (int n1=0;n1<M;n1++){
    for (int n2=0;n2<M;n2++){
      for (int w1=0;w1<W;w1++){
	P_1[w1][n1][n2]=0.;
	P_2[w1][n1][n2]=0.;
	for (int w2=0;w2<W;w2++){
	  P_1[w1][n1][n2]+=P_E[w1][w2][n1][n2];
	  P_2[w1][n1][n2]+=P_E[w2][w1][n1][n2];
	}
      }
    }
  }
  for (int n1=0;n1<M;n1++){
    for (int w=0;w<W;w++){
      P_1_mean[w][n1]=0.;
      P_2_mean[w][n1]=0.;
      for (int n2=0;n2<M;n2++){
	P_1_mean[w][n1]+=P_1[w][n1][n2]/M;
	P_2_mean[w][n1]+=P_2[w][n2][n1]/M;
      }
    }
  }
  for (int n1=0;n1<M;n1++)
    for (int n2=0;n2<M;n2++)
      for (int w=0;w<W;w++){
	double d1=fabs(P_1[w][n1][n2]-P_1_mean[w][n1]);
	double d2=fabs(P_2[w][n1][n2]-P_2_mean[w][n2]);
	if (d1>diff1) diff1=d1;
	if (d2>diff2) diff2=d2;
      }
  printf("%e  %e\n",diff1,diff2);

  /*********************************/


  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // 3 input parameters: P_E, W, M
  mxAssert(nrhs == 3, "Wrong number of input parameters");
  // 1 output parameters: P_E_reg
  mxAssert(nlhs == 1, "Wrong number of output parameters");

  const mxArray* m_P = prhs[0];
  const mxArray* m_W = prhs[1];
  const mxArray* m_M = prhs[2];

  mxAssert(mxIsScalar(m_M) && mxIsDouble(m_M) && !mxIsComplex(m_M),
	   "Wrong input parameter M");
  
  int M = (int)mxGetScalar(m_M);
  
  mxAssert(mxIsScalar(m_W) && mxIsDouble(m_W) && !mxIsComplex(m_W),
	   "Wrong input parameter W");
  
  int W = (int)mxGetScalar(m_W);

  mxAssert(mxGetNumberOfDimensions(m_P) == 4 &&
	   mxGetDimensions(m_P)[0] == W &&
	   mxGetDimensions(m_P)[1] == W &&
	   mxGetDimensions(m_P)[2] == M &&
	   mxGetDimensions(m_P)[3] == M &&
	   !mxIsSparse(m_P) && mxIsDouble(m_P) && !mxIsComplex(m_P),
	   "Wrong input 4D array P");

  mxArray* m_P_E = mxDuplicateArray(m_P);
  mxArray* m_P_E_prev = mxDuplicateArray(m_P);

  multi_array_ref<double, 4> P_E(mxGetPr(m_P_E),
				 extents[W][W][M][M],
				 boost::fortran_storage_order());

  multi_array_ref<double, 4> P_E_prev(mxGetPr(m_P_E_prev),
				      extents[W][W][M][M],
				      boost::fortran_storage_order());

  Iterations(P_E_prev, P_E, M, W);

  plhs[0] = m_P_E;
  
  mxDestroyArray(m_P_E_prev);

}
