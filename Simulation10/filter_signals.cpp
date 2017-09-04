#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>
#include <stdlib.h>

#define MACH_EPS 1e-15

const double ZERO=5e-3;

const int N_expe=1;

const int M=3;  // Number of measurements
const int W=2;  // Number of measurement outcomes

double **MM;

double P_E[N_expe][W][W][M][M];
double P_E_prev[N_expe][W][W][M][M];  // signaling distribution
double P_E_sign[N_expe][W][W][M][M];  // signaling distribution


void Remove_Signals(void);
void Check_Signalling(double P_E[][W][W][M][M]);

int main(int argc,char **argv)
{
	if (argc<2) {
		printf("Please, provide an input file name\n");
		return 0;
	}

	FILE *pIn=fopen(argv[1],"r");
	if (pIn==NULL){
		printf("Sorry, I cannot open the file\n");
		return 0;
	}

	MM=new double *[2*W*M*M+1];
	for (int n=0;n<2*W*M*M+1;n++)
		MM[n]=new double [2*W*M*M+1];

	/**** loading the data ****/
	double a;

	for (int n1=0;n1<M;n1++)
		for (int n2=0;n2<M;n2++)
			for (int w1=0;w1<W;w1++)
				for (int w2=0;w2<W;w2++){
					for (int expe=0;expe<N_expe;expe++){
						int ww=fscanf(pIn,"%lf",&a);
						if (ww==EOF){
							printf("Error, input file too short. Please, check variables N_expe, M and W\n"); 
							exit(0);
						}
						if (a<0) a=0;
						P_E_sign[expe][w1][w2][n1][n2]=a;
						P_E_prev[expe][w1][w2][n1][n2]=a;
				//		printf("%e  ",a);
					}
				//	printf("\n");
				}
	int ww=fscanf(pIn,"%lf",&a);
	if (ww!=EOF){
		printf("Error, input file too long. Please, check variables N_expe, M and W\n"); 
		exit(0);
	}


//	exit(0);
	/**************************/

	/**** Check normalization ****/
	for (int expe=0;expe<N_expe;expe++)
		for (int n1=0;n1<M;n1++)
			for (int n2=0;n2<M;n2++){
				double norm=0.;
				for (int w1=0;w1<W;w1++)
					for (int w2=0;w2<W;w2++)
						norm+=P_E_sign[expe][w1][w2][n1][n2];
#if 0
				if (fabs(norm-1)>1e-6){
					printf("error, the distributions are not normalized\n");
					printf("norm=%e\n",norm);
					exit(0);
				}
#endif
				for (int w1=0;w1<W;w1++)
					for (int w2=0;w2<W;w2++)
						P_E_sign[expe][w1][w2][n1][n2]/=norm;
		}
	/*****************************/

	Check_Signalling(P_E_sign);

	for (int ite=0;ite<10000;ite++){
		Remove_Signals();  // Questa funzione cambia P_E[]..[].


		double diff=0.,norm=0.;
		for (int expe=0;expe<N_expe;expe++){
			for (int n1=0;n1<2;n1++)
				for (int n2=0;n2<2;n2++){
					for (int w1=0;w1<W;w1++)
						for (int w2=0;w2<W;w2++){
							diff+=(P_E_prev[expe][w1][w2][n1][n2]-P_E[expe][w1][w2][n1][n2])
							*(P_E_prev[expe][w1][w2][n1][n2]-P_E[expe][w1][w2][n1][n2]);
							norm+=P_E[expe][w1][w2][n1][n2]*P_E[expe][w1][w2][n1][n2];
						}
				}
		}
		diff=sqrt(diff/norm);

	//	printf("%e\n",diff);
		if (diff<1e-14) break;

		for (int n1=0;n1<M;n1++)
			for (int n2=0;n2<M;n2++)
				for (int w1=0;w1<W;w1++)
					for (int w2=0;w2<W;w2++)
						for (int expe=0;expe<N_expe;expe++){
							P_E_prev[expe][w1][w2][n1][n2]=P_E[expe][w1][w2][n1][n2];
						}

	}
	for (int n1=0;n1<M;n1++)
		for (int n2=0;n2<M;n2++)
			for (int w1=0;w1<W;w1++)
				for (int w2=0;w2<W;w2++){
					for (int expe=0;expe<N_expe;expe++)
						printf("%e ",P_E[expe][w1][w2][n1][n2]);
					printf("\n");
				}

	Check_Signalling(P_E);


	return 0;
}

void ludcmp(double **a,int n,int *indx,double *d);
void lubksb(double **a,int n,int *indx,double *b);

void Remove_Signals(void)
{
	double d[1],b[2*W*M*M+1];
	int indx[2*W*M*M+1];

	for (int expe=0;expe<N_expe;expe++){ // START EXPE LOOP

	double P[W][W][M][M];
	double P1[W][M][M],P2[W][M][M];
	for (int m2=0;m2<M;m2++)
		for (int m1=0;m1<M;m1++){
			for (int w2=0;w2<W;w2++){
				P1[w2][m1][m2]=0.;
				P2[w2][m1][m2]=0.;
			}
			for (int w1=0;w1<W;w1++){
				for (int w2=0;w2<W;w2++){
					P[w1][w2][m1][m2]=P_E_prev[expe][w1][w2][m1][m2];
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
    
    
// Perform linear program
    // First decompose the matrix M into LU decomposition (as input in lubksb routine
    ludcmp(MM,2*W*M*M,indx,d);
    // Then perform linear program
    lubksb(MM,2*W*M*M,indx,b);
    // note that b[1..n] is the input as the right hand side vector B, and returns with the solution vector x.

//  Some kind of test Alberto induced for testing the lin Prog
#if 0
	double eta1[W][M][M],eta2[W][M][M];

	for (int m2=0;m2<M;m2++)
		for (int m1=0;m1<M;m1++)
			for (int w1=0;w1<W;w1++){
				int k1=(w1*M+m1)*M+m2;
				eta1[w1][m1][m2]=b[k1+1];
				eta2[w1][m1][m2]=b[k1+1+W*M*M];

				for (int m=0;m<M;m++){
					int k2=(w1*M+m1)*M+m;
					int k3=(w1*M+m)*M+m2;
					eta1[w1][m1][m2]-=b[k2+1]/M;
					eta2[w1][m1][m2]-=b[k3+1]/M;
				}

		//		printf("%e  %e\n",eta1[w1][m1][m2],b[k1+1]);
		//		printf("%e  %e\n",eta2[w1][m1][m2],b[k1+1+W*M*M]);
			}
#endif

	for (int m2=0;m2<M;m2++)
		for (int m1=0;m1<M;m1++){
			double norm=0.;
			for (int w1=0;w1<W;w1++)
				for (int w2=0;w2<W;w2++){
					int k1=(w1*M+m1)*M+m2;
					int k2=(w2*M+m1)*M+m2;
					P_E[expe][w1][w2][m1][m2]=P_E_prev[expe][w1][w2][m1][m2]*exp(b[k1+1]+b[k2+1+W*M*M]);
					norm+=P_E[expe][w1][w2][m1][m2];
				}
			for (int w1=0;w1<W;w1++)
				for (int w2=0;w2<W;w2++)
					P_E[expe][w1][w2][m1][m2]/=norm;
		}


	}      				    // END EXPE LOOP

	return;
}

void Check_Signalling(double P_E[][W][W][M][M])
{
	/****** Checking signalling ******/
	double P_1[W][M][M],P_2[W][M][M];
	double P_1_mean[W][M],P_2_mean[W][M];

	double diff1=0.;
	double diff2=0.;
	for (int expe=0;expe<N_expe;expe++){
		for (int n1=0;n1<M;n1++){
			for (int n2=0;n2<M;n2++){
				for (int w1=0;w1<W;w1++){
					P_1[w1][n1][n2]=0.;
					P_2[w1][n1][n2]=0.;
					for (int w2=0;w2<W;w2++){
						P_1[w1][n1][n2]+=P_E[expe][w1][w2][n1][n2];
						P_2[w1][n1][n2]+=P_E[expe][w2][w1][n1][n2];
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
	}
	printf("%e  %e\n",diff1,diff2);

	/*********************************/


	return;
}


