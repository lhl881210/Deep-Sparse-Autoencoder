#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include   <float.h>
    
 
 
 
 
 int ** CreateTwoDimensionalArray ( unsigned int  uiYsize, unsigned int uiXsize)
{
      int  **ppiHead = NULL;  
      
      unsigned int  uiYindex = 0;   
    
      ppiHead = (int **)malloc(uiYsize*sizeof(int *));  

      for(uiYindex = 0; uiYindex < uiYsize;uiYindex++) 
     {
          *(ppiHead + uiYindex) = (int *) malloc(uiXsize*sizeof(int)); 
     }
     return  ppiHead;  
}


double ** d_CreateTwoDimensionalArray ( unsigned int  uiYsize,unsigned  int uiXsize)
{
      double  **ppiHead = NULL;   
      
      unsigned int  uiYindex = 0;   
    
      ppiHead = (double **)malloc(uiYsize*sizeof(double *));  

      for(uiYindex = 0; uiYindex < uiYsize;uiYindex++) 
     {
          *(ppiHead + uiYindex) = (double *) malloc(uiXsize*sizeof(double)); 
     }
     return  ppiHead;   
}

char ** ch_CreateTwoDimensionalArray ( unsigned int  uiYsize,unsigned  int uiXsize)
{
      char  **ppiHead = NULL;   
      
      unsigned int  uiYindex = 0;   
    
      ppiHead = (char **)malloc(uiYsize*sizeof(char *));  

      for(uiYindex = 0; uiYindex < uiYsize;uiYindex++) 
     {
          *(ppiHead + uiYindex) = (char *) malloc(uiXsize*sizeof(char)); 
     }
     return  ppiHead;   
}

void DeleteTwoDimensionalArray(double **ppiHead, int  uiYsize,  int uiXsize)
{
       int uiYindex;
    
      for(uiYindex = 0;uiYindex < uiYsize;uiYindex++)
     {
          free(*(ppiHead + uiYindex)); 
     }
	  free(ppiHead); 
}

 int *CreateONEDimensionalArray ( unsigned int  uisize){
        int  *piHead = NULL;
        piHead = (int *) malloc(uisize*sizeof(int));
      
        return  piHead;

}
 double *d_CreateONEDimensionalArray ( unsigned int  uisize){
        double  *piHead = NULL;
        piHead = (double *) malloc(uisize*sizeof(double));
      
        return  piHead;

}

 char *ch_CreateONEDimensionalArray ( unsigned int  uisize){
        char  *piHead;
        piHead = (char *) malloc(uisize*sizeof(char));
      
        return  piHead;

}

 void clear_vector(double *w,int size)
{
	int i;
#pragma omp parallel for private(i) schedule(dynamic)
	for(i=0;i<size;i++){
		w[i]=0;
	}
}
 void clear_vector_i(int *w,int size)
{
	int i;
#pragma omp parallel for private(i) schedule(dynamic)
	for(i=0;i<size;i++){
		w[i]=0;
	}
}

void clear_matrix(double **L,int sizeX,int sizeY)
{
	int i,j;
#pragma omp parallel for private(i,j) schedule(dynamic)
	for(i=0;i<sizeY;i++){
		for(j=0;j<sizeX;j++){
			L[i][j]=0;
		}
	}
}


void AVG_matrix(double **X,double *X_AVG,int XD,int N){
       int i=0,j=0;
       clear_vector(X_AVG,XD);
	   
       for(i=0;i<XD;i++)
                          
                      for(j=0;j<N;j++)
						X_AVG[i]=X_AVG[i]+X[j][i]/N;
       
       
       }


void transposed_matrix( double **a,double **b,int Y,int X){
       int i=0, j=0;
	//int temp;
       

//double **b = d_CreateTwoDimensionalArray(Y,X); 

for(j=0;j<Y;j++)  for(i=0;i<X;i++)  b[j][i]=a[i][j];
   //return b;            
      
       }


double **Autocorrelation_matrix( int **a, int **b,int d,int N){
       int i=0,j=0,k=0;
       double **c=NULL;
       c=d_CreateTwoDimensionalArray(N,N); 
       for(i=0;i<N;i++)
                          for(j=0;j<N;j++)
                                             for(k=0;k<d;k++) c[i][j]+=(double)(a[i][k]*b[k][j])/N;
                                                             
                                                                  
       return c;                           
       }                    
void MXM_multiplication( double **a, double **b,double **c,int D,int N){
       int i=0,j=0,k=0;
      // double **c=NULL;
      // c=d_CreateTwoDimensionalArray(N,D); 
       for(i=0;i<N;i++)
                          for(j=0;j<D;j++) c[i][j]=0.0;
       for(i=0;i<N;i++)
                          for(j=0;j<D;j++)
                                             for(k=0;k<N;k++) c[i][j]+=a[i][k]*b[k][j];
                                
       }                    
void MXM2(double **a,double **b,double **c,int D,int N){
	int i=0,j=0,k=0;
       //double **c=d_CreateTwoDimensionalArray(D,D); 
		clear_matrix(c,D,D);
		 for(i=0;i<D;i++)
                          for(j=0;j<D;j++)
                                             for(k=0;k<N;k++) c[i][j]+=a[k][i]*b[k][j];
                                                             
                                                                  
       //return c;                           

}
void MXM_PCA(double **a,double **b,double **c,int XD,int lowD,int N){
	int i=0,j=0,k=0;
       //double **c=d_CreateTwoDimensionalArray(D,D); 
		clear_matrix(c,lowD,N);
		 for(i=0;i<N;i++)
                          for(j=0;j<lowD;j++)
                                             for(k=0;k<XD;k++) c[i][j]+=a[i][k]*b[k][j];
                                                             
                                                                  
       //return c;                           

}
void MXMT(double **a,double **b,double **c,int D,int N){
	int i=0,j=0,k=0;
       //double **c=d_CreateTwoDimensionalArray(D,D); 
		clear_matrix(c,N,N);
		 for(i=0;i<N;i++)
                          for(j=0;j<N;j++)
                                             for(k=0;k<D;k++) c[i][j]+=a[i][k]*b[k][j];
                                                             
                                                                  
       //return c;                           

}
void WXV(double **W,double **V,double **H,int VD,int VN,int HD){
	int a=0,b=0,c=0;
    clear_matrix(H,HD,VN);
	for(a=0;a<VN;a++){
			for(b=0;b<HD;b++){
				for(c=0;c<VD;c++)
					H[a][b]+=W[b][c]*V[a][c];
				// printf("%f\t",H[a][b]);
			}
	// printf("\n");
	}
}
void MXV_multiplication( double **a, double *b,double *c,int N){
       int i=0,k=0;
       //double *c=d_CreateONEDimensionalArray(N); 
      // for(i=0;i<N;i++) printf("b=%f\t",b[i]);
      for(i=0;i<N;i++)
                           c[i]=0.0;
       for(i=0;i<N;i++)
                         
                          for(k=0;k<N;k++) {c[i]=c[i]+a[i][k]*b[k];
	   //printf("a=%f\t",a[i][k]);
	   }                 
                                                                  
       //return c;                           
       }
void MPM(double **a,double **b,double **c,int D,int N){
	int i=0,j=0,k=0;
       //double **c=d_CreateTwoDimensionalArray(D,N); 
	   for(i=0;i<D;i++)
                       for(k=0;k<N;k++) c[i][k]=a[i][k]+b[i][k];
      // return c;                           

}
void MmM(double **a,double **b,double **c,int D,int N){
	int i=0,j=0,k=0;
      // double **c=d_CreateTwoDimensionalArray(D,N); 
	   for(i=0;i<D;i++)
                       for(k=0;k<N;k++) c[i][k]=a[i][k]-b[i][k];
       //return c;                           

}

void INV_matrix(double **a,double **inv_a,int size){
      // double **inv_a=d_CreateTwoDimensionalArray(size,size);
       double temp;
       int i,j,k;
       
       for(i=0;i<size;i++){
                        for(j=0;j<size;j++){
                                         if(i==j)inv_a[i][j]=1;
                                         else inv_a[i][j]=0;
                                         }
                        }
       
       for(i=0;i<size;i++){
                        temp=1/a[i][i];
                        for(j=0;j<size;j++){
                                         a[i][j]*=temp;
                                         inv_a[i][j]*=temp;
                                         }
                        for(j=0;j<size;j++){
                                         if(i!=j){
                                                  temp=a[j][i];
                                                  for(k=0;k<size;k++){
                                                                   a[j][k]-=a[i][k]*temp;
                                                                   inv_a[j][k]-=inv_a[i][k]*temp;
                                                                   }
                                                  }
                                         }
                        }
      
      /* printf("INV-----------------\n");
      for(i=0;i<size;i++){
                        for(j=0;j<size;j++){
                                         printf(" %f",inv_a[i][j]);
                                         }
                        printf("\n");
                        }*/
      
      //return inv_a;
       }

 
       
       
       
       



void Cholesky(double **A,double **TL,int size)
{
	
	double temp;
	double **L=d_CreateTwoDimensionalArray(size,size);
	clear_matrix(L,size,size);

	int i,k,m;

	double *w=d_CreateONEDimensionalArray(size);
	//double **TL=d_CreateTwoDimensionalArray(size,size);
	clear_vector(w,size);
	for(i=0;i<size;i++){
		for(m=i;m<size;m++){
			w[m]=A[i][m];
		}

		for(k=0;k<i;k++){
			temp=L[k][i];
			if(temp!=0){
				for(m=i;m<size;m++){
					w[m] -= temp*L[k][m];
				}
			}
		}

		w[i]=sqrt(w[i]);
		for(m=i+1;m<size;m++){
			w[m] /=w[i];
		}

		for(m=i;m<size;m++){
			L[i][m]=w[m];
		}
		
		clear_vector(w,size);
	}
	 transposed_matrix(L,TL,size,size);
	free(w);
	DeleteTwoDimensionalArray(L,size,size);
	
	// return TL;
}





double **Normalization(double **X,double **Y,int n,int d){
       int i,j;
       double *max=d_CreateONEDimensionalArray(d);
       for(j=0;j<d;j++)max[j]=(-LDBL_MAX);
       double *min=d_CreateONEDimensionalArray(d);
       for(j=0;j<d;j++)min[j]=LDBL_MAX;
       for(j=0;j<d;j++)
                       for(i=0;i<n;i++){
                                        if(X[i][j]>max[j])max[j]=X[i][j];
                                        if(X[i][j]<min[j])min[j]=X[i][j];
                                                     
                                        }
       
         for(j=0;j<d;j++){
                          if(max[j]==min[j]) for(i=0;i<n;i++) Y[i][j]=0;
                    else   for(i=0;i<n;i++){
                                        Y[i][j]=(X[i][j]-min[j])/(max[j]-min[j]);
                                        Y[i][j]=((0.5-0.001)/0.5)*(Y[i][j])+0.001;
                                        Y[i][j]=Y[i][j]*2 -1;
                                        }
                       }
       free(max);
       free(min);
       return Y;
       }



double **Normalization2(double **X,double **Y,int n,int d){
       int i,j;
       double *mean=d_CreateONEDimensionalArray(d);
       //for(j=0;j<d;j++)max[j]=(-LDBL_MAX);
       double *sig=d_CreateONEDimensionalArray(d);
       //for(j=0;j<d;j++)min[j]=LDBL_MAX;
       for(j=0;j<d;j++)
                       for(i=0,mean[j]=0;i<n;i++)
						   mean[j]=mean[j]+X[i][j];

	   for(j=0;j<d;j++)
		   mean[j]=mean[j]/n;

	   for(j=0;j<d;j++)
                       for(i=0,sig[j]=0;i<n;i++)
						   sig[j]=sig[j]+(X[i][j]-mean[j])*(X[i][j]-mean[j]);

	   for(j=0;j<d;j++)
		   sig[j]=sqrt(sig[j]/n);


         for(j=0;j<d;j++){
                          if(sig[j]==0) for(i=0;i<n;i++) Y[i][j]=0;
                    else   for(i=0;i<n;i++){
                                        Y[i][j]=(X[i][j]-mean[j])/sig[j];
                                   
                                        }
                       }
       free(mean);
       free(sig);
       return Y;
       }

double **Normalization2_for_test(double **X,double **Y,double **T,int n,int d,int Tn){
       int i,j;
       double *mean=d_CreateONEDimensionalArray(d);
       //for(j=0;j<d;j++)max[j]=(-LDBL_MAX);
       double *sig=d_CreateONEDimensionalArray(d);
       //for(j=0;j<d;j++)min[j]=LDBL_MAX;
       for(j=0;j<d;j++)
                       for(i=0,mean[j]=0;i<n;i++)
						   mean[j]=mean[j]+X[i][j];

	   for(j=0;j<d;j++)
		   mean[j]=mean[j]/n;

	   for(j=0;j<d;j++)
                       for(i=0,sig[j]=0;i<n;i++)
						   sig[j]=sig[j]+(X[i][j]-mean[j])*(X[i][j]-mean[j]);

	   for(j=0;j<d;j++)
		   sig[j]=sqrt(sig[j]/n);


         for(j=0;j<d;j++){
                          if(sig[j]==0) for(i=0;i<n;i++) Y[i][j]=0;
                    else   for(i=0;i<Tn;i++){
                                        Y[i][j]=(T[i][j]-mean[j])/sig[j];
                                   
                                        }
                       }
       free(mean);
       free(sig);
       return Y;
       }





double **Normalization_kernel(double **X,double **Y,int n,int d){
       int i,j;
    
         for(j=0;j<d;j++){
                          
                      for(i=0;i<n;i++){
                                        Y[i][j]=(X[i][j])/(X[i][i]);
                                       // Y[i][j]=((0.5-0.001)/0.5)*(Y[i][j])+0.001;
                                        Y[i][j]=Y[i][j]*2 -1;
                                        }
                       }
       return Y;
       }      
       
double **Re_Normalization(double **X,double **Y,double **RE_X,int n,int d){
       int i,j;
       double *max=d_CreateONEDimensionalArray(d);
       for(j=0;j<d;j++)max[j]=(-LDBL_MAX);
       double *min=d_CreateONEDimensionalArray(d);
       for(j=0;j<d;j++)min[j]=LDBL_MAX;
       for(j=0;j<d;j++)
                       for(i=0;i<n;i++){
                                        if(X[i][j]>max[j])max[j]=X[i][j];
                                        if(X[i][j]<min[j])min[j]=X[i][j];
                                                     
                                        }
       for(j=0;j<d;j++)
                       for(i=0;i<n;i++){
                                        Y[i][j]=((Y[i][j]-0.001)*0.5)/(0.5-0.001);
                                        RE_X[i][j]=Y[i][j]*(max[j]-min[j])+min[j];
                                      }
       free(max);
       free(min);
       return RE_X;
       }
        
double *M_TO_V(double **M,double *V,int Y,int X){
       int i,j;
       for(i=0;i<Y;i++)
                       for(j=0;j<X;j++)
                                       V[i*X+j]=M[i][j];
       
       return V;
       
       }       
double **V_TO_M(double *V,double **M,int Y,int X){
       int i,j;
       for(i=0;i<Y;i++)
                       for(j=0;j<X;j++)
                                       M[i][j]=V[i*X+j];
       
       return M;
       }

/*Code is far away from bug with the animal protecting　
*　　　　　　　　┏┓　　　┏┓+ +
*　　　　　　　┏┛┻━━━┛┻┓ + +
*　　　　　　　┃　　　　　　　┃ 　
*　　　　　　　┃　　　━　　　┃ ++ + + +
*              ┃　┳┛　┗┳　┃
*　　　　　　　┃　　　　　　　┃ +
*　　　　　　　┃　　　┻　　　┃
*　　　　　　　┃　　　　　　　┃ + +
*　　　　　　　┗━┓　　　┏━┛
*　　　　　　　　　┃　　　┃　　　　　　　　　　　
*　　　　　　　　　┃　　　┃ + + + +
*　　　　　　　　　┃　　　┃　　　　　　　　　　
*　　　　　　　　　┃　　　┃ + 　　　　　　
*　　　　　　　　　┃　　　┃
*　　　　　　　　　┃　　　┃　　+　　　　　　　　　
*　　　　　　　　　┃　 　 ┗━━━┓ + +
*　　　　　　　　　┃ 　　　 　　　┣┓
*　　　　　　　　　┃ 　　　　 　　┏┛
*　　　　　　　　　┗┓┓┏━┳┓┏┛ + + + +
*　　　　　　　　　　┃┫┫　┃┫┫
*　　　　　　　　　　┗┻┛　┗┻┛+ + + +
*/