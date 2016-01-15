#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <windows.h>
#endif

 int rndn (int r){
  int rndno;
  while((rndno=((int)rand()/(RAND_MAX+1.0))*r)==r);
  return rndno;
}

double d_rndn (double r){
  double rndno;
   
  while((rndno=((double)rand() / (RAND_MAX + 1.0))*r)==r);
  return rndno;
}
 
   
  double  Guassian2(double   mu,double   sigma,int ppp)  
  {  
    double   r1,r2,y,out;  
   
    r1=d_rndn(0.99999999)+0.000000001;  
    r2=d_rndn(0.99999999)+0.000000001;  
   
#pragma omp parallel if (2)	
	if(ppp==1)y=sqrt(-2*log(r1))*cos(2*M_PI*r2);
	else if(ppp==0)y=sqrt(-2*log(r1))*sin(2*M_PI*r2);


   out=y*sigma+mu;
   //printf("Noies=%f\n",y*sigma);
    return   out;  
   
  }  
   

   
