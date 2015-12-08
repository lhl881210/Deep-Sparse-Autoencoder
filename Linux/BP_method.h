
#include <algorithm>
#include <time.h>
#include "my_conio.h"

double LR_h(double **V,double **V_denoising,double **Wh,double **Wr,double **dH_Wh,int VN,int VD,int HD,double eck,double p,double A,double B,double step){
       int i,j,k;
       double a=0,check=LDBL_MAX;
      double mintemp=LDBL_MAX,mina=0;
      double e_all,SP=0,L2_Wh;
      double temp=0,temp1=0;
       double *h=d_CreateONEDimensionalArray(HD+1);
       double *r=d_CreateONEDimensionalArray(VD+1);
       double *p_hat=d_CreateONEDimensionalArray(HD);
double **t_Wh=d_CreateTwoDimensionalArray(VD+1,HD);

      double SSS=step*0.1;
      while(check>=SSS){

                                      a=a+step;

									  temp1=temp;
									  clear_vector(p_hat,HD);
                                      for (j = 0; j < VD+1; j++)
                                          for (k = 0; k < HD; k++)
								              t_Wh[j][k]=Wh[j][k]-a*dH_Wh[j][k];


                                      for(i=0,e_all=0;i<VN;i++){
                                       h=hidden(V_denoising[i],h,t_Wh,VD,HD);
                                         for(j=0;j<HD;j++)
                                                            p_hat[j]=p_hat[j]+h[j];

                                              r=reconstruction(h,r,Wr,VD,HD);

                                              e_all=e_all+error(V[i],r,VD);
                                              }
                                              for(j=0;j<HD;j++)
                                                            p_hat[j]=(p_hat[j]/VN+1)/2;
                                               SP=sparce(p,p_hat,HD);
                                                for(L2_Wh=0,j=0;j<VD;j++)
                                                               for(k=0;k<HD;k++)
                                                                             L2_Wh=L2_Wh+t_Wh[j][k]*t_Wh[j][k];
                                               temp=e_all/VN+A*0.5*(L2_Wh)+B*SP;

                                     if(temp<mintemp){
                                                      mintemp=temp;

                                                    mina=a;
                                                      }


                                      check=fabs(temp1-temp);


									   if(temp>temp1) step=(-0.5)*step;





                                      }
       DeleteTwoDimensionalArray(t_Wh,VD+1,HD);
       free(h);
       free(r);
       free(p_hat);


       return mina;
}

double LR_r(double **V,double **V_denoising,double **Wh,double **Wr,double **dH_Wr,int VN,int VD,int HD,double eck,double p,double A,double B,double step){
       int i,j,k;
       double a=0,check=LDBL_MAX,L2_Wr=0;
      double mintemp=LDBL_MAX,mina=0;
double e_all;
      double temp=0,temp1=0;
       double **H=d_CreateTwoDimensionalArray(VN,HD+1);
       double *r=d_CreateONEDimensionalArray(VD+1);

       double **t_Wr=d_CreateTwoDimensionalArray(HD+1,VD);

      double SSS=step*0.1;
      for(i=0;i<VN;i++){
                                       H[i]=hidden(V_denoising[i],H[i],Wh,VD,HD);

                                                            }
      while(check>=SSS){

                                      a=a+step;

									  temp1=temp;

                                      for (j = 0; j < HD+1; j++)
                                          for (k = 0; k < VD; k++)
								              t_Wr[j][k] =Wr[j][k]-a*dH_Wr[j][k];


                                       for(i=0,e_all=0;i<VN;i++){
                                              r=reconstruction(H[i],r,t_Wr,VD,HD);

                                              e_all=e_all+error(V[i],r,VD);
                                              }
                                       for(L2_Wr=0,j=0;j<HD;j++)
                                                               for(k=0;k<VD;k++)
                                                                             L2_Wr=L2_Wr+t_Wr[j][k]*t_Wr[j][k];

                                               temp=e_all/VN+A*0.5*L2_Wr;

                                     if(temp<mintemp){
                                                      mintemp=temp;

                                                      mina=a;
                                                      }


                                      check=fabs(temp1-temp);


									  if(temp>temp1) step=(-0.5)*step;




                                      }
        DeleteTwoDimensionalArray(t_Wr,HD+1,VD);
       DeleteTwoDimensionalArray(H,VN,HD+1);
       free(r);


       return mina;
}

////////////////////////////////////BP():1???P??////////////////////////////////////////////////////////////////
void BP(double **V,double **V_denoising,double **Wh,double **Wr,int VN,int VD,int HD,double eck,double p,double A,double B,double step){
int i,j,k;
double **d_Wh=d_CreateTwoDimensionalArray(VD+1,HD);
double **d_Wr=d_CreateTwoDimensionalArray(HD+1,VD);

//double *h=d_CreateONEDimensionalArray(HD+1);
double **H = d_CreateTwoDimensionalArray(VN, HD + 1);
double *r=d_CreateONEDimensionalArray(VD+1);
double *p_hat=d_CreateONEDimensionalArray(HD);
double *d_r=d_CreateONEDimensionalArray(VD);
double *d_h=d_CreateONEDimensionalArray(HD);
double temp;


clear_matrix(d_Wr,VD,HD+1);
clear_matrix(d_Wh,HD,VD+1);
clear_vector(d_r,VD);
clear_vector(d_h,HD);
clear_vector(p_hat,HD);

for(i=0;i<VN;i++){
			H[i] = hidden(V_denoising[i], H[i], Wh, VD, HD);
                                 for(j=0;j<HD;j++)
									 p_hat[j] = p_hat[j] + H[i][j];

                                                      }

                            for(j=0;j<HD;j++){
                                        p_hat[j]=(p_hat[j]/VN+1)/2;

                                        }

                            for(i=0;i<VN;i++){

                                              //h=hidden(V_denoising[i],h,Wh,VD,HD);


											r = reconstruction(H[i], r, Wr, VD, HD);

                                             // e_all=e_all+error(V[i],r,VD);



                                             for (j = 0; j < VD; j++)
                                                  d_r[j] =(-(1 - r[j]*r[j]) * (V[i][j] - r[j]));

                                              for (j = 0; j < HD; j++)
                                                  for (k = 0; k < VD; k++)
													  d_Wr[j][k] += H[i][j] * d_r[k];
                                              for (k = 0; k < VD; k++)
                                                     d_Wr[HD][k] += d_r[k];





                                             for (j = 0; j < HD; j++) {
                                                   for (k = 0,temp=0; k < VD; k++)
                                                           temp = temp + Wr[j][k] * d_r[k];
												   d_h[j] = (temp + B*(((1 - p) / (1 - p_hat[j])) - (p / p_hat[j]))) * (1 - H[i][j] * H[i][j]);
												  }



                                             for (j = 0; j < VD; j++)
                                                  for (k = 0; k < HD; k++)
                                                       d_Wh[j][k] +=V[i][j] * d_h[k];
                                             for (k = 0; k < HD; k++)
                                                             d_Wh[VD][k]+= d_h[k];


                                             }



                   for (j = 0; j < HD; j++)
                               for (k = 0; k < VD; k++){
                                   d_Wr[j][k] = d_Wr[j][k]/VN+A*Wr[j][k];

                                   }
							    for (k = 0; k < VD; k++){
                                   d_Wr[HD][k] = d_Wr[HD][k]/VN;

                                   }
                   for (j = 0; j < VD; j++)
                               for (k = 0; k < HD; k++){
								   d_Wh[j][k] =  d_Wh[j][k]/VN+A*Wh[j][k];

                                   }
							   for (k = 0; k < HD; k++){
							    d_Wh[VD][k]= d_Wh[VD][k]/VN;

                               }




//printf("-------aa_r---------\n");
double aa_r=LR_r(V,V_denoising,Wh,Wr,d_Wr,VN,VD,HD,eck,p,A,B,step);
//printf("aa_r=%f\n",aa_r);
//printf("-------aa_h---------\n");
double aa_h=LR_h(V,V_denoising,Wh,Wr,d_Wh,VN,VD,HD,eck,p,A,B,step);
//printf("aa_h=%f\n",aa_h);


for (j = 0; j < HD+1; j++)
                               for (k = 0; k < VD; k++){
                                   Wr[j][k] = Wr[j][k] - aa_r*d_Wr[j][k];

                                   }
 for (j = 0; j < VD+1; j++)
                               for (k = 0; k < HD; k++)  {
								   Wh[j][k] = Wh[j][k] - aa_h*d_Wh[j][k];
                                  }

DeleteTwoDimensionalArray(d_Wh,VD+1,HD);
DeleteTwoDimensionalArray(d_Wr,HD+1,VD);
DeleteTwoDimensionalArray(H, VN, HD + 1);
//free(h);
free(r);
free(p_hat);
free(d_r);
free(d_h);
}






////////////////////////////////////BP_ALL():???????P??////////////////////////////////////////////////////////////////
void BP_ALL(double **V,double **Wh,double **Wr,int VN,int VD,int HD,double eck,double p,double A,double B){
int i,j,k;
double errorck=LDBL_MAX;

double SP=0,L2_Wh=0,L2_Wr=0,e_all=1;
double *p_hat=d_CreateONEDimensionalArray(HD+1);
double *r=d_CreateONEDimensionalArray(VD+1);
 double *h=d_CreateONEDimensionalArray(HD+1);

clock_t start, finish;
   double  duration;


double ertemp=0;
 double step=1.0;
double Var=1.0;
min_Var=setting[5];
 int loop=0;
 double **V_denoising=d_CreateTwoDimensionalArray(VN,VD+1);
while(errorck>= eck){
start = clock();
ertemp=e_all;
//?p?????[?^//
p=setting[0];//Sparse?̖ڕW?l
A=setting[1];//L2?m?????̋???
B=setting[2];//Sparse?̋???
eck=setting[3];//???????
denoising_rate=setting[4];//noise?̃p?Z???g


////////////////////////////Denoising/////////////////////////////////////////////
for(i=0;i<VN;i++){
	int ppp=1;
	for(j=0;j<VD;j++){
		// srand( (unsigned)time( NULL ) );
		double tmp = d_rndn(1);//printf("tmp=%f\t",tmp);
					if (tmp < denoising_rate&&tmp>0) {
						//printf("[%d][%d]=",i,j);
						V_denoising[i][j]=Gaussian2(V[i][j],Var,ppp);
						if(ppp==1)ppp=0;
						else if(ppp==0)ppp=1;
					}
					else V_denoising[i][j]=V[i][j];
					//printf("[%d][%d]=%f\n",i,j,V_denoising[i][j]);
	 }
	 V_denoising[i][VD]=1.0;
}

///////////////////////////////////////////////////////////////////////////////



BP(V,V_denoising,Wh,Wr,VN,VD,HD,eck,p,A,B,step);//BP():1???P??

 clear_vector(p_hat,HD);
for(i=0,e_all=0;i<VN;i++){

                                              h=hidden(V[i],h,Wh,VD,HD);

                                              for(j=0;j<HD;j++)
                                                            p_hat[j]=p_hat[j]+h[j];

                                              r=reconstruction(h,r,Wr,VD,HD);

                                              e_all=e_all+error(V[i],r,VD);

}
for(j=0;j<HD;j++){
                                        p_hat[j]=(p_hat[j]/VN+1)/2;
                                       if(p_hat[j]==0.0)p_hat[j]=0.0000000001;
                                        }
for(L2_Wh=0,j=0;j<VD;j++)
                                                               for(k=0;k<HD;k++)
                                                                             L2_Wh=L2_Wh+Wh[j][k]*Wh[j][k];
                                     for(L2_Wr=0,j=0;j<HD;j++)
                                                                for(k=0;k<VD;k++)
                                                                              L2_Wr=L2_Wr+Wr[j][k]*Wr[j][k];


FILE *BFGS;
BFGS=fopen(Gradient_filename,"a+");
e_all=e_all/VN;
SP=sparce(p,p_hat,HD);
printf("[%d]\te_all=%f\t",loop++,e_all);
fprintf(BFGS,"%f\t",e_all);
Var=noise_variance(e_all);
Var=std::max(min_Var,Var);

e_all=e_all+A*0.5*(L2_Wh+L2_Wr)+B*SP;

printf("E_all=%f\t",e_all);
fprintf(BFGS,"%f\t",e_all);

printf("Var=%f\t",Var);
fprintf(BFGS,"%f\t",Var);

errorck=fabs(ertemp-e_all);
step=errorck+0.1;
finish = clock();
   duration = (double)(finish - start) / CLOCKS_PER_SEC;
   printf( "%f seconds\n", duration );
   fprintf(BFGS,"%f\n",duration);
   fclose(BFGS);

////////////////////////?P?????ꂽ?B??(????)?A?e?d?ݍs???o?C?A?X??o??///////////////////////
if(_kbhit()){
FILE *Setting_file;
if((Setting_file=fopen("setting.txt","r")) == NULL){
        printf("Cannot open file setting.txt.\n");
		//system("pause");
        //exit(1);
    }

else{
  fscanf(Setting_file,"%*[^\n]%*c");
 setting=READ_FILE1(Setting_file,setting,6);
 fclose(Setting_file);
 //?p?????[?^//
p=setting[0];//Sparse?̖ڕW?l
A=setting[1];//L2?m?????̋???
B=setting[2];//Sparse?̋???
eck=setting[3];//???????
denoising_rate=setting[4];//noise?̃p?Z???g
min_Var=setting[5];
}

printf("//////////////Parameters//////////////\n");
printf("Sparse taget: p=%f\n",p);
printf("L2norm intensity: A=%f\n",A);
printf("Sparse intensity: B=%f\n",B);
printf("Convergence condition: eck=%f\n",eck);
printf("Noise rate: noise_rate=%f(%%)\n",denoising_rate);
printf("Variance of Gaussian noise=%f\n",Var);
printf("Min_Var=%f(%%)\n",min_Var=setting[5]);
printf("//////////////////////////////////////\n");
/////////////
////////////////?B??(????)??o??////////////////////
FILE *H_OUT;



                            if((H_OUT=fopen(H_OUT_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",H_OUT_filename);
							//system("pause");

						   }
							else{
                             for(i=0;i<VN;i++){
                                              h=hidden(V[i],h,Wh,VD,HD);
                                               for(k=0;k<HD;k++){
                                                               if(k<HD-1) fprintf(H_OUT,"%f\t",h[k]);
															   else fprintf(H_OUT,"%f",h[k]);
											   }
                                               fprintf(H_OUT,"\n");
                                          }

                              fclose(H_OUT);
							}
////////////////encoder?̏d?ݍs???o?C?A?X??o??////////////////////
FILE *W_H;


                            if((W_H=fopen(W_H_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_H_filename);
							//system("pause");

						   }
							else{
                             for(i=0;i<VD+1;i++){

                                               for(k=0;k<HD;k++){
                                                                if(k<HD-1)fprintf(W_H,"%f\t",Wh[i][k]);
																else fprintf(W_H,"%f",Wh[i][k]);
											   }
                                               fprintf(W_H,"\n");
                                          }

                              fclose(W_H);
							}
////////////////decoder?̏d?ݍs???o?C?A?X??o??////////////////////
FILE *W_R;


                           if((W_R=fopen(W_R_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_R_filename);
							//system("pause");

						   }
						   else{
                             for(i=0;i<HD+1;i++){

                                               for(k=0;k<VD;k++){
                                                                if(k<(VD-1))fprintf(W_R,"%f\t",Wr[i][k]);
																else fprintf(W_R,"%f",Wr[i][k]);
											   }
                                               fprintf(W_R,"\n");
                                          }

                              fclose(W_R);
						   }


////////////////Reconstruction??o??////////////////////
FILE *V_Re;

 if((V_Re=fopen(V_Re_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_R_filename);
							//system("pause");

						   }
			else{

				for(i=0;i<VN;i++){
								 h=hidden(V[i],h,Wh,VD,HD);
								r=reconstruction(h,r,Wr,VD,HD);
						   for(k=0;k<VD;k++){
														        if(k<(VD-1))fprintf(V_Re,"%f\t",r[k]);
																else fprintf(V_Re,"%f",r[k]);
											   }
                                               fprintf(V_Re,"\n");

				}
			}
				 fclose(W_R);


printf("\n-----------------------------------------------------\n");
printf("  Hidden values and weight matrixes were saved in:   \n");
printf("./%s\n",H_OUT_filename);
printf("./%s\n",W_H_filename);
printf("./%s\n",W_R_filename);
printf("./%s\n",V_Re_filename);
printf("-----------------------------------------------------\n");


///////////////////////////////////////////////////////////////
 int need_nomalization=0;
double **READ=d_CreateTwoDimensionalArray(XN,XD);
if((X_data=fopen(X_data_name,"r")) == NULL){
        fprintf(stderr, "Cannot open file [%s].\n", X_data_name);
		system("pause");
        exit(1);
    }
READ=READ_FILE2(X_data,READ,XD,XN);
 fclose(X_data);

for(i=0;i<XN;i++)
	for(j=0;j<XD;j++)
		if(READ[i][j]>0.998||READ[i][j]<(-0.998)) need_nomalization=1;

if(need_nomalization==1){
	 printf("---------------------\nNeed to do nomalization\n---------------------\n");
	 //////////////////////////////Normalization//////////////////////////////////////
		double **X=d_CreateTwoDimensionalArray(XN,XD);

		X=Normalization(READ,X,XN,XD);

		for(i=0;i<VN;i++){
			for(j=0,e=0;j<Wsize;j++)
				for(q=0;q<XD;q++,e++)

					V[i][e]=X[i+j][q];
			        V[i][VD]=1.0;
        }
	DeleteTwoDimensionalArray(X,XN,XD);
 }


/*********????SAE?ł?idowing(?K?i???Ȃ?)*******/
if(need_nomalization!=1) {
 printf("---------------------\nNeed not to do nomalization\n---------------------\n");
 for(i=0;i<VN;i++){
	 for(j=0,e=0;j<Wsize;j++)
		for(q=0;q<XD;q++,e++)

	 V[i][e]=READ[i+j][q];

        V[i][VD]=1.0;
        }
}

DeleteTwoDimensionalArray(READ,XN,XD);


 FILE *V_OUT;
 V_OUT=fopen(Windowing_DATA_filename,"w");
 for(i=0;i<VN;i++){
	 for(int k=0;k<VD;k++){
        if(k<VD-1)fprintf(V_OUT,"%f\t",V[i][k]);
		else fprintf(V_OUT,"%f",V[i][k]);
	 }
        fprintf(V_OUT,"\n");
                   }
fclose(V_OUT);


while(_getch()!=32);
}

}
////////////////隠れ層(特徴)を出力/////////////////////
FILE *H_OUT;



                            if((H_OUT=fopen(H_OUT_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",H_OUT_filename);
							//system("pause");
						    
						   }    
							else{ 
                             for(i=0;i<VN;i++){
                                              h=hidden(V[i],h,Wh,VD,HD);
                                               for(k=0;k<HD;k++){
                                                               if(k<HD-1) fprintf(H_OUT,"%f\t",h[k]);
															   else fprintf(H_OUT,"%f",h[k]);
											   }
                                               fprintf(H_OUT,"\n");
                                          }
						   
                              fclose(H_OUT);
							}
////////////////encoderの重み行列とバイアスを出力/////////////////////
FILE *W_H;


                            if((W_H=fopen(W_H_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_H_filename);
							//system("pause");
						    
						   }  
							else{  
                             for(i=0;i<VD+1;i++){
                                             
                                               for(k=0;k<HD;k++){
                                                                if(k<HD-1)fprintf(W_H,"%f\t",Wh[i][k]);
																else fprintf(W_H,"%f",Wh[i][k]);
											   }
                                               fprintf(W_H,"\n");
                                          }
						   
                              fclose(W_H);
							}
////////////////decoderの重み行列とバイアスを出力/////////////////////
FILE *W_R;


                           if((W_R=fopen(W_R_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_R_filename);
							//system("pause");
						    
						   }
						   else{  
                             for(i=0;i<HD+1;i++){
                                             
                                               for(k=0;k<VD;k++){
                                                                if(k<(VD-1))fprintf(W_R,"%f\t",Wr[i][k]);
																else fprintf(W_R,"%f",Wr[i][k]);
											   }
                                               fprintf(W_R,"\n");
                                          }
						   
                              fclose(W_R);
						   }
							

////////////////Reconstructionを出力/////////////////////
FILE *V_Re;	
				
 if((V_Re=fopen(V_Re_filename,"w"))==NULL){
						   printf("Cannot open file %s.\n",W_R_filename);
							//system("pause");
						    
						   }
			else{  
				
				for(i=0;i<VN;i++){						   
								 h=hidden(V[i],h,Wh,VD,HD);
								r=reconstruction(h,r,Wr,VD,HD);
						   for(k=0;k<VD;k++){
														        if(k<(VD-1))fprintf(V_Re,"%f\t",r[k]);
																else fprintf(V_Re,"%f",r[k]);
											   }
                                               fprintf(V_Re,"\n");

				}
			}
				 fclose(W_R);


printf("\n-----------------------------------------------------\n");
printf("  Hidden values and weight matrixes were saved in:   \n");
printf("./%s\n",H_OUT_filename);
printf("./%s\n",W_H_filename);
printf("./%s\n",W_R_filename);
printf("./%s\n",V_Re_filename);
printf("-----------------------------------------------------\n");

DeleteTwoDimensionalArray(V_denoising,VN,VD+1);
free(p_hat);
free(r);
free(h);
}


