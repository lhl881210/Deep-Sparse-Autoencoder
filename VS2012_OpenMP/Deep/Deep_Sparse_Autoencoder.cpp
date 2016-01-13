//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88
//                      (| -_- |)
//                      0\  =  /0
//                    ___/`---'\___
//                  .' \\|     |// '.
//                 / \\|||  :  |||// \
//                / _||||| -:- |||||- \
//               |   | \\\  -  /// |   |
//               | \_|  ''\---/''  |_/ |
//               \  .-\__  '-'  ___/-. /
//             ___'. .'  /--.--\  `. .'___
//          ."" '<  `.___\_<|>_/___.' >' "".
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
//         \  \ `_.   \_ __\ /__ _/   .-` /  /
//     =====`-.____`.___ \_____/___.-`___.-'=====
//                       `=---='
//
//
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//      Plead for Buddha bless, NO BUG
//
// ---------------------------------------------------------------------------
// (Deep) Sparse Denoising Autoencoder
//
// Programming by HaiLong LIU
// ---------------------------------------------------------------------------

#include <math.h>
#include "random.h"
#include "matrix.h"
#include "tools.h"
#include   <float.h> 
#include <time.h>
#include <stdio.h>  
#include <stdlib.h> 
#include <conio.h>
 char H_OUT_filename[250] = {'\0'};
 char W_H_filename[250] = {'\0'};
 char W_R_filename[250] = {'\0'};
 char V_Re_filename[250] = {'\0'};
 char Gradient_filename[250] = {'\0'};
 char Windowing_DATA_filename[250] = {'\0'};
double p;//Sparse target
double A;//L2 parameter
double B;//Sparse parameter
double eck;//error check(convergence condition)
double denoising_rate;//denoising rate in [0,1]
double *setting=d_CreateONEDimensionalArray(7); //Parameters reading from setting file
int XD=0,XN=0,i,j,e,q;
int VN=0,VD=0;
 FILE *X_data;
 char X_data_name[80];
 int Wsize=0;
 double min_Var;
#include "BP_function.h"
#include "BP_method.h"

///////////////////MAIN////////////////////////
int main(int argc, char *argv[])
{
srand( (unsigned)time( NULL ) );
printf("|---------------------------------------------------------|\n");
printf("|---------------------------------------------------------|\n");
printf("|           (Deep) Sparse Denoising Autoencoder           |\n");
printf("|               Programming by HaiLong LIU                |\n");
printf("|            Email: liu@em.ci.ritsumei.ac.jp              |\n");
printf("|                  Update: 2015.12.10                     |\n");
printf("|---------------------------------------------------------|\n\n");




 
 int need_nomalization=0;
 FILE *Setting_file;
if((Setting_file=fopen("setting.txt","r")) == NULL){
        printf("Cannot open file setting.txt.\n");
		system("pause");
        exit(1);
    }


fscanf(Setting_file,"%*[^\n]%*c");
setting=READ_FILE1(Setting_file,setting,6);
fclose(Setting_file);


//////////////////////////Parameters check///////////////////
p=setting[0];//Sparse target
A=setting[1];//L2 parameter
B=setting[2];//Sparse parameter
eck=setting[3];//error check(convergence condition)
denoising_rate=setting[4];//denoising rate in [0,1] 
printf("|---------------------------------------------------------|\n");
printf("|        Parameters reading from [.//setting.txt]         |\n");
printf("|---------------------------------------------------------|\n");
printf("Sparse taget: p=%f\n",setting[0]);
printf("L2norm intensity: A=%f\n",setting[1]);
printf("Sparse intensity: B=%f\n",setting[2]);
printf("Convergence condition: eck=%f\n",setting[3]);
printf("Noise rate: noise_rate=%f\n",setting[4]);
printf("Min_Var=%f\n\n",min_Var=setting[5]);



////////////////////////////     TYPE     ///////////////////////////
//     ,----------------------------------------------------,            
//     | [][][][][]  [][][][][]  [][][][]  [][__]  [][][][] |           
//     |                                                    |             
//     |  [][][][][][][][][][][][][][_]    [][][]  [][][][] |
//     |  [_][][][][][][][][][][][][][ |   [][][]  [][][][] |
//     | [][_][][][][][][][][][][][][]||     []    [][][][] |
//     | [__][][][][][][][][][][][][__]    [][][]  [][][]|| |
//     |   [__][________________][__]              [__][]|| |
//     `----------------------------------------------------'
//                            | |
//                            | |
//                          \+++++/
//                           \+++/
//                            \+/
//     					       +
 
printf("-----------------------------------------------------------\n");   
printf("|                 Please choose the mod                   |\n");
printf("-----------------------------------------------------------\n");
printf("|         [1] Train the new model                         |\n");
printf("|         [2] Continue to train the old model             |\n");
printf("-----------------------------------------------------------\n");

int Model=0;
printf("Model>>>");
scanf("%d",&Model);
if(Model!=1&&Model!=2){
printf("-----------------------------------------------------------\n");
printf("|                 Please choose the mod                   |\n");
printf("-----------------------------------------------------------\n");
printf("|         [1] Train the new model                         |\n");
printf("|         [2] Continue to train the old model             |\n");
printf("-----------------------------------------------------------\n");
printf("Model>>>");
scanf("%d",&Model);

}

printf("Learning set X File name ->>");
 scanf("%s",X_data_name);
 ////////////////////Check data's quantity of input(X)///////////////////////
  if((X_data=fopen(X_data_name,"r")) == NULL){
        fprintf(stderr, "Cannot open file [%s].\n", X_data_name);
		system("pause");
        exit(1);
    }
 XN=Check_NO(X_data);
 fclose(X_data);

////////////////////Check data's dimensions of input(X)///////////////////////
  if((X_data=fopen(X_data_name,"r")) == NULL){
        fprintf(stderr, "Cannot open file [%s].\n", X_data_name);
		system("pause");
        exit(1);
    }
 XD=Check_D(X_data);
 fclose(X_data);

 //////////////////////////////READ DATA//////////////////////////////////////
 double **READ=d_CreateTwoDimensionalArray(XN,XD);
 if((X_data=fopen(X_data_name,"r")) == NULL){
        fprintf(stderr, "Cannot open file [%s].\n", X_data_name);
		system("pause");
        exit(1);
    }
READ=READ_FILE2(X_data,READ,XD,XN);
 fclose(X_data);


///////////Is data need to be normalize? Check!////////////////////
for(i=0;i<XN;i++)
	for(j=0;j<XD;j++)
		if(READ[i][j]>0.998||READ[i][j]<(-0.998)) need_nomalization=1;


/////////////////////Set size of window(Wsize)///////////////////////////////
 
 printf("Size of time window ->>");
 scanf("%d",&Wsize);
 
 ///////////////////Create Visual layer(V)////////////////////////
 VN=XN-Wsize+1;printf("VN=%d\t",VN);
 VD=XD*Wsize;printf("VD=%d\n",VD);
 double **V=d_CreateTwoDimensionalArray(VN,VD+1);


 /*********Normalization and Windowing*******/
 if(need_nomalization==1){
	 printf("---------------------\nNomalization\n---------------------\n");
	 //////////////////////////////Normalization//////////////////////////////////////
		double **X=d_CreateTwoDimensionalArray(XN,XD);
		
		X=Normalization(READ,X,XN,XD);
		
		for(i=0;i<VN;i++){
			for(j=0,e=0;j<Wsize;j++)
				for(q=0;q<XD;q++,e++)

					V[i][e]=X[i+j][q];
			        V[i][VD]=1;
        }
	DeleteTwoDimensionalArray(X,XN,XD);
 }


/*********Windowing with out normalization*******/
if(need_nomalization!=1) {
 printf("---------------------\nNeed not to do nomalization\n---------------------\n");
 for(i=0;i<VN;i++){
	 for(j=0,e=0;j<Wsize;j++)
		for(q=0;q<XD;q++,e++)

	 V[i][e]=READ[i+j][q];

        V[i][VD]=1.0;
        }
}
/*************************************************/   

sprintf(Windowing_DATA_filename,"Widowing_DATA_%d.txt", VD);
sprintf(V_Re_filename,"Reconstruction_DATA_%d.txt", VD);



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
						
 ///////////////////Set dimensions of Hidden layer(H)////////////////////////
int HD;
 printf("Dimensionas of Hidden layer ->>");
 scanf("%d",&HD);


  sprintf(H_OUT_filename, "%d.txt",HD);
  sprintf(W_H_filename, "W_H_%d.txt", HD);
  sprintf(W_R_filename, "W_R_%d.txt", HD);
  sprintf(Gradient_filename, "Gradient_%d_to_%d.txt",VD,HD);

  FILE *BFGS;

if(Model==1){
  BFGS=fopen(Gradient_filename,"w");
  fclose(BFGS);
  BFGS=fopen(Gradient_filename,"a+");
  fprintf(BFGS,"e_all\tE_all\ttime(sec)\n");
  fclose(BFGS);
}

double *h=d_CreateONEDimensionalArray(HD+1);//Hidden layer

double**Wh=d_CreateTwoDimensionalArray(VD+1,HD);//Weights of encoder
double **Wr=d_CreateTwoDimensionalArray(HD+1,VD);//Weights of decoder

if(Model==1){
for(int ppp=1,i=0;i<VD+1;i++)
	  for(j=0;j<HD;j++){
					
		  Wh[i][j]=Guassian2(0,0.2,ppp);//d_rndn(99)/10000+0.0001; //Initialize weights of encoder by Guassian
		  if(ppp==1)ppp=0;
		  else if(ppp==0)ppp=1;
	 // printf("%f\n",Wh[i][j]);system("pause");
	  }

for(int ppp=1,i=0;i<HD+1;i++)
	  for(j=0;j<VD;j++){
		  Wr[i][j]=Guassian2(0,0.2,ppp);//d_rndn(99)/10000+0.0001;//Initialize weights of decoder by Guassian
#pragma omp parallel if (2)
		  if(ppp==1)ppp=0;
		  else if(ppp==0)ppp=1;
	  }
}

if(Model==2){
	FILE *READ_WH,*READ_WR;
	READ_WH=fopen(W_H_filename,"r");
	READ_WR=fopen(W_R_filename,"r");
	Wh=READ_FILE2(READ_WH,Wh,HD,VD+1);
	Wr=READ_FILE2(READ_WR,Wr,VD,HD+1);
	fclose(READ_WH);
	fclose(READ_WR);
}





///////////////////Training////////////////////////


BP_ALL(V,Wh,Wr,VN,VD,HD,eck,p,A,B);//Train SAE by BP method
			
  
 system("pause"); 
	return 0;
}
