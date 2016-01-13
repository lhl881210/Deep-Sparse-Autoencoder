double noise_variance(double x)
{
	double noise=tanh(x);   
	return pow(noise,4.0);
}



double error(double *v,double *r,int VD){
        int i;
        double e=0;
        for(i=0;i<VD;i++)
        {
			double temp=r[i]-v[i];

                e+=(temp*temp);
        }
        return e*0.5;
}

double *hidden(double *v,double *h,double **Wh,int VD,int HD){
       int i,j;
        //printf("h=\t");
        for(i=0;i<HD;i++)
        {
                h[i] = 0;

                for(j=0;j<VD+1;j++)
                {
                        h[i] += v[j]*Wh[j][i];
                }
                h[i] = tanh(h[i]);
                //zprintf("%f\t",h[i]);
        }
        h[HD]=1.0;
        //printf("\n");
        return h;
}

double *reconstruction(double *h,double *r,double **Wr,int VD,int HD){
       int i,j;
       //printf("r=\n");
       for(i=0;i<VD;++i)
        {
                r[i] = 0.0;
                for(j=0;j<HD+1;++j)
                {
                        r[i] += h[j]*Wr[j][i];
                }
                r[i] = tanh(r[i]);
                //printf("%f\t",r[i]);
        }
        //printf("\n");
	    r[VD]=1.0; 
        return r;
       }

double error_all(double **V,double *h,double *r,double **Wh,double **Wr,int VN,int VD,int HD){
       int i;
       double e_all=0;
	   double **R = d_CreateTwoDimensionalArray(VN, VD+1);

		#pragma omp parallel for private(i) schedule(dynamic)
	    for(i=0;i<VN;i++){
          h=hidden(V[i],h,Wh,VD,HD);
          R[i]=reconstruction(h,r,Wr,VD,HD);
          }
	   
for (i = 0; i<VN; i++)  e_all = e_all + error(V[i], R[i], VD);

DeleteTwoDimensionalArray(R, VN, VD);
  return e_all;
   }


double sparce(double p,double *p_hat,int VD){
       double s;
       int i;
       for(s=0,i=0;i<VD;i++){
                            
                            s=s+p*log(p/p_hat[i])+(1-p)*log((1-p)/(1-p_hat[i]));
                            }
      //for(s=0,i=0;i<VD;i++)
     // printf("p_had=%f\t",p_hat[i]);
     // printf("\n");
      
      
      
       return s;
       }
