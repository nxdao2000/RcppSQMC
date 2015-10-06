

#include "Model_Neuro.h"




double logfactorial(int n)
{
	int c;
	double fact = 0;
	
 	for (c = 2; c <= n; c++)
	{
    		fact+=log(c);
	}
	return(fact);
}



void potential_Neuro(double *y, double *x, double *xh, int dy,  int dx, int time, int N, double *theta, double *result)
{
	int i,j,k;
	int npx=pow(dx,2);
	double lambda, fy;
	
	
	for(j=0;j<dy;j++)
	{
		fy=logfactorial((int)y[dy*time+j]);

		for(i=0;i<N;i++)
		{
			if(j==0)
			{
				result[i]=0;
			}
			lambda=theta[dx+j];
			for(k=0;k<dx;k++)
			{
				lambda+=theta[dx+dy+npx+j*dx+k]*x[dx*i+k];
			}
	
			if(lambda>700)
			{
				result[i]+=-exp(700)+y[dy*time+j]*lambda-fy;
			}
			else
			{
				result[i]+=-exp(lambda)+y[dy*time+j]*lambda-fy;
			}
		}
	}

}





void simTransition_Neuro(gsl_rng *random, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i,k,s, npx=pow(dx,2), ddx=dx/2;
	double muI;
	double *mu=NULL, *wSim=NULL;

	double **M = (double**)malloc (sizeof (double)*N);

	for(i=0;i<N;i++)
	{
		M[i] =  (double*)malloc (sizeof (double)*dx);
	}

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx/2,dx/2);

	//Compute  A*X+mux

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx/2;k++)
		{
			M[i][k]=param[k];
			for(s=0;s<dx;s++)
			{
				M[i][k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}

		for(k=dx/2;k<dx;k++)
		{
			x[dx*i+k]=param[k];
			for(s=0;s<dx;s++)
			{
				x[dx*i+k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}
	}
		
	//Choletsky decomposition of Sigma

	for(k=0;k<dx/2;k++)
	{
		for(s=0;s<dx/2;s++)
		{
			gsl_matrix_set(SIGMA,k,s,param[dx+npx+dx*k/2+s]);
		}
	}


	//Simulation

	mu=(double*)malloc(sizeof(double)*(dx/2));
	wSim=(double*)malloc(sizeof(double)*(dx/2));

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx/2;k++)
		{
			mu[k]=M[i][k];
		}

		rmnorm(random,ddx, mu, SIGMA, wSim);

		for(k=0;k<dx/2;k++)
		{
			x[dx*i+k]=wSim[k];
		}
	}
	
	free(mu);
	mu=NULL;
	free(wSim);
	wSim=NULL;
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;



	for (i = 0; i < N; i++)
	{
  		free(M[i]);
	}
	free (M);
	
}


void simTransitionQMC_Neuro(double *qseq, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i,k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx;
	double muI;
	double **M = (double**)malloc (sizeof (double)*N);

	for(i=0;i<N;i++)
	{
		M[i] =  (double*)malloc (sizeof (double)*dx);
	}
	
	//Compute  A*X+mux

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx/2;k++)
		{
			M[i][k]=param[k];
			for(s=0;s<dx;s++)
			{
				M[i][k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}

		for(k=dx/2;k<dx;k++)
		{
			x[dx*i+k]=param[k];
			for(s=0;s<dx;s++)
			{
				x[dx*i+k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}
	}

	for(i=0;i<N;i++)			//sample first dimension
	{
		x[dx*i]=param[start1]*gsl_cdf_ugaussian_Pinv(qseq[(dx+1)*i+1])+M[i][0];
	}


	for(k=1;k<dx/2;k++)	//sample other dimensions
	{
		for(i=0;i<N;i++)
		{
			x[dx*i+k]=param[start1+k]*gsl_cdf_ugaussian_Pinv(qseq[(dx+1)*i+1+k])+M[i][k];
		}	
	}
	
	
	for (i = 0; i < N; i++)
	{
  		free(M[i]);
		M[i]=NULL;
	}
	free (M);
	M=NULL;


}






double* paramTrans_Neuro(double *t, int dx, int dy)
{
	int i,k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx;
	gsl_matrix *SIGMA=gsl_matrix_calloc(dx/2,dx/2);

	double *param=(double*)malloc(sizeof(double)*(start1+npx/4));

	//Store mux 

	for(k=0;k<dx;k++)
	{
		param[k]=t[k];
	}

	//store A matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=t[dx+dy+dx*k+s];
		}
	}

	//store parameters for the variance

	for(k=0;k<dx/2;k++)
	{
		gsl_matrix_set(SIGMA,k,k,t[dx+dy+npx+mpy+dx*k+k]);
		
	}

	gsl_linalg_cholesky_decomp(SIGMA); 
 

	for(k=0;k<dx/2;k++)
	{
		for(s=0;s<dx/2;s++)
		{
			param[start1+dx*k/2+s]=gsl_matrix_get(SIGMA,k,s);
		}
	}
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;


	return(param);

}


double* paramTransQMC_Neuro(double *t, int dx, int dy)
{
	int k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx, start2=start1+dx/2;
	double var;

	double *param=(double*)malloc(sizeof(double)*(start2));


	//Store mux vector

	for(k=0;k<dx;k++)
	{
		param[k]=t[k];
	}

	//store A matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=t[dx+dy+dx*k+s];
		}
	}

	//store parameters for the variance

	for(k=0;k<dx/2;k++)
	{
		param[start1+k]=pow(t[dx+dy+npx+mpy+dx*k+k],0.5);
	}

	

	return(param);
}













