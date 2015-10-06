

#include "Model_LG.h"


/******************************************************************************************************
		PRELIMIARIES (USED FOR GAUSSIAN PRIOR DISTRIBUTION)
********************************************************************************************************/

void simLinearGaussian(gsl_rng *random, int dx, int start1, int start2, int N, double *theta, double *x)
{
	int i, s,k;
	
	double *mu=(double*)malloc(sizeof(double)*dx);
	double *wSim=(double*)malloc(sizeof(double)*dx);

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);

	//mean

	for(k=0;k<dx;k++)
	{
		mu[k]=theta[start1+k];
	}

	//Sigma 

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[start2+dx*k+s]);
		}
	}

	gsl_linalg_cholesky_decomp(SIGMA); 
 
	for(k=0;k<dx-1;k++)
	{
		for(s=k+1;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,0);
		}
	}

	//Simulations

	for(i=0;i<N;i++)
	{
		rmnorm(random,dx, mu, SIGMA,wSim);

		for(k=0;k<dx;k++)
		{
			x[dx*i+k]=wSim[k];
		}
	}
	
	free(wSim);
	wSim=NULL;
	free(mu);
	mu=NULL;
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;

}



void simLinearGaussianQMC(double* qseq, int dx, int start1, int start2, int N, double *theta, double *x)
{
	int i, s,k;
	double var, var1,muI;	

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);

	gsl_matrix *winv=NULL, *work1=NULL;
	gsl_vector *v1=NULL, *v2=NULL;
	gsl_permutation *p=NULL;

	
	//Sigma 

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[start2+dx*k+s]);
		}
	}
	
	var=pow(gsl_matrix_get(SIGMA,0,0),0.5);

	for(i=0;i<N;i++)			//sample first dimension
	{
		x[dx*i]=var*gsl_cdf_ugaussian_Pinv(qseq[dx*i])+theta[start1];
	}

	for(k=1;k<dx;k++)			//sample other dimensions
	{
		v1=gsl_vector_alloc(k);
		v2=gsl_vector_alloc(k);
		p=gsl_permutation_alloc(k);
		winv=gsl_matrix_alloc(k,k); 
		work1=gsl_matrix_alloc(k,k); 

		for(s=0;s<k;s++)
		{
			gsl_vector_set(v1,s,gsl_matrix_get(SIGMA,k,s));		//vector of covariance
		}

		gsl_matrix_view cov= gsl_matrix_submatrix(SIGMA, 0, 0, k,k);	//variance-covariance matrix
		gsl_matrix_memcpy(work1,&cov.matrix);

		gsl_linalg_LU_decomp(work1, p, &s);            
		gsl_linalg_LU_invert(work1, p, winv);				//inverse of the variance-covariance matrix

		gsl_blas_dsymv (CblasUpper, 1, winv, v1, 0, v2);		//Sima_{i,i}^{-1}Sigma_{i,j}

		gsl_blas_ddot (v2, v1, &var1);				        //Sigma_{j,i}Sima_{i,i}^{-1}Sigma_{i,j}
		
		var=pow(gsl_matrix_get(SIGMA,k,k)-var1,0.5);			//variance of x_k|x_{1:k-1}
			
		for(i=0;i<N;i++)
		{
			muI=theta[start1+k];

			for(s=0;s<k;s++)
			{
				muI+=gsl_vector_get(v2,s)*(x[dx*i+s]-theta[start1+s]);
			}
		
			x[dx*i+k]=var*gsl_cdf_ugaussian_Pinv(qseq[dx*i+k])+muI;
		}
		
		gsl_vector_free(v1);
		v1=NULL;
		gsl_vector_free(v2);
		v2=NULL;
		gsl_matrix_free(winv);
		winv=NULL;
		gsl_matrix_free(work1);
		work1=NULL;
		gsl_permutation_free(p);
		p=NULL;
	}
	
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;


}

/******************************************************************************************************
			GAUSSIAN PRIOR DISTRIBUTIONS
********************************************************************************************************/

void simPrior_LG(gsl_rng *random, int dx, int dy, int N, double *theta, double *x)
{
	int start2=dx+dy+2*pow(dx,2)+dy*dx+pow(dy,2);
	
	simLinearGaussian(random, dx, 0, start2, N, theta, x);
}

void simPriorQMC_LG(double* qseq, int dx, int dy, int N, double *theta, double *x)
{
	int start2=dx+dy+2*pow(dx,2)+dy*dx+pow(dy,2);
	
	simLinearGaussianQMC(qseq, dx, 0, start2, N, theta, x);

}



/******************************************************************************************************
			LINEAR GAUSSIAN MARKOV KERNEL
********************************************************************************************************/

void simTransition_LG(gsl_rng *random, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i,k,s, npx=pow(dx,2);
	double muI;
	double *mu=NULL, *wSim=NULL;

	double **M = (double**)malloc (sizeof (double)*N);

	for(i=0;i<N;i++)
	{
		M[i] =  (double*)malloc (sizeof (double)*dx);
	}

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);

	//Compute  A*X+mux

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx;k++)
		{
			M[i][k]=param[k];
			for(s=0;s<dx;s++)
			{
				M[i][k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}
	}
		
	//Choletsky decomposition of Sigma

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,param[dx+npx+dx*k+s]);
		}
	}


	//Simulation

	mu=(double*)malloc(sizeof(double)*dx);
	wSim=(double*)malloc(sizeof(double)*dx);

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx;k++)
		{
			mu[k]=M[i][k];
		}

		rmnorm(random,dx, mu, SIGMA, wSim);

		for(k=0;k<dx;k++)
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


void simTransitionQMC_LG(double *qseq, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i,k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx, start2=start1+(dx-1)*dx*0.5;
	double muI;
	double **M = (double**)malloc (sizeof (double)*N);

	for(i=0;i<N;i++)
	{
		M[i] =  (double*)malloc (sizeof (double)*dx);
	}
	
	//Compute  A*X+mux

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx;k++)
		{
			M[i][k]=param[k];
			for(s=0;s<dx;s++)
			{
				M[i][k]+=param[dx*(k+1)+s]*(xh[dx*i+s]-param[s]);
			}
		}
	}

	for(i=0;i<N;i++)			//sample first dimension
	{
		x[dx*i]=param[start2]*gsl_cdf_ugaussian_Pinv(qseq[(dx+1)*i+1])+M[i][0];
	}


	for(k=1;k<dx;k++)	//sample other dimensions
	{
		for(i=0;i<N;i++)
		{
			muI=M[i][k];

			for(s=0;s<k;s++)
			{
				muI+=param[start1+((k-1)*k)/2+s]*(x[dx*i+s]-M[i][s]);
			}
			if(!((k==(dx-1) && (i==(N-1)))))
			x[dx*i+k]=param[start2+k]*gsl_cdf_ugaussian_Pinv(qseq[(dx+1)*i+1+k])+muI;
			
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





/*************************************************************************************
		PARAMETERS FOR THE GAUSSIAN LINEAR PROPOSAL
*************************************************************************************/


double* paramTrans_LG(double *theta, int dx, int dy)
{
	int i,k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx;
	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);

	double *param=(double*)malloc(sizeof(double)*(start1+npx));

	//Store mux 

	for(k=0;k<dx;k++)
	{
		param[k]=theta[k];
	}

	//store A matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=theta[dx+dy+dx*k+s];
		}
	}

	//store parameters for the variance

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+dy+npx+mpy+dx*k+s]);
		}
	}

	gsl_linalg_cholesky_decomp(SIGMA); 
 
	for(k=0;k<dx-1;k++)
	{
		for(s=k+1;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,0);
		}
	}

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[start1+dx*k+s]=gsl_matrix_get(SIGMA,k,s);
		}
	}
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;


	return(param);

}


double* paramTransQMC_LG(double *theta, int dx, int dy)
{
	int k,s, npx=pow(dx,2), mpy=dy*dx;
	int start1=dx+npx, start2=start1+(dx-1)*dx*0.5;
	double var;

	double *param=(double*)malloc(sizeof(double)*(start2+dx));

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);
	gsl_matrix *winv=NULL, *work=NULL;
	gsl_vector *v1=NULL, *v2=NULL;
	gsl_permutation *p=NULL;
	

	//Store mux vector

	for(k=0;k<dx;k++)
	{
		param[k]=theta[k];
	}

	//store A matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=theta[dx+dy+dx*k+s];
		}
	}

	//store parameters for the variance

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+dy+npx+mpy+dx*k+s]);
		}
	}

	param[start2]=pow(gsl_matrix_get(SIGMA,0,0),0.5);

	for(k=1;k<dx;k++)			
	{
		v1=gsl_vector_alloc(k);
		v2=gsl_vector_alloc(k);
		p=gsl_permutation_alloc(k);
		winv=gsl_matrix_alloc(k,k);
		work=gsl_matrix_alloc(k,k); 

		for(s=0;s<k;s++)
		{
			gsl_vector_set(v1,s,gsl_matrix_get(SIGMA,k,s));		//vector of covariance
		}

		gsl_matrix_view cov= gsl_matrix_submatrix(SIGMA, 0, 0, k,k);	//variance-covariance matrix
		gsl_matrix_memcpy(work,&cov.matrix);

		gsl_linalg_LU_decomp(work, p, &s);            
		gsl_linalg_LU_invert(work, p, winv);				//inverse of the variance-covariance matrix

		gsl_blas_dsymv (CblasUpper, 1, winv, v1, 0, v2);		//Sima_{i,i}^{-1}Sigma_{i,j}

		gsl_blas_ddot (v2, v1, &var);				        //Sigma_{j,i}Sima_{i,i}^{-1}Sigma_{i,j}
		
		param[start2+k]=pow(gsl_matrix_get(SIGMA,k,k)-var,0.5);		//variance of x_k|x_{1:k-1}

		for(s=0;s<k;s++)
		{
			param[start1+((k-1)*k)/2+s]=gsl_vector_get(v2,s);
		}	

		gsl_vector_free(v1);
		v1=NULL;
		gsl_vector_free(v2);
		v2=NULL;
		gsl_matrix_free(winv);
		winv=NULL;
		gsl_matrix_free(work);
		work=NULL;
		gsl_permutation_free(p);
		p=NULL;
	}
	
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;

	return(param);
}









