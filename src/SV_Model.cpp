

#include "Model_SV.h"
#include <R.h>


double dmnorm_SV(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, const gsl_matrix *XVAR);


/******************************************************************************************************
			PRIOR DISTRIBUTIONS
********************************************************************************************************/

//// [[Rcpp::export]]
void simPrior_SV(gsl_rng *random, int dx, int dy, int N, double *theta, double *x)
{
	int start=dx+5*pow(dx,2)+1;
	
	simLinearGaussian(random, dx, start, start+dx, N, theta, x);
}
//// [[Rcpp::export]]
void simPriorQMC_SV(double* qseq, int dx, int dy, int N, double *theta, double *x)
{
	int start=dx+5*pow(dx,2)+1;
	
	simLinearGaussianQMC(qseq, dx, start, start+dx, N, theta, x);
}



/******************************************************************************************************
				MEASURE EQUATION
********************************************************************************************************/

//// [[Rcpp::export]]
void potential_SV(double *y, double *x, double *xh, int dy,  int dx, int time, int N, double *theta, double *result)
{
	int i,k,j,s, npx=pow(dx,2);
	double work1, work2, work3;

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);;
	gsl_matrix *WORK=gsl_matrix_calloc(dx,dx);
	gsl_vector *v2=gsl_vector_alloc(dx);
	gsl_vector *v3=gsl_vector_alloc(dx);

	
	//compute SigmaY
	
	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+2*npx+dx*k+s]);	
		}
	}

	//point where the density is computed

	for(k=0;k<dx;k++)
	{
		gsl_vector_set(v3,k,y[dx*time+k]);
	}


	if(time==0)
	{
		for(i=0;i<N;i++)
		{
			for(k=0;k<dx;k++)	//compute A*(h_{t-1}-mu)
			{
				work1=0;
				for(j=0;j<dx;j++)
				{
					work1+=theta[dx+npx+dx*k+j]*(x[dx*i+j]-theta[j]);
				}

				gsl_vector_set(v2,k,work1*exp(0.5*x[dx*i+k]));
				gsl_matrix_set(WORK,k,k,exp(-0.5*x[dx*i+k]));
			}
		
			result[i]=dmnorm_SV(dy, v3, v2, SIGMA, WORK)-0.5*theta[dx+5*npx];

			for(k=0;k<dx;k++)
			{
				result[i]+= -0.5*x[dx*i+k];
			}

		}
	}
	else
	{
		for(i=0;i<N;i++)
		{
			
			for(k=0;k<dx;k++)
			{
				work1=0;
				for(j=0;j<dx;j++)
				{
					work1+=theta[dx+npx+dx*k+j]*(x[dx*i+j]-theta[j]);
				}

				work2=0;

				for(j=0;j<dx;j++)
				{
					work2+=theta[dx+dx*k+j]*(xh[dx*i+j]-theta[j]);
				}

				gsl_vector_set(v2, k, exp(0.5*x[dx*i+k])*(work1-work2));
	
				gsl_matrix_set(WORK,k,k, exp(-0.5*x[dx*i+k]));
			}

			result[i]=dmnorm_SV(dy, v3, v2, SIGMA, WORK)-0.5*theta[dx+5*npx]; 

			for(k=0;k<dx;k++)
			{
				result[i]+= -0.5*x[dx*i+k];
			}


		}
	}
	
	gsl_vector_free(v2);
	v2=NULL;
	gsl_vector_free(v3);
	v3=NULL;
	gsl_matrix_free(SIGMA);
	SIGMA=NULL;
	gsl_matrix_free(WORK);
	WORK=NULL;


}



/******************************************************************************************************
				PROPOSAL DISTRIBUTION
********************************************************************************************************/


//// [[Rcpp::export]]
double* paramTrans_SV(double *theta, int dx, int dy)
{
	int i,k,s, npx=pow(dx,2);
	int start1=dx+npx;

	double *param=(double*)malloc(sizeof(double)*(start1+npx));
	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);

	//Store mux 

	for(k=0;k<dx;k++)
	{
		param[k]=theta[k];
	}

	//store Phi matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=theta[dx+4*npx+dx*k+s];
		}
	}

	
	//store parameters for the variance 

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+3*npx+dx*k+s]);
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

	
//// [[Rcpp::export]]
double* paramTransQMC_SV(double *theta, int dx, int dy)
{
	int k,s, npx=pow(dx,2);
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

	//store Phi matrix

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			param[dx+dx*k+s]=theta[dx+4*npx+dx*k+s];
		}
	}


	//store parameters for the variance

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+3*npx+dx*k+s]);
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







//// [[Rcpp::export]]
double dmnorm_SV(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *winv, const gsl_matrix *XVAR)
{
	int s;
	double ax,ay;
	gsl_vector *ym, *xm;
	gsl_matrix *work = gsl_matrix_alloc(n,n);

	gsl_blas_dsymm (CblasLeft, CblasUpper, 1, XVAR, winv, 0, work);

	gsl_blas_dsymm (CblasRight, CblasUpper, 1, XVAR, work, 0, work);

	xm = gsl_vector_alloc(n);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, mean );
	ym = gsl_vector_alloc(n);
	gsl_blas_dsymv(CblasUpper,1.0,work,xm,0.0,ym);
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);
	gsl_matrix_free(work);
	
	ay = -0.5*ay-0.5*n*log(2*M_PI);

	return ay;
}





//// [[Rcpp::export]]
void transition_SV(double *x, double *xh, int dx, int time, int N, double *theta, double *wei)
{
	int i,k,s, npx=pow(dx,2), mpy=dx*dx;
	double work;

	gsl_matrix *SIGMA=gsl_matrix_alloc(dx,dx);
	gsl_vector *mean=gsl_vector_alloc(dx);
	gsl_vector *point=gsl_vector_alloc(dx);
	
	for(k=0;k<dx;k++)
	{
		gsl_vector_set(mean,k,0);
	}

	//Variance

	for(k=0;k<dx;k++)
	{
		for(s=0;s<dx;s++)
		{
			gsl_matrix_set(SIGMA,k,s,theta[dx+3*npx+dx*k+s]);
		}
	}

	for(i=0;i<N;i++)
	{
		for(k=0;k<dx;k++)
		{
			work=theta[k];
			for(s=0;s<dx;s++)
			{	
				work+=theta[dx+4*npx+dx*k+s]*(xh[dx*i+s]-theta[s]);
			}	
			gsl_vector_set(point,k,x[k]-work);
		}

		wei[i]=dmnorm_C(dx, point, mean, SIGMA);
	}

	gsl_matrix_free(SIGMA);
	SIGMA=NULL;
	gsl_vector_free(mean);
	mean=NULL;
	gsl_vector_free(point);
	point=NULL;
	

}







































