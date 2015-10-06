#include<Rcpp.h>
#include <R.h>
#include "SMC_Back.h"

/****************************************************************************************************
				 FORWARD-BACKWARD SMC TO TARGET SMOOTHING DISTRIBUTION
****************************************************************************************************/

void SMC_ForBack(double *y, int dy, int dx, int T, double *theta, int N, int Nb, gsl_rng *random, ParamTransitionPF paramTransition, 
ResamplingPF resampling, SimInitPF simInit, TransitionPF transition, SimTransitionPF simTransition, PotentialPF potential,int *par, 
double *lik, double *expx)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	
	double *W=(double*)malloc(sizeof(double)*N);
	double *storeW=(double*)malloc(sizeof(double)*(N*T));

	double *x=(double*)malloc(sizeof(double)*n);
	double *storeX=(double*)malloc(sizeof(double)*(T*n));
	
	//FORWARD PASS

	SMC_Forward(y, dy, dx, T, theta, N,  random, paramTransition, resampling, simInit, simTransition, potential,
	par, lik, storeX, storeW, x, W);

	//BACKWARD PASS

	double *xT=(double*)malloc(sizeof(double)*(Nb*dx));

	SMC_Back(dx, T, 0, theta, N, Nb, random, resampling, transition, x, storeX, W, storeW, expx, xT);

	free(x);
	x=NULL;
	free(xT);
	xT=NULL;
	free(W);
	W=NULL;
	free(storeW);
	storeW=NULL;
	free(storeX);
	storeX=NULL;
	
}


/****************************************************************************************************
				 FORWARD-BACKWARD SMC TO TARGET SMOOTHING DISTRIBUTION
****************************************************************************************************/

void SMC_ForBackMarg(double *y, int dy, int dx, int T, double *theta, int N, int Nb, gsl_rng *random, ParamTransitionPF paramTransition, 
ResamplingPF resampling, SimInitPF simInit, TransitionPF transition, SimTransitionPF simTransition, PotentialPF potential,int *par, 
double *lik, double *expx)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	
	double *W=(double*)malloc(sizeof(double)*N);
	double *storeW=(double*)malloc(sizeof(double)*(N*T));

	double *x=(double*)malloc(sizeof(double)*n);
	double *storeX=(double*)malloc(sizeof(double)*(T*n));
	
	//FORWARD PASS

	SMC_Forward(y, dy, dx, T, theta, N,  random, paramTransition, resampling, simInit, simTransition, potential,
	par, lik, storeX, storeW, x, W);

	//BACKWARD PASS

	double *xT=(double*)malloc(sizeof(double)*(Nb*dx));

	SMC_BackMarg(dx, T, 0, theta, N, Nb, random, resampling, transition, x, storeX, W, storeW, expx, xT);

	free(x);
	x=NULL;
	free(xT);
	xT=NULL;
	free(W);
	W=NULL;
	free(storeW);
	storeW=NULL;
	free(storeX);
	storeX=NULL;
	
}

/***********************************************************************************************************************************/

/****************************************************************************************************
				BACKWARD STEP OF SMC (TARGET SMOOTHING DISTRIBUTION) 
****************************************************************************************************/


void SMC_Back(int dx, int T, int T_in, double *theta, int N, int Nb, gsl_rng *random, ResamplingPF resampling, TransitionPF transition, 
double *x, double *storeX, double *W, double *storeW, double *expx, double *xh)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	
	double *xhh=NULL, *wei=NULL;
	
	iter= T-1;

	for(i=0;i<T*dx;i++)
	{
		expx[i]=0;
	}

	(*resampling)(random, x, dx, N, Nb, W, xh);   

	for(k=0;k<dx;k++)
	{
		for(j=0;j<Nb;j++)
		{
			expx[dx*iter+k]+=xh[dx*j+k];
		}
		expx[dx*iter+k]/=(double)(Nb);
		
	}

	xhh=(double*)malloc(sizeof(double)*dx);
	wei=(double*)malloc(sizeof(double)*N);

	for(iter= T-2; iter>=T_in; iter--)
	{
		for(k=0;k<dx;k++)
		{
			for(i=0; i<N; i++)
			{
				x[dx*i+k]=storeX[n*iter+dx*i+k];
			}
		}

		for(j=0; j<Nb; j++)
		{

			for(k=0;k<dx;k++)
			{
				xhh[k]=xh[dx*j+k];
			}
			
			(*transition)(xhh, x, dx, iter+1, N, theta, wei);

			for(i=0;i<N;i++)
			{
				wei[i]+=storeW[N*iter+i];
			}

			cc=weight(wei,W,N); 

			//sample a new particle

			ResampleBack(gsl_rng_uniform_pos(random), x, dx, N, W, j, xh);   

			for(k=0;k<dx;k++)
			{
				expx[dx*iter+k]+=xh[dx*j+k]/((double) Nb);
			}

		}
	}


	free(xhh);
	xhh=NULL;
	free(wei);
	wei=NULL;
}


/****************************************************************************************************
				BACKWARD STEP OF SQMC (TARGET MARGINAL SMOOTHING DISTRIBUTIONS) 
****************************************************************************************************/

void SMC_BackMarg(int dx, int T, int T_in, double *theta, int N, int Nb, gsl_rng *random, ResamplingPF resampling, TransitionPF transition, 
double *x, double *storeX, double *W, double *storeW, double *expx, double *xh)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	
	double *xhh=NULL, *wei=NULL, *Wb=NULL;
	

	//BACKWARD PASS
	iter= T-1;

	for(i=0;i<T*dx;i++)
	{
		expx[i]=0;
	}

	Wb=(double*)malloc(sizeof(double)*N);

	ResampleMargBackMC(random, x, dx, N, Nb, W,  Wb, xh);
	
	for(k=0; k<dx; k++)
	{
		for(j=0;j<Nb;j++)
		{
			expx[dx*iter+k]+=xh[dx*j+k];
		}
		expx[dx*iter+k]/=(double)(Nb);
	}
	
	for(i=0;i<N;i++)
	{
		Wb[i]=storeW[N*iter+i];
	}

	cc=weight(Wb,Wb,N);

	xhh=(double*)malloc(sizeof(double)*dx);
	wei=(double*)malloc(sizeof(double)*N);

	//2.3 Iterate

	for(iter= T-2; iter>=T_in; iter--)
	{	
		for(i=0; i<N; i++)
		{
			for(k=0;k<dx;k++)
			{
				x[dx*i+k]=storeX[n*iter+dx*i+k];
			}
			W[i]=0;
		}

		for(j=0;j<N;j++)
		{
			for(k=0;k<dx;k++)		
			{
				xhh[k]=storeX[n*(iter+1)+dx*j+k];
			}

			(*transition)(xhh, x, dx, iter+1, N, theta, wei);

			for(i=0;i<N;i++)
			{
				wei[i]+=storeW[N*iter+i];
			}

			cc=weight(wei,wei,N);

			for(i=0;i<N;i++)
			{	
				W[i]+=exp(log(Wb[j])+log(wei[i]));
			}
		}

		for(i=0;i<N;i++)
		{
			Wb[i]=W[i];
		}

		//2.3.1 Sample particle

		ResampleMargBackMC(random, x, dx, N, Nb, W,  W, xh);

		//2.3.2 compute expectaion of x_{iter} w.r.t. the smoothing distribution at time "iter"

		for(k=0; k<dx; k++)
		{
			for(j=0;j<Nb;j++)
			{
				expx[dx*iter+k]+=xh[dx*j+k];
			}
			expx[dx*iter+k]/=(double)(Nb);
		}
		
		
	}
	free(xhh);
	xhh=NULL;
	free(wei);
	wei=NULL;
	free(Wb);
	Wb=NULL;
	
}


void SMC_2F(double *y, int tstar, int dy, int dx, int T, double *theta, double *theta2F, int N, int Nb,  gsl_rng *random,
 ParamTransitionPF paramTransition, ResamplingPF resampling, SimInitPF simInit, TransitionPF transition, SimTransitionPF simTransition, 
PotentialPF potential, ParamTransitionPF paramTransition2F, 
SimInitPF simInit2F, TransitionPF transition2F, SimTransitionPF simTransition2F, PotentialPF potential2F,  int* par, double *lik, double *expx)
{
	int i, j, k, iter, n=N*dx;
	
	double cc;

	//1. Forward step of SQMC
	
	double *W=(double*)malloc(sizeof(double)*N); 
	double *storeW=(double*)malloc(sizeof(double)*(N*T));

	double *x=(double*)malloc(sizeof(double)*n);
	double *storeX=(double*)malloc(sizeof(double)*(T*n));

	SMC_Forward(y, dy, dx, T, theta, N,  random, paramTransition, resampling, simInit, simTransition, potential,
	par, lik, storeX, storeW, x, W);

	
	//2. Estimation of Q_{tstar|T} using marginal backward step of SQMC
	
	double *xstar=(double*)malloc(sizeof(double)*(N*dx));

	SMC_BackMarg(dx, T, tstar, theta, N, N, random, resampling, transition, x, storeX, W, storeW, expx, xstar);


	//3. Generalized backward information filter

	double *W2F=(double*)malloc(sizeof(double)*N);
	double *storeW2F=(double*)malloc(sizeof(double)*(N*T));

	double *x2F=(double*)malloc(sizeof(double)*n);
	double *storeX2F=(double*)malloc(sizeof(double)*(T*n));
	
	SMC_BIF(y, dy, dx, T, theta2F, N,  random, paramTransition2F, resampling, simInit2F, simTransition2F, potential2F, par, 
	lik,  storeX2F, storeW2F, x2F, W2F);

	free(x2F);	
	x2F=NULL;
	free(W2F);
	W2F=NULL;

	//4. Sampling from the approximation of the smoothing distribution

	double *xh=(double*)malloc(sizeof(double)*(dx*Nb));
	double *xhSave=(double*)malloc(sizeof(double)*(dx*Nb));

	for(i=0; i<T*dx; i++)
	{
		expx[i]=0;
	}


	// for t=tstar

	iter=tstar;

	for(i=0;i<N;i++)
	{
		W[i]=(double)(1.0/(N));	
	}

	(*resampling)(random, xstar, dx, N, Nb, W, xh); 
	

	free(xstar);
	xstar=NULL;

	for(k=0; k<dx; k++)
	{
		for(j=0;j<Nb;j++)
		{
			expx[dx*iter+k]+=xh[dx*j+k];
		}
		expx[dx*iter+k]/=(double)(Nb);
	}

	for(k=0;k<dx;k++)
	{
		for(j=0;j<Nb;j++)
		{
			xhSave[dx*j+k]=xh[dx*j+k];
		}
	}

	
	double *xhh=(double*)malloc(sizeof(double)*dx);
	double *wei=(double*)malloc(sizeof(double)*N);
		
	// for t=tstar-1 until t=0

	for(iter= tstar-1; iter>=0; iter--)
	{
		//Set of particles generated at time iter by the forward setp
		for(k=0;k<dx;k++)
		{
			for(i=0; i<N; i++)
			{
				x[dx*i+k]=storeX[n*iter+dx*i+k];
			}
		}

		for(j=0; j<Nb; j++)
		{	
			//Set of particles \tilde{x}_{iter+1}^j

			for(k=0;k<dx;k++)
			{
				xhh[k]=xh[dx*j+k];
			}

			//2.3.1 Compute backward weights at time iter	

			(*transition)(xhh, x, dx, iter+1, N, theta, wei);

			for(i=0;i<N;i++)
			{	
				wei[i]+=storeW[N*iter+i];
			}

			cc=weight(wei,W,N); 

			//2.3.2 Sample particle

			ResampleBack(gsl_rng_uniform_pos(random), x, dx, N,  W, j, xh);   

			//2.3.3 compute expectaion of x_{iter} w.r.t. the smoothing distribution at time "iter"

			for(k=0;k<dx;k++)
			{
				expx[dx*iter+k]+=xh[dx*j+k]/((double) Nb);
			}
		}
	}
	free(storeW);
	storeW=NULL;

	free(storeX);
	storeX=NULL;

	// for t=tstar+1 until t=T-1

	for(k=0;k<dx;k++)
	{
		for(j=0;j<Nb;j++)
		{
			xh[dx*j+k]=xhSave[dx*j+k];
		}
	}

	free(xhSave);
	xhSave=NULL;

	for(iter= tstar+1; iter< T; iter++)
	{
		//Set of particles generated at time iter by the forward setp
		for(k=0;k<dx;k++)
		{
			for(i=0; i<N; i++)
			{
				x[dx*i+k]=storeX2F[n*iter+dx*i+k];
			}
		}

		for(j=0; j<Nb; j++)
		{	
			//Set of particles \tilde{x}_{iter+1}^j

			for(k=0;k<dx;k++)
			{
				xhh[k]=xh[dx*j+k];
			}

			//2.3.1 Compute backward weights at time iter	

			(*transition2F)(xhh, x, dx, T-iter, N, theta2F, wei);

			for(i=0;i<N;i++)
			{	
				wei[i]+=storeW2F[N*iter+i];
			}

			cc=weight(wei,W,N); 

			//2.3.2 Sample particle

			ResampleBack(gsl_rng_uniform_pos(random), x, dx, N,  W, j, xh);   

			//2.3.3 compute expectaion of x_{iter} w.r.t. the smoothing distribution at time "iter"

			for(k=0;k<dx;k++)
			{
				expx[dx*iter+k]+=xh[dx*j+k]/((double) Nb);
			}
		}
		
	}

	free(storeW2F);
	storeW2F=NULL;
	free(storeX2F);
	storeX2F=NULL;
	

	free(W);
	W=NULL;	
	free(x);
	x=NULL;

	free(xh);
	xh=NULL;
	free(xhh);
	xhh=NULL;

	free(wei);
	wei=NULL;

}


