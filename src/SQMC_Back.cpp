
#include <R.h>
#include "SQMC_Back.h"




/****************************************************************************************************
				BACKWARD STEP OF SQMC (TARGET SMOOTHING DISTRIBUTION) 
****************************************************************************************************/


void SQMC_Back(int dx, int T,  int T_in, double *theta, int N, int Nb, int qmcB, ResamplingBack_QMC resampling, TransitionPF transition, 
double *parPsi, double *x, double *storeX, double *W, double *storeW, double *expx, double *xh)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	
	double *xhh=NULL, *sim=NULL, *wei=NULL;
	int *J=NULL;

	for(i=0; i<T*dx; i++)
	{
		expx[i]=0;
	}

	//Generate a RQMC point set in (0,1)^T (nested scrambling is used)

	sim=(double*)malloc(sizeof(double)*(Nb*T));
	
	Scrambled *os1 = Scrambled_Create(qmcB, T, Nb);

	Scrambling(os1);

	getPoints(os1, T, Nb, sim);

	Scrambled_Destroy(os1);

	os1=NULL;
	
	//Sample  particles  at time T-1

	iter=T-1;

	J=(*resampling)(sim, T, parPsi, x, dx, N, Nb, W, xh);	

	for(k=0; k<dx; k++)
	{
		for(j=0;j<Nb;j++)
		{
			expx[dx*iter+k]+=xh[dx*j+k]/((double)Nb);
		}
	}


	free(J);
	J=NULL;
	
	xhh=(double*)malloc(sizeof(double)*dx);
	wei=(double*)malloc(sizeof(double)*N);
		
	//Start iterate

	for(iter= T-2; iter>=T_in; iter--)
	{
		//Set of particles generated at time iter by the forward setp, sorted according to their Hilbert index

		for(k=0;k<dx;k++)
		{
			for(i=0; i<N; i++)
			{
				x[dx*i+k]=storeX[n*iter+dx*i+k];
			}
		}

		//Sample  particles  at time iter

		for(j=0; j<Nb; j++)
		{	
			//Particle \tilde{x}_{iter+1}^j

			for(k=0;k<dx;k++)
			{
				xhh[k]=xh[dx*j+k];
			}

			//Compute backward weights at time iter	

			(*transition)(xhh, x, dx, iter+1, N, theta, wei);

			for(i=0;i<N;i++)
			{	
				wei[i]+=storeW[N*iter+i];
			}

			cc=weight(wei,W,N); 

			//Sample particle j

			ResampleBack(sim[T*j+T-1-iter], x, dx, N,  W, j, xh);   

			//Compute smoothing expectation

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
	free(sim);
	sim=NULL;
}



/****************************************************************************************************
				BACKWARD STEP OF SQMC (TARGET MARGINAL SMOOTHING DISTRIBUTION) 
****************************************************************************************************/


void SQMC_BackMarg(int dx, int T, int T_in, double *theta, int N, int Nb, int qmc, ResamplingBack_QMC resampling, TransitionPF transition, 
double *parPsi, double *x, double *storeX, double *W, double *storeW, double *expx, double *xh)
{
	int i, j, k, iter, n=N*dx;
	double cc;
	int dimSeq=1;

	double *xhh=NULL, *sim=NULL, *wei=NULL, *Wb=NULL;
	
	sim=(double*)malloc(sizeof(double)*(Nb*dimSeq));
	
	Scrambled *os1 = Scrambled_Create(qmc, dimSeq, Nb);

	for(i=0; i<T*dx; i++)
	{
		expx[i]=0;
	}

	iter=T-1;

	//Generate RQMC points in (0,1) (nested scrambling is used)

	Scrambling(os1);

	getPoints(os1, dimSeq, Nb, sim);

	//Sample particles at tim T-1
	
	Wb=(double*)malloc(sizeof(double)*N);

	int *J=(*resampling)(sim, dimSeq, parPsi, x, dx, N, Nb, W, xh);	

	free(J);
	J=NULL;

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
		//Generate RQMC points in (0,1) (nested scrambling is used)

		Scrambling(os1);

		getPoints(os1, dimSeq, Nb, sim);

		//Reweights the particles generated at time iter in the forward pass

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

		//Sample particles at time iter

		ResampleMargBackQMC(sim, x, dx, N, Nb, W,  W, xh);

		//Compute Smoothing expectation

		for(k=0; k<dx; k++)
		{
			for(j=0;j<Nb;j++)
			{
				expx[dx*iter+k]+=xh[dx*j+k];
			}
			expx[dx*iter+k]/=(double)(Nb);
		}

	}

	Scrambled_Destroy(os1);
	os1=NULL;

	free(sim);
	sim=NULL;
	free(xhh);
	xhh=NULL;
	free(Wb);
	Wb=NULL;
	free(wei);
	wei=NULL;
}

/***********************************************************************************************************************************/



/****************************************************************************************************
				SQMC GENERALIZED INFORMATION FILTER	
****************************************************************************************************/


void SQMC_2F(double *y, int tstar, int dy, int dx, int T, double *theta, double *theta2F, int N, int Nb, int qmc,
 ParamTransitionPF paramTransition, ResamplingBack_QMC resamplingFor, ResamplingBack_QMC resamplingBack, SimInitPF_QMC simInit, 
TransitionPF transition, SimTransitionPF_QMC simTransition, PotentialPF potential, ParamTransitionPF paramTransition2F, 
SimInitPF_QMC simInit2F, TransitionPF transition2F, SimTransitionPF_QMC simTransition2F, PotentialPF potential2F,  int* par, 
double *parPsi, double *lik, double *expx)
{
	int i, j, k, iter, n=N*dx;

	double cc;

	//1. Forward step of SQMC
	
	double *W=(double*)malloc(sizeof(double)*N); 
	double *storeW=(double*)malloc(sizeof(double)*(N*T));

	double *x=(double*)malloc(sizeof(double)*n);
	double *storeX=(double*)malloc(sizeof(double)*(T*n));

	SQMC_Forward(y, dy, dx, T, theta, N,  qmc, paramTransition, resamplingFor,simInit, simTransition,potential, 
	par, parPsi, lik, storeX, storeW, x, W);
	
	//2. Estimation of Q_{tstar|T} using marginal backward step of SQMC
	
	double *xstar=(double*)malloc(sizeof(double)*(N*dx));

	SQMC_BackMarg(dx, T, tstar, theta, N, N, qmc, resamplingBack, transition, parPsi, x, storeX, W, storeW,expx, xstar);

	//3. Generalized backward information filter

	double *W2F=(double*)malloc(sizeof(double)*N);
	double *storeW2F=(double*)malloc(sizeof(double)*(N*T));

	double *x2F=(double*)malloc(sizeof(double)*n);
	double *storeX2F=(double*)malloc(sizeof(double)*(T*n));
	
	SQMC_BIF(y, dy, dx, T, theta2F, N,  qmc, paramTransition2F, resamplingFor, simInit2F, simTransition2F, potential2F, 
	par, parPsi, lik, storeX2F, storeW2F, x2F, W2F);

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

	
	double *sim=(double*)malloc(sizeof(double)*(Nb*T));
	
	Scrambled *os1 = Scrambled_Create(qmc, T, Nb);

	Scrambling(os1);

	getPoints(os1, T, Nb, sim);

	Scrambled_Destroy(os1);

	os1=NULL;

	// for t=tstar

	iter=tstar;

	for(i=0;i<N;i++)
	{
		W[i]=(double)(1.0/(N));	
	}

	int *J=(*resamplingFor)(sim, T, parPsi, xstar, dx, N, Nb, W, xh);

	free(J);
	J=NULL;
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

	int count=1;

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

			ResampleBack(sim[T*j+count], x, dx, N,  W, j, xh);   

			//2.3.3 compute expectaion of x_{iter} w.r.t. the smoothing distribution at time "iter"

			for(k=0;k<dx;k++)
			{
				expx[dx*iter+k]+=xh[dx*j+k]/((double) Nb);
			}
		}
		count++;
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

			ResampleBack(sim[T*j+count], x, dx, N,  W, j, xh);   

			//2.3.3 compute expectaion of x_{iter} w.r.t. the smoothing distribution at time "iter"

			for(k=0;k<dx;k++)
			{
				expx[dx*iter+k]+=xh[dx*j+k]/((double) Nb);
			}
		}
		count++;
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

	free(sim);
	sim=NULL;

}

































































