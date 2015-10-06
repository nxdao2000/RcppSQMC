
#include<Rcpp.h>
#include "SMC.h"
#include <R.h>
#include "functions.h"


/****************************************************************************************************
				SEQUENTIAL MONTE CARLO

Output: Estimate of the log-likelihood function at each time steps ans filtering expectations
******************************************************************************************************/
void SMC(double *y, int dy, int dx, int T, double *theta, int N, gsl_rng *random, ParamTransitionPF paramTransition, ResamplingPF resampling, SimInitPF simInit, SimTransitionPF simTransition, PotentialPF potential, int *par, double *lik, double *expx)
{
    int i, k,iter;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *W=(double*)malloc(sizeof(double)*N);
    double *x=(double*)malloc(sizeof(double)*(N*dx));
    double *xh=(double*)malloc(sizeof(double)*(N*dx));
    double *param=NULL;

    iter=0;

    //Sample from initial distribution

    (*simInit)(random, dx, dy, N, theta, x);

    //Evaluate the potential function

    (*potential)(y, x, xh, dy, dx, iter, N, theta, wei);

    //Normalize the weights and compute log-likelihood

    lik[0]=weight(wei,W,N)-log(N);

    //Compute filtering expectation

    if(par[0]==1)
    {
        for(k=0;k<dx;k++)
        {
            cc=0;
            for(i=0;i<N;i++)
            {
                cc+=W[i]*x[dx*i+k];
            }
            expx[k]=cc;
        }
    }

    //Compute parameters of the Markov transition

    param=(*paramTransition)(theta,dx,dy);

    for(iter=1; iter< T; iter++)
    {
        //Resampling

        (*resampling)(random, x, dx, N, N, W, xh);

        //Mutation

        (*simTransition)(random, y, dy, dx, iter, N, param, xh, x);

        //Evaluate the potential function

        (*potential)(y, x, xh, dy, dx, iter, N, theta, wei);

        //Nomalize the weights and compute log-likelihood

        lik[iter]=lik[iter-1]+weight(wei,W,N)-log(N);

        //Compute filtering expectation

        if(par[0]==1)
        {
            for(k=0;k<dx;k++)
            {
                cc=0;
                for(i=0;i<N;i++)
                {
                    cc+=W[i]*x[dx*i+k];
                }
                expx[dx*iter+k]=cc;
            }
        }

    }

    free(param);
    param=NULL;
    free(wei);
    wei=NULL;
    free(W);
    W=NULL;
    free(x);
    x=NULL;
    free(xh);
    xh=NULL;

}

/****************************************************************************************************
FORWARD STEP OF THE FORWARD-BACKWARD SMOOTHING ALGORITHM
******************************************************************************************************/


void SMC_Forward(double *y, int dy, int dx, int T, double *theta, int N,  gsl_rng *random, ParamTransitionPF paramTransition, ResamplingPF resampling, SimInitPF simInit, SimTransitionPF simTransition, PotentialPF potential, int *par, double *lik,  double *storeX, double *storeW, double *x, double *W)
{
    int i, j, k, iter, n=N*dx;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *xh=(double*)malloc(sizeof(double)*n);

    double *param=NULL;

    iter=0;

    //Sample from initial distribution

    (*simInit)(random, dx, dy, N, theta, x);

    //Evaluate the potential function

    (*potential)(y, x, xh, dy, dx, iter, N, theta, wei);

    //Normalize the weights and compute log-likelihood.

    lik[0]=weight(wei, W, N)-log(N);

    //Save weights and particles

    for(i=0;i<N;i++)
    {
        storeW[i]=wei[i];
    }

    for(i=0;i<N;i++)
    {
        for(k=0;k<dx;k++)
        {
            storeX[dx*i+k]=x[dx*i+k];
        }
    }

    //Compute parameters of the Markov transition

    param=(*paramTransition)(theta,dx,dy);

    //start iterate

    for(iter=1;iter<T;iter++)
    {
        //Resampling

        (*resampling)(random, x, dx, N, N, W, xh);

        //Mutation

        (*simTransition)(random, y, dy, dx, iter, N, param, xh, x);

        //Evaluate the potential function

        (*potential)(y, x, xh, dy, dx, iter, N, theta, wei);

        //Nomalize the weights and compute log-likelihood

        lik[iter]=lik[iter-1]+weight(wei, W, N)-log(N);

        //Save weights and particles

        for(i=0;i<N;i++)
        {
            storeW[N*iter+i]=wei[i];
        }

        for(i=0;i<N;i++)
        {
            for(k=0;k<dx;k++)
            {
                storeX[n*iter+dx*i+k]=x[dx*i+k];
            }
        }

    }

    free(wei);
    wei=NULL;
    free(xh);
    xh=NULL;
    free(param);
    param=NULL;
}


/****************************************************************************************************
SQMC GENERALIZED INFORMATION FILTER
****************************************************************************************************/


void SMC_BIF(double *y, int dy, int dx, int T, double *theta, int N,  gsl_rng *random, ParamTransitionPF paramTransition,
             ResamplingPF resampling, SimInitPF simInit, SimTransitionPF simTransition, PotentialPF potential, int *par, double *lik,  double *storeX,
             double *storeW, double *x, double *W)
{
    int i, j, k, iter, n=N*dx;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *xh=(double*)malloc(sizeof(double)*n);

    double *param=NULL;


    //1.1 Initialization

    //Sample from initial distribution

    (*simInit)(random, dx, dy, N, theta, x);

    //Evaluate the potential function

    (*potential)(y, x, xh, dy, dx, T-1, N, theta, wei);

    //Normalize the weights and compute log-likelihood.

    lik[0]=weight(wei, W, N)-log(N);

    //Save weights and particles

    for(i=0;i<N;i++)
    {
        storeW[N*(T-1)+i]=wei[i];
    }

    for(i=0;i<N;i++)
    {
        for(k=0;k<dx;k++)
        {
            storeX[n*(T-1)+dx*i+k]=x[dx*i+k];
        }
    }

    //Compute parameters of the Markov transition

    param=(*paramTransition)(theta,dx,dy);

    //1.2. Iterate

    for(iter=1;iter<T;iter++)
    {

        (*resampling)(random, x, dx, N, N, W, xh);

        //Mutation

        (*simTransition)(random, y, dy, dx, T-1-iter, N, param, xh, x);

        //Evaluate the potential function

        (*potential)(y, x, xh, dy, dx, T-1-iter, N, theta, wei);

        //1.2.3 Store weights and particles

        lik[iter]=lik[iter-1]+weight(wei, W, N)-log(N);

        //Save weights and particles


        for(i=0;i<N;i++)
        {
            for(k=0;k<dx;k++)
            {
                storeX[n*(T-iter-1)+dx*i+k]=x[dx*i+k];
            }
        }

        for(i=0;i<N;i++)
        {
            storeW[N*(T-iter-1)+i]=wei[i];
        }


    }

    free(wei);
    wei=NULL;
    free(xh);
    xh=NULL;
    free(param);
    param=NULL;

}








