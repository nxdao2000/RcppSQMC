#include <Rcpp.h>
#include <R.h>
#include "generate_RQMC.hpp"
#include "SQMC.h"



/****************************************************************************************************
				SEQUENTIAL quasi-MONTE CARLO

Output: Estimate of the log-likelihood function at each time steps ans filtering expectations
******************************************************************************************************/


void SQMC(double *y, int dy, int dx, int T, double *theta, int N, int qmc, ParamTransitionPF paramTransition, ResamplingPF_QMC resampling,
          SimInitPF_QMC simInit,  SimTransitionPF_QMC simTransition, PotentialPF potential, int* par, double *parPsi, double *lik, double *expx)
{
    int i,j,k,iter;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *W=(double*)malloc(sizeof(double)*N);
    double *x=(double*)malloc(sizeof(double)*(N*dx));
    double *xh=(double*)malloc(sizeof(double)*(N*dx));
    double *sim=(double*)malloc(sizeof(double)*(N*(dx+1)));
    double *param=NULL;

    Scrambled * os1=NULL;
    DigitalNetGenerator *os2=NULL;

    //Generate a (R)QMC point set on (0,1)^d

    if(par[1]==1)
    {
        os1 = Scrambled_Create(qmc, dx+1,N);

        Scrambling(os1);

        getPoints(os1, dx, N, sim);
    }
    else
    {
        os2=DigitalNetGenerator_Create(qmc, dx+1);
        DigitalNetGenerator_GetPoints(os2, dx, N, sim);
    }

    //Sample from initial distribution

    (*simInit)(sim, dx, dy, N, theta, x);

    //Evaluate the potential function

    (*potential)(y, x, xh, dy, dx, 0, N, theta, wei);

    //Normalize the weights and compute log-likelihood

    lik[0]=weight(wei, W, N)-log(N);

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

    //Generate QMC point set in (0,1)^d if deterministiv SQMC is used

    if(par[1]==0)
    {
        DigitalNetGenerator_GetPoints(os2, dx+1, N, sim);
    }


    for(iter=1;iter<T;iter++)
    {
        //Generate RQMC points in (0,1)^{d+1}

        if(par[1]==1)
        {
            Scrambling(os1);

            getPoints(os1, dx+1, N, sim);
        }

        //Hilbert Resampling

        (*resampling)(sim, parPsi, x, dx, N, W, xh);

        //Mutation

        //(*simTransition)(sim, y, dy, dx, iter, N, param, xh, x);

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

    if(par[1]==1)
    {
        Scrambled_Destroy(os1);
        os1=NULL;
    }
    else
    {
        DigitalNetGenerator_Destroy(os2);
        os2=NULL;

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
    free(sim);
    sim=NULL;



}

/****************************************************************************************************
FORWARD STEP FOR THE SQMC FORWARD-BACKWARD SMOOTHING ALGORITHM
****************************************************************************************************/


void SQMC_Forward(double *y, int dy, int dx, int T, double *theta, int N,  int qmc, ParamTransitionPF paramTransition,
                  ResamplingBack_QMC resampling,SimInitPF_QMC simInit, SimTransitionPF_QMC simTransition, PotentialPF potential,
                  int* par, double *parPsi, double *lik, double *storeX, double *storeW, double *x, double *W)
{
    int i, j, k, iter, n=N*dx, nx=dx+1;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *xh=(double*)malloc(sizeof(double)*n);
    double *sim=(double*)malloc(sizeof(double)*(N*nx));;

    int *J=NULL;
    double *param=NULL;

    Scrambled * os1=NULL;

    //Generate a (R)QMC point set on (0,1)^d

    os1 = Scrambled_Create(qmc, dx+1,N);

    Scrambling(os1);

    getPoints(os1, dx, N, sim);


    //Sample from initial distribution

    (*simInit)(sim, dx, dy, N, theta, x);

    //Evaluate the potential function

    (*potential)(y, x, xh, dy, dx, 0, N, theta, wei);

    //Normalize the weights and compute the log-liklihood function

    lik[0]=weight(wei, W, N)-log(N);

    //Compute parameters for the MArkov transition

    param=(*paramTransition)(theta,dx,dy);

    for(iter=1;iter<T;iter++)
    {
        //Generate RQMC points in (0,1)^{d+1}

        Scrambling(os1);

        getPoints(os1, dx+1, N, sim);

        //Hilbert Resampling

        J=(*resampling)(sim, dx+1, parPsi, x, dx, N, N, W, xh);

        //Store weights and particles at time iter-1, sorted according to their Hilbert index

        for(i=0;i<N;i++)
        {
            for(k=0;k<dx;k++)
            {
                storeX[n*(iter-1)+dx*i+k]=x[dx*J[i]+k];
            }
        }

        for(i=0;i<N;i++)
        {
            storeW[N*(iter-1)+i]=wei[J[i]];
        }

        //Mutation

        (*simTransition)(sim, y, dy, dx, iter, N, param, xh, x);

        //Evaluate the potential function

        (*potential)(y, x, xh, dy, dx, iter, N, theta, wei);

        //Normalize the weights and compute the log-liklihood function

        lik[iter]=lik[iter-1]+weight(wei, W, N)-log(N);

        free(J);
        J=NULL;

    }

    //Store weights and particles at time T-1, sorted according to their Hilbert index

    J=(*resampling)(sim, nx, parPsi, x, dx, N, N, W, xh);

    for(i=0;i<N;i++)
    {
        for(k=0;k<dx;k++)
        {
            storeX[n*(T-1)+dx*i+k]=x[dx*J[i]+k];
        }
    }

    for(i=0;i<N;i++)
    {
        storeW[N*(T-1)+i]=wei[J[i]];
    }
    free(J);
    J=NULL;

    Scrambled_Destroy(os1);
    os1=NULL;

    free(sim);
    sim=NULL;
    free(xh);
    xh=NULL;
    free(wei);
    wei=NULL;
    free(param);
    param=NULL;

}

/****************************************************************************************************
SQMC GENERALIZED INFORMATION FILTER
****************************************************************************************************/



void SQMC_BIF(double *y, int dy, int dx, int T, double *theta, int N,  int qmc, ParamTransitionPF paramTransition,
              ResamplingBack_QMC resampling,SimInitPF_QMC simInit, SimTransitionPF_QMC simTransition, PotentialPF potential,
              int* par, double *parPsi, double *lik, double *storeX, double *storeW, double *x, double *W)
{
    int i, j, k, iter, n=N*dx, nx=dx+1;
    double cc;

    double *wei=(double*)malloc(sizeof(double)*N);
    double *xh=(double*)malloc(sizeof(double)*n);
    double *sim=(double*)malloc(sizeof(double)*(N*nx));;

    int *J=NULL;
    double *param=NULL;

    Scrambled * os1=NULL;


    //1.1 Initialization

    //1.1.1 generate RQMC points in (0,1)^d

    os1 = Scrambled_Create(qmc, dx+1,N);

    Scrambling(os1);

    getPoints(os1, dx, N, sim);


    //1.1.2 sample from the prior distribution

    (*simInit)(sim, dx, dy, N, theta, x);

    //1.1.3 correction

    (*potential)(y, x, xh, dy, dx, T-1, N, theta, wei);

    cc=weight(wei, W, N);

    param=(*paramTransition)(theta,dx,dy);

    //1.2. Iterate

    for(iter=1;iter<T;iter++)
    {
        //1.2.1 generate RQMC points in (0,1)^{d+1}

        Scrambling(os1);

        getPoints(os1, dx, N, sim);

        //1.2.2 Hilbert Resampling

        J=(*resampling)(sim, nx, parPsi, x, dx, N, N, W, xh);

        //1.2.3 Store weights and particles at time iter-1, sorted according to their Hilbert index

        for(i=0;i<N;i++)
        {
            for(k=0;k<dx;k++)
            {
                storeX[n*(T-iter)+dx*i+k]=x[dx*J[i]+k];
            }
        }

        for(i=0;i<N;i++)
        {
            storeW[N*(T-iter)+i]=wei[J[i]];
        }

        //1.2.4 Mutation

        (*simTransition)(sim, y, dy, dx, T-1-iter, N, param, xh, x);

        //1.2.5 correction

        (*potential)(y, x, xh, dy, dx, T-1-iter, N, theta, wei);

        cc=weight(wei, W, N);

        free(J);
        J=NULL;

    }

    //Store weights and particles at time 0, sorted according to their Hilbert index

    J=(*resampling)(sim, nx, parPsi, x, dx, N, N, W, xh);

    for(i=0;i<N;i++)
    {
        for(k=0;k<dx;k++)
        {
            storeX[dx*i+k]=x[dx*J[i]+k];
        }
    }

    for(i=0;i<N;i++)
    {
        storeW[i]=wei[J[i]];
    }
    free(J);
    J=NULL;

    Scrambled_Destroy(os1);
    os1=NULL;
    free(sim);
    sim=NULL;
    free(xh);
    xh=NULL;
    free(wei);
    wei=NULL;
    free(param);
    param=NULL;

}




