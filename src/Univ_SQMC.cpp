
#include <R.h>


#include "HilbertResamplerQMC.h"
#include "Resampler.h"
#include "functionsCC.hpp"

#include "Model_Univ.h"

#include "SQMC.h"
#include "SQMC_Back.h"
#include "SMC.h"
#include "SMC_Back.h"

void SQMC_Univ(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int qmc, int src, int computeExp,
double *lik, double *expx, double *expx2)
{
	int i,j;
	int par[2]={computeExp, src};
	double *parPsi=NULL;

	if(par[0]==1)
	{
		for(j=0;j<T;j++)
		{
			expx[j]=0;
			expx2[j]=0;
		}
	}

	gsl_rng *random=gsl_rng_alloc (gsl_rng_mt19937);

	if(seed>0)
	{
		gsl_rng_set(random, seed);
	}
	else
	{
		gsl_rng_set(random, getSeed());
	}

	double *storeExp=(double*)malloc(sizeof(double)*(T));
	double *v1=(double*)malloc(sizeof(double)*(T));

	for(i=0;i<ns;i++)
	{
		if(qmc!=0)
		{
			SQMC(y, dy, dx, T, theta, N, qmc, paramTrans_Univ, HilbertResampler,
			simPriorQMC_Univ, simTransitionQMC_Univ, potential_Univ, par,parPsi, v1,storeExp);
		}
		else
		{
			SMC(y, dy, dx, T, theta, N, random, paramTrans_Univ, ResamplePF, simPrior_Univ, simTransition_Univ,
			potential_Univ, par,v1,storeExp);
		}

		for(j=0;j<T;j++)
		{
			lik[T*i+j]=v1[j];
		}

		if(par[0]==1)
		{
			for(j=0;j<T;j++)
			{
				expx[j]+=storeExp[j]/(ns);
				expx2[j]+=pow(storeExp[j],2)/(ns);
			}
		}
	}

	free(storeExp);
	storeExp=NULL;
	free(v1);
	v1=NULL;
	gsl_rng_free(random);
	random=NULL;
}

// SQMC_Univ

RcppExport SEXP is2_SQMC_Univ(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP qmcSEXP, SEXP srcSEXP, SEXP computeExpSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;

    Rcpp::NumericVector y(ySEXP);
    std::vector<double> y1 = Rcpp::as<std::vector<double> >(y);

    Rcpp::traits::input_parameter< int >::type dy(dySEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);

    Rcpp::NumericVector theta(thetaSEXP);
    std::vector<double> theta1 = Rcpp::as<std::vector<double> >(theta);

    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type qmc(qmcSEXP);
    Rcpp::traits::input_parameter< int >::type src(srcSEXP);
    Rcpp::traits::input_parameter< int >::type computeExp(computeExpSEXP);

    Rcpp::NumericVector lik(likSEXP);
    std::vector<double> lik1 = Rcpp::as<std::vector<double> >(lik);

    Rcpp::NumericVector expx(expxSEXP);
    std::vector<double> expx1 = Rcpp::as<std::vector<double> >(expx);

    Rcpp::NumericVector expx2(expx2SEXP);
    std::vector<double> expx21 = Rcpp::as<std::vector<double> >(expx2);

    SQMC_Univ(&y1[0], dy, dx, T, &theta1[0], seed, ns, N, qmc, src, computeExp, &lik1[0], &expx1[0], &expx21[0]);

    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;
    END_RCPP
}


void SQMCBack_Univ(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int Nb, int qmc, int qmcb,  int M, double *lik,
double *expx, double *expx2)
{
	int i, j, k;
	double cc, *parPsi=NULL;

	double *W=(double*)malloc(sizeof(double)*(N));
	double *storeW=(double*)malloc(sizeof(double)*(N*(T)));

	double *x=(double*)malloc(sizeof(double)*(N));
	double *storeX=(double*)malloc(sizeof(double)*((T)*(N)));

	double *xT=(double*)malloc(sizeof(double)*(Nb));

	int par[2]={1, 1};

	if(ns!=1)
	{
		for(j=0;j<T;j++)
		{
			for(k=0;k<dx;k++)
			{
				expx[dx*j+k]=0;
				expx2[dx*j+k]=0;
			}
		}
	}

	gsl_rng *random=gsl_rng_alloc (gsl_rng_mt19937);

	if(seed>0)
	{
		gsl_rng_set(random, seed);
	}
	else
	{
		gsl_rng_set(random, getSeed());
	}

	double *storeExp=(double*)malloc(sizeof(double)*(dx*(T)));
	double *v1=(double*)malloc(sizeof(double)*(T));

	for(i=0;i<ns;i++)
	{
		if(qmc!=0)
		{
			SQMC_Forward(y, dy, dx, T, theta, N,  qmc, paramTrans_Univ, ForHilbertResampler,simPriorQMC_Univ,
			simTransitionQMC_Univ, potential_Univ, par, parPsi, v1, storeX, storeW, x, W);

			if(M==0)
			{
				if(qmcb==0)
				{
					SMC_Back(dx, T, 0, theta, N, Nb, random, ResamplePF, transition_Univ, x, storeX, W, storeW, storeExp, xT);
				}
				else
				{
					SQMC_Back(dx, T, 0, theta, N, Nb, qmcb, BackHilbertResampler, transition_Univ, parPsi, x,
					storeX, W, storeW,  storeExp , xT);
				}
			}
			else
			{
				if(qmcb==0)
				{
					SMC_BackMarg(dx, T, 0, theta, N, Nb, random, ResamplePF, transition_Univ, x, storeX, W, storeW, storeExp, xT);
				}
				else
				{
					SQMC_BackMarg(dx, T, 0, theta, N, Nb, qmcb, BackHilbertResampler, transition_Univ, parPsi, x,
					storeX, W, storeW,  storeExp , xT);
				}
			}
		}
		else
		{
			if(M==0)
			{
				SMC_ForBack(y, dy, dx, T, theta, N, Nb, random, paramTrans_Univ, ResamplePF, simPrior_Univ,
				transition_Univ,  simTransition_Univ, potential_Univ, par, v1,storeExp);
			}
			else
			{
				SMC_ForBackMarg(y, dy, dx, T, theta, N, Nb, random, paramTrans_Univ, ResamplePF, simPrior_Univ,
				transition_Univ,  simTransition_Univ, potential_Univ, par, v1,storeExp);
			}
		}

		for(j=0;j<T;j++)
		{
			lik[T*i+j]=v1[j];
		}

		for(j=0;j<T;j++)
		{
			expx[j]+=storeExp[j]/(ns);
			expx2[j]+=pow(storeExp[j],2)/(ns);
		}
	}

	free(W);
	W=NULL;
	free(storeW);
	storeW=NULL;
	free(x);
	x=NULL;
	free(storeX);
	storeX=NULL;
	free(xT);
	xT=NULL;

	free(storeExp);
	storeExp=NULL;
	free(v1);
	v1=NULL;
	gsl_rng_free(random);
	random=NULL;


}
// SQMCBack_Univ

RcppExport SEXP is2_SQMCBack_Univ(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP NbSEXP, SEXP qmcSEXP, SEXP qmcbSEXP, SEXP MSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;

    Rcpp::NumericVector y(ySEXP);
    std::vector<double> y1 = Rcpp::as<std::vector<double> >(y);

    Rcpp::traits::input_parameter< int >::type dy(dySEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);

    Rcpp::NumericVector theta(thetaSEXP);
    std::vector<double> theta1 = Rcpp::as<std::vector<double> >(theta);


    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Nb(NbSEXP);
    Rcpp::traits::input_parameter< int >::type qmc(qmcSEXP);
    Rcpp::traits::input_parameter< int >::type qmcb(qmcbSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);


    Rcpp::NumericVector lik(likSEXP);
    std::vector<double> lik1 = Rcpp::as<std::vector<double> >(lik);


    Rcpp::NumericVector expx(expxSEXP);
    std::vector<double> expx1 = Rcpp::as<std::vector<double> >(expx);


    Rcpp::NumericVector expx2(expx2SEXP);
    std::vector<double> expx21 = Rcpp::as<std::vector<double> >(expx2);


    SQMCBack_Univ(&y1[0], dy, dx, T, &theta1[0], seed, ns, N, Nb, qmc, qmcb, M, &lik1[0], &expx1[0], &expx21[0]);


    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;

    END_RCPP
}




