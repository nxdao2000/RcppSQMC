


#include <R.h>


#include "HilbertResamplerQMC.h"
#include "Resampler.h"
#include "functionsCC.hpp"

#include "Model_LG.h"
#include "Model_Neuro.h"

#include "SQMC.h"
#include "SMC.h"




void SQMC_Neuro(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int qmc, int src, int computeExp,
double *parPsi, double *lik, double *expx, double *expx2)
{
	int i,j,k;
	int par[2]={computeExp, src};

	if(par[0]==1 && ns!=1)
	{
		for(j=0;j<T;j++)
		{
			expx[j]=0;
			expx2[j]=0;
		}
	}

	gsl_rng *random=gsl_rng_alloc (gsl_rng_taus);

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
		if(ns==1)
		{
			if(qmc!=0)
			{
				SQMC(y, dy, dx, T, theta, N, qmc, paramTransQMC_Neuro, HilbertResampler,
				simPriorQMC_LG, simTransitionQMC_Neuro,potential_Neuro, par, parPsi, lik,expx);
			}
			else
			{
				SMC(y, dy, dx, T, theta, N, random, paramTrans_Neuro, ResamplePF, simPrior_LG,
				simTransition_Neuro, potential_Neuro, par,lik,expx);
			}
		}
		else
		{
			if(qmc!=0)
			{
				SQMC(y, dy, dx, T, theta, N, qmc, paramTransQMC_Neuro, HilbertResampler,
				simPriorQMC_LG, simTransitionQMC_Neuro,potential_Neuro, par, parPsi, v1,storeExp);
			}
			else
			{
				SMC(y, dy, dx, T, theta, N, random, paramTrans_Neuro, ResamplePF, simPrior_LG,
				simTransition_Neuro, potential_Neuro, par,v1,storeExp);
			}

			for(j=0;j<T;j++)
			{
				lik[T*i+j]=v1[j];
			}

			if(par[0]==1)
			{
				for(j=0;j<T;j++)
				{
					for(k=0;k<dx;k++)
					{
						expx[dx*j+k]+=storeExp[dx*j+k]/(ns);
						expx2[dx*j+k]+=pow(storeExp[dx*j+k],2)/(ns);
					}
				}
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

// SQMC_Neuro
RcppExport SEXP is2_SQMC_Neuro(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP qmcSEXP, SEXP srcSEXP, SEXP computeExpSEXP, SEXP parPsiSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
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


    Rcpp::NumericVector parPsi(parPsiSEXP);
    std::vector<double> parPsi1 = Rcpp::as<std::vector<double> >(parPsi);


    Rcpp::NumericVector lik(likSEXP);
    std::vector<double> lik1 = Rcpp::as<std::vector<double> >(lik);


    Rcpp::NumericVector expx(expxSEXP);
    std::vector<double> expx1 = Rcpp::as<std::vector<double> >(expx);


    Rcpp::NumericVector expx2(expx2SEXP);
    std::vector<double> expx21 = Rcpp::as<std::vector<double> >(expx2);


    SQMC_Neuro(&y1[0], dy, dx, T, &theta1[0], seed, ns, N, qmc, src, computeExp, &parPsi1[0], &lik1[0], &expx1[0], &expx21[0]);


    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "parPsi" ) = parPsi1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;
    END_RCPP
}














