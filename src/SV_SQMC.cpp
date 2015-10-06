

#include <R.h>


#include "HilbertResamplerQMC.h"
#include "Resampler.h"
#include "functionsCC.hpp"

#include "Model_LG.h"
#include "Model_SV.h"

#include "SMC.h"
#include "SMC_Back.h"
#include "SQMC.h"
#include "SQMC_Back.h"

/*******************************************************************************************************************************************/

//void SQMC_SV(double*, int, int, int, double*, int, int, int, int,  int, int,double*, double*, double*, double*);

//void SQMCBack_SV(double*, int, int, int, double*, int, int, int, int, int, int, int,  double*, double*, double*, double*);

//void SQMC2F_SV(double*, int, int, int, int, double*, double*, int, int, int, int, int, double*, double*, double*, double*);
/************************* PARTICLE FILTER ***************************************************************************************************/

void SQMC_SV(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int qmc, int src, int computeExp,
double *parPsi, double *lik, double *expx, double *expx2)
{
	int i,j,k;
	int par[2]={computeExp, src};

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

	double *storeExp=(double*)malloc(sizeof(double)*(dx*(T)));
	double *v1=(double*)malloc(sizeof(double)*(T));

	for(i=0;i<ns;i++)
	{
		if(qmc!=0)
		{
			SQMC(y, dy, dx, T, theta, N, qmc, paramTransQMC_SV, HilbertResampler, simPriorQMC_SV, simTransitionQMC_LG,potential_SV, &par[0], parPsi, v1,storeExp);
		}
		else
		{
			SMC(y, dy, dx, T, theta, N, random, paramTrans_SV, ResamplePF, simPrior_SV, simTransition_LG, potential_SV, par,v1,storeExp);
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

	free(storeExp);
	storeExp=NULL;
	free(v1);
	v1=NULL;
	gsl_rng_free(random);
	random=NULL;
}

// SQMC_SV
RcppExport SEXP is2_SQMC_SV(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP qmcSEXP, SEXP srcSEXP, SEXP computeExpSEXP, SEXP parPsiSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
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

    SQMC_SV(&y1[0], dy, dx, T, &theta1[0], seed, ns, N, qmc, src, computeExp, &parPsi1[0], &lik1[0], &expx1[0], &expx21[0]);

    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "parPsi" ) = parPsi1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;

    END_RCPP
}


//// [[Rcpp::export]]
void SQMCBack_SV(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int Nb, int qmc, int qmcB,  int Marg, double *parPsi,
double *lik, double *expx, double *expx2)
{
	int i, j, k,  n=(N)*(dx);
	double cc;

	double *W=(double*)malloc(sizeof(double)*(N));
	double *storeW=(double*)malloc(sizeof(double)*(N*(T)));

	double *x=(double*)malloc(sizeof(double)*(n));
	double *storeX=(double*)malloc(sizeof(double)*((T)*n));

	double *xT=(double*)malloc(sizeof(double)*(Nb*(dx)));

	int par[2]={1, 1};

	for(j=0;j<T;j++)
	{
		for(k=0;k<dx;k++)
		{
			expx[dx*j+k]=0;
			expx2[dx*j+k]=0;
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

			SQMC_Forward(y, dy, dx, T, theta, N,  qmc, paramTransQMC_SV, ForHilbertResampler,simPriorQMC_SV,
			simTransitionQMC_LG, potential_SV, par, parPsi, v1, storeX, storeW, x, W);

			if(Marg==0)
			{
				if(qmcB==0)
				{
					SMC_Back(dx, T, 0, theta, N, Nb, random, ResamplePF, transition_SV, x, storeX, W, storeW, storeExp, xT);
				}
				else
				{
					SQMC_Back(dx, T, 0, theta, N, Nb, qmcB, BackHilbertResampler, transition_SV, parPsi, x,
					storeX, W, storeW,  storeExp , xT);
				}
			}else
			{
				if(qmcB==0)
				{
					SMC_BackMarg(dx, T, 0, theta, N, Nb, random, ResamplePF, transition_SV, x, storeX, W, storeW, storeExp, xT);
				}
				else
				{
					SQMC_BackMarg(dx, T, 0, theta, N, Nb, qmcB, BackHilbertResampler, transition_SV, parPsi, x,
					storeX, W, storeW,  storeExp , xT);
				}
			}
		}
		else
		{
			if(Marg==0)
			{
				SMC_ForBack(y, dy, dx, T, theta, N, Nb, random, paramTrans_SV, ResamplePF, simPrior_SV,
				transition_SV,  simTransition_LG, potential_SV, par, v1,storeExp);
			}
			else
			{
				SMC_ForBackMarg(y, dy, dx, T, theta, N, Nb, random, paramTrans_SV, ResamplePF, simPrior_SV,
				transition_SV,  simTransition_LG, potential_SV, par, v1,storeExp);
			}
		}


		for(j=0;j<T;j++)
		{
			lik[T*i+j]=v1[j];
		}

		for(j=0;j<T;j++)
		{
			for(k=0;k<dx;k++)
			{
				expx[dx*j+k]+=storeExp[dx*j+k]/(ns);
				expx2[dx*j+k]+=pow(storeExp[dx*j+k],2)/(ns);
			}
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

RcppExport SEXP is2_SQMCBack_SV(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP NbSEXP, SEXP qmcSEXP, SEXP qmcBSEXP, SEXP MargSEXP, SEXP parPsiSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
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
    Rcpp::traits::input_parameter< int >::type qmcB(qmcBSEXP);
    Rcpp::traits::input_parameter< int >::type Marg(MargSEXP);

    Rcpp::NumericVector parPsi(parPsiSEXP);
    std::vector<double> parPsi1 = Rcpp::as<std::vector<double> >(parPsi);

    Rcpp::NumericVector lik(likSEXP);
    std::vector<double> lik1 = Rcpp::as<std::vector<double> >(lik);

    Rcpp::NumericVector expx(expxSEXP);
    std::vector<double> expx1 = Rcpp::as<std::vector<double> >(expx);

    Rcpp::NumericVector expx2(expx2SEXP);
    std::vector<double> expx21 = Rcpp::as<std::vector<double> >(expx2);

    SQMCBack_SV(&y1[0], dy, dx, T, &theta1[0], seed, ns, N, Nb, qmc, qmcB, Marg, &parPsi1[0], &lik1[0], &expx1[0], &expx21[0]);

    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "parPsi" ) = parPsi1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;
    END_RCPP
}





//// [[Rcpp::export]]
void SQMC2F_SV(double *y, int dy, int dx, int tstar, int T, double *theta, double *theta2F, int seed, int ns, int N, int Nb, int qmc,
double *parPsi, double *lik, double *expx, double *expx2)
{
	int i,j,k;
	int par[2]={1, 1};

	for(j=0;j<T;j++)
	{
		for(k=0;k<dx;k++)
		{
			expx[dx*j+k]=0;
			expx2[dx*j+k]=0;
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
			SQMC_2F(y, tstar, dy, dx, T, theta, theta2F, N, Nb, qmc, paramTransQMC_SV, ForHilbertResampler,
			BackHilbertResampler, simPriorQMC_SV, transition_SV, simTransitionQMC_LG,
			potential_SV, paramTransQMC_SV,  simPriorQMC_SV, transition_SV, simTransitionQMC_LG,
			potential_SV,  par, parPsi, v1, storeExp);
		}
		else
		{
			SMC_2F(y, tstar, dy, dx, T, theta, theta2F, N, Nb,  random, paramTrans_SV, ResamplePF, simPrior_SV,
			transition_SV, simTransition_LG, potential_SV, paramTrans_SV, simPrior_SV, transition_SV,
			simTransition_LG, potential_SV,  par, v1, storeExp);
		}

		for(j=0;j<T;j++)
		{
			lik[T*i+j]=v1[j];
		}

		for(j=0;j<T;j++)
		{
			for(k=0;k<dx;k++)
			{
				expx[dx*j+k]+=storeExp[dx*j+k]/(ns);
				expx2[dx*j+k]+=pow(storeExp[dx*j+k],2)/(ns);
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

// SQMC2F_SV
RcppExport SEXP is2_SQMC2F_SV(SEXP ySEXP, SEXP dySEXP, SEXP dxSEXP, SEXP tstarSEXP, SEXP TSEXP, SEXP thetaSEXP, SEXP theta2FSEXP, SEXP seedSEXP, SEXP nsSEXP, SEXP NSEXP, SEXP NbSEXP, SEXP qmcSEXP, SEXP parPsiSEXP, SEXP likSEXP, SEXP expxSEXP, SEXP expx2SEXP) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;

    Rcpp::NumericVector y(ySEXP);
    std::vector<double> y1 = Rcpp::as<std::vector<double> >(y);

    Rcpp::traits::input_parameter< int >::type dy(dySEXP);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type tstar(tstarSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);


    Rcpp::NumericVector theta(thetaSEXP);
    std::vector<double> theta1 = Rcpp::as<std::vector<double> >(theta);

    Rcpp::NumericVector theta2F(theta2FSEXP);
    std::vector<double> theta2F1 = Rcpp::as<std::vector<double> >(theta2F);

    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Nb(NbSEXP);
    Rcpp::traits::input_parameter< int >::type qmc(qmcSEXP);


    Rcpp::NumericVector parPsi(parPsiSEXP);
    std::vector<double> parPsi1 = Rcpp::as<std::vector<double> >(parPsi);

    Rcpp::NumericVector lik(likSEXP);
    std::vector<double> lik1 = Rcpp::as<std::vector<double> >(lik);

    Rcpp::NumericVector expx(expxSEXP);
    std::vector<double> expx1 = Rcpp::as<std::vector<double> >(expx);

    Rcpp::NumericVector expx2(expx2SEXP);
    std::vector<double> expx21 = Rcpp::as<std::vector<double> >(expx2);

    SQMC2F_SV(&y1[0], dy, dx, tstar, T, &theta1[0], &theta2F1[0], seed, ns, N, Nb, qmc, &parPsi1[0], &lik1[0], &expx1[0], &expx21[0]);

    return Rcpp::List::create(Rcpp::Named( "y" ) = y1,
                              Rcpp::Named( "theta" ) = theta1,
                              Rcpp::Named( "parPsi" ) = parPsi1,
                              Rcpp::Named( "lik" ) = lik1,
                              Rcpp::Named( "expx" ) = expx1,
                              Rcpp::Named( "expx2" ) = expx21) ;
    END_RCPP
}












