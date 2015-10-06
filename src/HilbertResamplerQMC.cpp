
#include <Rcpp.h>
using namespace Rcpp;
#include "HilbertResamplerQMC.h"
#include "qsort_index.h"



/****************************************************/


void HilbertResampler(double *quasi, double *parPhi, double *x, int dx, int N, double *W, double* xh)
{
    int i,k;
    uint64_t *h=NULL;
    int *J1=(int*)malloc(sizeof(int)*N);

    //FILE* fw=NULL;

    if(dx==1)
    {
        J1=qsort_index(x, N, sizeof(double), compare_double);
    }
    else
    {
        h=(uint64_t*)malloc(sizeof(uint64_t)*N);

        HilbertIndex(x, parPhi, &dx, &N , h);

        J1=qsort_index(h, N, sizeof(uint64_t), compare_uint64);


        /*fw = fopen("hval", "w");


         for(i=0;i<N;i++)
         {
         fprintf(fw, "%f", (double)(h[J1[i]]));
         fprintf(fw, "\n ");

         }


         fprintf(fw, "),\n nrow=%d)\n\n", dx);
         fclose(fw);*/

        free(h);
        h=NULL;
    }

    quasi_Resample(quasi, &N, W, J1, x, xh);

    free(J1);
    J1=NULL;
}
RcppExport SEXP is2_HilbertResampler(SEXP quasiSEXP, SEXP parPhiSEXP, SEXP xSEXP, SEXP dxSEXP, SEXP NSEXP, SEXP WSEXP, SEXP xhSEXP) {
    BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::NumericVector quasi(quasiSEXP);
    std::vector<double> quasi1 = Rcpp::as<std::vector<double> >(quasi);
    Rcpp::NumericVector parPhi(parPhiSEXP);
    std::vector<double> parPhi1 = Rcpp::as<std::vector<double> >(parPhi);
    Rcpp::NumericVector x(xSEXP);
    std::vector<double> x1 = Rcpp::as<std::vector<double> >(x);
    Rcpp::traits::input_parameter< int >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::NumericVector W(WSEXP);
    std::vector<double> W1 = Rcpp::as<std::vector<double> >(W);
    Rcpp::NumericVector xh(xhSEXP);
    std::vector<double> xh1 = Rcpp::as<std::vector<double> >(xh);
    HilbertResampler(&quasi1[0], &parPhi1[0], &x1[0], dx, N, &W1[0], &xh1[0]);
    return Rcpp::List::create(Rcpp::Named( "quasi" ) = quasi1,
                              Rcpp::Named( "parPhi" ) = parPhi1,
                              Rcpp::Named( "x" ) = x1,
                              Rcpp::Named( "W" ) = W1,
                              Rcpp::Named( "xh" ) = xh1) ;
    END_RCPP
}


/****************************************************/

void quasi_Resample(double *quasi,  int *N, double *W, int *J1, double *x, double *J)
{
	int i,j,k;
	double s;

    	s=0;
	j=0;

	for(i=0; i<*N; i++)
	{
		s+=W[J1[i]];

		while(quasi[j]<=s && j<*N)
		{
			J[j]=J1[i]+1;

			j++;

		}
		if(j==*N)
		{
			break;
		}

	}
}

int* ForHilbertResampler(double *quasi, int dimQuasi, double *parPhi, double *x, int dx, int N, int Nb, double *W, double* xh)
{
    int i,k;
    uint64_t *h=NULL;
    int *J1=(int*)malloc(sizeof(int)*N);

    //FILE* fw=NULL;

    if(dx==1)
    {
        J1=qsort_index(x, N, sizeof(double), compare_double);
    }
    else
    {
        h=(uint64_t*)malloc(sizeof(uint64_t)*N);

        HilbertIndex(x, parPhi, &dx, &N , h);

        J1=qsort_index(h, N, sizeof(uint64_t), compare_uint64);


        /*fw = fopen("hval", "w");


         for(i=0;i<N;i++)
         {
         fprintf(fw, "%f", (double)(h[J1[i]]));
         fprintf(fw, "\n ");

         }


         fprintf(fw, "),\n nrow=%d)\n\n", dx);
         fclose(fw);*/

        free(h);
        h=NULL;
    }
    quasi_ResampleFor(quasi, dimQuasi, dx, N, Nb, W, J1, x, xh);

    return(J1);

}




void quasi_ResampleFor(double *quasi, int dimQuasi, int dx,  int N, int Nb, double *W, int *J1, double *x, double *xh)
{
    int i,j,k;
    double s;

    s=0;
    j=0;

    for(i=0; i<N; i++)
    {
        s+=W[J1[i]];

        while(quasi[dimQuasi*j]<=s && j<Nb)
        {
            for(k=0;k<dx;k++)
            {
                xh[dx*j+k]=x[dx*J1[i]+k];
            }

            j++;
        }
        if(j==Nb)
        {
            break;
        }

    }
}

/*****************************************************************************************************************
RESAMPLER ALGORITHMS FOR  BACKWARD STEP
******************************************************************************************************************/


int* BackHilbertResampler(double *quasi, int dimQuasi, double *parPhi, double *x, int dx, int N, int Nb, double *W, double* xh)
{
    int i;

    uint64_t *h=NULL;

    int *J1=(int*)malloc(sizeof(int)*N);
    int *J2=(int*)malloc(sizeof(int)*Nb);

    if(dx==1)
    {
        J1=qsort_index(x, N, sizeof(double), compare_double);
    }
    else
    {
        h=(uint64_t*)malloc(sizeof(uint64_t)*N);

        HilbertIndex(x, parPhi, &dx, &N , h);

        J1=qsort_index(h, N, sizeof(uint64_t), compare_uint64);


        /*fw = fopen("hval", "w");


         for(i=0;i<N;i++)
         {
         fprintf(fw, "%f", (double)(h[J1[i]]));
         fprintf(fw, "\n ");

         }


         fprintf(fw, "),\n nrow=%d)\n\n", dx);
         fclose(fw);*/

        free(h);
        h=NULL;
    }


    quasi_ResampleBack(quasi, dimQuasi, dx, N, Nb, W, J1, J2, x, xh);

    free(J1);
    J1=NULL;
    return(J2);

}


void quasi_ResampleBack(double *quasi, int dimQuasi, int dx,  int N, int Nb, double *W, int *J1, int *J2, double *x, double *xh)
{
    int i,j,k;
    double s;

    s=0;
    j=0;

    for(i=0; i<N; i++)
    {
        s+=W[J1[i]];

        while(quasi[dimQuasi*j]<=s && j<Nb)
        {
            for(k=0;k<dx;k++)
            {
                xh[dx*j+k]=x[dx*J1[i]+k];
            }
            J2[j]=J1[i];
            j++;
        }
        if(j==Nb)
        {
            break;
        }

    }
}






