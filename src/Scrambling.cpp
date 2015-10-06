#include "Scrambling.h"
#include <Rcpp.h>

void Scrambling(Scrambled* os )
{
    Scrambled_Randomize(os);
}


void getPoints(Scrambled* os, int d, int n, double *seq)
{
    Scrambled_GetPoints(os, d, n, seq);
}



void Digital_getPoints(DigitalNetGenerator* os, int d, int n, double *seq)
{
    DigitalNetGenerator_GetPoints(os, d, n, seq);
}

/****************************************************************************************************/
// [[Rcpp::export]]
Rcpp::NumericVector RQMC(int qmc, int N, int dx, int ns, Rcpp::NumericVector x)
{
	int i,k;

	Scrambled *os1=Scrambled_Create(qmc, dx, N);
	std::vector<double> x1 = Rcpp::as<std::vector<double> >(x);

	double *work=(double*)malloc(sizeof(double)*(N*(dx)));

	for(i=0;i<ns;i++)
	{
		Scrambling(os1);

		getPoints(os1, dx, N, work);

		for(k=0;k< N*(dx);k++)
		{
			x1[N*(dx)*i+k]=work[k];
		}
	}

	free(work);
	work=NULL;
	Scrambled_Destroy(os1);
	os1=NULL;
	return Rcpp::wrap(x1);
}

// [[Rcpp::export]]
Rcpp::NumericVector QMC (int qmc, int N, int dx,  Rcpp::NumericVector x)
{

    printf("Hello World!\n");


    DigitalNetGenerator *os1=DigitalNetGenerator_Create(qmc, dx);

    std::vector<double> x1 = Rcpp::as<std::vector<double> >(x);
  	DigitalNetGenerator_GetPoints(os1, dx, N, &x1[0]);

	  DigitalNetGenerator_Destroy(os1);
  	os1=NULL;
  	return Rcpp::wrap(x1);

}
























