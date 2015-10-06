

#include "Model_Univ.h"

/*****************************************************************************************************************************/


void simPrior_Univ(gsl_rng* random, int dx, int dy, int N, double *theta, double *res)
{
	int i=0;

	for(i=0;i<N; i++)
	{
		res[i]=theta[6]+gsl_ran_gaussian(random,theta[7]);
	}
}


void simPriorQMC_Univ(double* qseq, int dx, int dy, int N, double *theta, double *res)
{
	int i;
	for(i=0;i<N;i++)
	{
		res[i]=theta[6]+theta[7]*gsl_cdf_ugaussian_Pinv(qseq[i]);
	}
}



double *paramTrans_Univ(double *theta, int dx, int dy)
{
	int k;
	double *param=(double*)malloc(sizeof(double)*(8));
	
	for(k=0;k< 8;k++)
	{
		param[k]=theta[k];
	}
	return(param);
}


void  simTransition_Univ(gsl_rng *random, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i=0;
	for(i=0;i<N; i++)
	{
		x[i]=param[0]*xh[i]+param[1]*(double)xh[i]/(1+pow(xh[i],2))+param[2]*cos(param[3]*time)+gsl_ran_gaussian(random,param[4]);
	}
}	


void simTransitionQMC_Univ(double *qseq, double *y, int dy, int dx, int time, int N, double *param, double *xh, double *x)
{
	int i=0;

	for(i=0;i<N; i++)
	{
		x[i]=param[0]*xh[i]+param[1]*(double)xh[i]/(1+pow(xh[i],2))+param[2]*cos(param[3]*time)+
		param[4]*gsl_cdf_ugaussian_Pinv(qseq[(dx+1)*i+1]);
	}
}





/************************MEASURE********************************************************************/

void potential_Univ(double *y, double *x, double *xh, int dy,  int dx, int time, int N, double *theta, double *res)
{
	int i=0;
	for(i=0;i<N;i++)
	{
		res[i]=-0.5*log(2*M_PI)-0.5*pow(y[time]-(double)(pow(x[i],2)/theta[5]),2);
	}
}	

void transition_Univ(double *x, double *xh, int dx, int time, int N, double *theta, double *wei)
{
	int i;
	double work;

	for(i=0;i<N; i++)
	{
		work=theta[0]*xh[i]+theta[1]*(double)xh[i]/(1+pow(xh[i],2))+theta[2]*cos(theta[3]*time);
		wei[i]=-0.5*pow((x[0]-work)/theta[4],2);
	}
}
































