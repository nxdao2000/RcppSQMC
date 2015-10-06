


#include "Resampler.h"


//systematic resampling

void ResamplePF(gsl_rng *r, double *x, int dx, int N, int Nb, double *prob, double *xh)
{
	int i,j,k;
	double s;
	double u=(double)gsl_rng_uniform(r)/(Nb);
	
    	s=0;
	j=0;
	for(i=0;i<N;i++)
	{
		s+=prob[i];
		while(u<=s && j<Nb)
		{
			for(k=0;k<dx;k++)
			{
				xh[dx*j+k]=x[dx*i+k];
			}
			j++;
			u+=(double)1/(Nb);	
		}

		if(j==Nb)
		{
			break;
		}
	}

}


int* ResamplePF_index(gsl_rng *r, double *x, int dx, int N, int Nb, double *prob, double *xh)
{
	int i,j,k;
	double s;
	int *J=(int*)malloc(sizeof(int)*Nb);

	double u=(double)gsl_rng_uniform(r)/(Nb);
	
    	s=0;
	j=0;
	for(i=0;i<N;i++)
	{
		s+=prob[i];
		while(u<=s && j<Nb)
		{
			for(k=0;k<dx;k++)
			{
				xh[dx*j+k]=x[dx*i+k];
			}
			J[j]=i;
			j++;
			u+=(double)1/(Nb);	
		}

		if(j==Nb)
		{
			break;
		}
	}
	return(J);

}


//Resample on point 

void ResampleBack(double u, double *x, int dx, int N, double *prob, int index, double *xh)
{
	int i,j,k;
	double s;
	
    	s=0;
	j=0;
	for(i=0;i<N;i++)
	{
		s+=prob[i];
		if(u<=s && j<1)
		{
			for(k=0;k<dx;k++)
			{
				xh[dx*index+k]=x[dx*i+k];
			}
			j++;
		}
		if(j==1)
		{
			break;
		}

	}
	

}



//Resample Nb points  (use for QMC marginal  Backward smoothing)

void ResampleMargBackQMC(double *quasi, double *x, int dx, int N, int Nb, double *prob,  double *Wb, double *xh)
{
	int i,j,k;
	double s;
	
	qsort(quasi, Nb, sizeof(double), compare_double); 

    	s=0;
	j=0;
	for(i=0;i<N;i++)
	{
		s+=prob[i];
		while(quasi[j]<=s && j<Nb)
		{
			for(k=0;k<dx;k++)
			{
				xh[dx*j+k]=x[dx*i+k];
			}
			j++;
		}
		if(j==Nb)
		{
			break;
		}
	}


}


//Resample Nb points  (use for MC marginal  Backward smoothing)

void ResampleMargBackMC(gsl_rng *r, double *x, int dx, int N, int Nb, double *prob,  double *Wb, double *xh)
{
	int i,j,k;
	double s;
	
	double u=(double)(gsl_rng_uniform(r)/(Nb));
	
    	s=0;
	j=0;
	for(i=0;i<N;i++)
	{
		s+=prob[i];
		while(u<=s && j<Nb)
		{
			for(k=0;k<dx;k++)
			{
				xh[dx*j+k]=x[dx*i+k];
			}
			u+=(double)(1.0/(Nb));
			j++;
		}
		if(j==Nb)
		{
			break;
		}
	}


}



	



	

