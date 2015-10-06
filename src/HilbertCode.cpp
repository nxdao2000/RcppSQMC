
#include"Hilbert.hpp"

#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>


#define PRES 64	//MAXIMUM ADMITED VALUE IS 64. CAN BE REDUCED TO ACCELERATE THE HILBERT SORT



typedef CFixBitVec FV;
typedef CBigBitVec BV;



#include"HilbertCode.hpp"



//// [[Rcpp::export]]
void HilbertIndex(double *x, double *par, int *dx, int *N, uint64_t *d)
{
	int i,k, count;
	double val1, val2,factor,xx;
	int u=floor(PRES/(*dx));
	int dd=floor((u-1)*log(2)/log(10));

	double *work = new double[*N*(*dx)];

	val1=-log(1+exp(-(x[0]-par[0])/(par[1]-par[0])));

	for(i=0;i<*N;i++)
	{
		for(k=0;k<*dx;k++)
		{
			work[*dx*i+k]=-log(1+exp(-(x[*dx*i+k]-par[2*k])/(par[2*k+1]-par[2*k])));

			if(work[*dx*i+k]>val1);
			{
				val1=work[*dx*i+k];
			}
		}
	}

	val2=exp(val1);
	count = 0;

	while(static_cast<long>(val2) % 10 == 0)
	{
    		val2*= 10;
    		++count;
	}

	factor=(dd+count)*log(10);

	BV *p;
	BV h(PRES);

	p= new BV[*dx];

	for(i=0;i<*N;i++)
	{
		for(k=0;k<*dx;k++)
		{
			p[k]=(uint32_t)exp(work[*dx*i+k]+factor);
		}

		Hilbert::CoordsToIndex<BV,BV>(p, u, *dx, h);
		d[i]=h.Rack();
	}
	delete [] p;
	p=NULL;
	delete [] work;
	work=NULL;
}




















