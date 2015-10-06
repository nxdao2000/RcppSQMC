###############################################################################################################################################

# RANDOM WALK PARTICLE MARGINAL METROPOLIS HASTINGS ALGORITHM

#parameters:

#y=(y_0,...,Y_{T-1}): dy*T vector of data
#dx		    : dimensions of of the hidden process
#m		    : length of the Markov chain
#theta0	 	    : starting point of the Markov chain
#prior		    : density function of the prior distribution
#parPrior	    : parameters of "prior"
#N		    : number of particles for the filter
#qmc		    : qmc=1 for Sobol, qmc=2 for Niederreiter in base 2 and qmc=3 for Niederreiter-Xing. For other values of qmc SMC is used
#src		    : scrambling method: src=0 for linear scrambling and src=1 for full (i.e. Owen's) scrambling
#parFilter	    : functions that puts "theta0" in the right format to be used in the filter
#filter		    : filter algorithm

################################################################################################################################################

library("mnormt")

Marginal_PMCMC=function(y, dx, m, c, theta0, prior, parPrior, N, qmc, src, parFilter, filter)
{
	n=nrow(y)
	dy=ncol(y)
	T=nrow(y)
	np=length(theta0)
	theta=matrix(0,np,m)
	

	theta[,1]=theta0
	d=0

	work=parFilter(theta[,1],dx,dy)
	p0=prior(theta[,1], dx, dy, parPrior)+filter(y, dx, work, N, seed=-1, computeExp=0, qmc, src,ns=1)$L[T]

	p1=0

	for(i in 2:m)
	{
		theta1=theta[,(i-1)]+rmnorm(1,rep(0,np),c)
		
		p1=prior(theta1, dx, dy, parPrior)

		if(p1==-Inf)
		{
			theta[,i]=theta[,(i-1)]
		}else
		{
			work=parFilter(theta1,dx,dy)
			p1=p1+filter(y, dx, work, N, seed=-1, computeExp=0, qmc, src, ns=1)$L[T]

			if(log(runif(1))>p1-p0 | p1=="NaN")
			{

				theta[,i]=theta[,(i-1)]

			}else
			{
				theta[,i]=theta1
				p0=p1
				d=d+1
			}	
		}	
	}
	return(list(Y=theta, P=d/m))
}


























