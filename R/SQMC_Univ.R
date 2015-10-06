############################################################################################################################################
#		NON LINEAR AND NON STATIONNARY UNIVARIARE MODEL
############################################################################################################################################


############################################################################################################################################
#SQMC_Univ: FILTERING AND LIKLIHOOD ESTIMATION

#INPUTS:

#y=(y_0,...,Y_{T-1}): vector of data
#theta	 	    : vector of parameters
#N		    : number of particles
#seed		    : used for SMC (if seed<0, seed is taken at random)
#computeExp	    : computeExp=1 if the algorithm compute E[X_t|y_{0:t}]
#qmc		    : qmc=0 for SMC, qmc=1 for Sobol and qmc=2 for Niederreiter in base 2.
#src		    : src=1 to scramble the sequence,  src=0 otherwise
#ns		    : number of SMC/SQMC runs

#OUTPUTS

#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : vector of size T giving estimate of E[X_t|Y_{0:t}] at time t=0,...,T-1
#EX2		    : vector of size T giving estimate of E[X^2_t|Y_{0:t}] at time t=0,...,T-1

############################################################################################################################################


SQMC_Univ1=function(y, theta, N, seed=-1, computeExp=0, qmc=0, src=1, ns=1)
{
	if(qmc>2 || qmc<0)
	{
		print("Error: Not allowed: qmc should be 0,1 or 2")

	}else if( (qmc<3) && (src>1 || src<0))
	{
		print("Error: Not allowed: src should be 0 or 1")
	}else
	{
		T=length(y)
		dx=1
		dy=1

		filter=SQMC_Univ(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(T), as.double(theta), as.integer(seed), as.integer(ns), as.integer(N),as.integer(qmc), as.integer(src),
		as.integer(computeExp), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))

		return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))
	}
}

############################################################################################################################################
	#SQMC_Back_SV: QMC FORWARD-BACKWARD SMOOTHING ALGORITHM

#INPUTS:

#y=(y_0,...,Y_{T-1}): vector of data
#theta	 	    : vector of parameters
#N		    : number of particles used in the forward step
#Nb		    : number of particles used in the backward step
#seed		    : used for SMC (if seed<0, seed is taken at random)
#qmc		    : qmc=0 for SMC, qmc=1 for scrambled Sobol and qmc=2 for scrambled Niederreiter in base 2 (forward step)
#qmcb		    : qmc=0 for SMC, qmc=1 for  scrambled Sobol, qmc=2 for  scrambled Niederreiter in base 2 (backward step)
#Marg		    : Marg=1 if target marginal smoothing distributions, 0 otherwise
#ns		    : number of SMC/SQMC runs


#OUTPUTS
#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : matrix of size (T times dx) giving estimate of E[X_t|Y_{0:T-1}] at time t=0,...,T-1
#EX2		    : matrix of size (T times dx) giving estimate of E[X^2_t|Y_{0:T-1}] at time t=0,...,T-1

############################################################################################################################################


SQMCBack_Univ1=function(y, theta, N, Nb, seed=-1, qmc=0,  qmcb=0, Marg=0, ns=1)
{

	if(qmc>2 || qmc<0)
	{
		print("Error: Not allowed: qmc should be 0,1 or 2")
	}else
	{

		T=length(y)
		dx=1
		dy=1

		filter=SQMCBack_Univ(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(T),
		as.double(theta), as.integer(seed), as.integer(ns), as.integer(N),as.integer(Nb), as.integer(qmc), as.integer(qmcb),
		as.integer(Marg), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))

		return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))


	}

}











