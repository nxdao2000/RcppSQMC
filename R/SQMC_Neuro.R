############################################################################################################################################
#		NEURO DECODING MODEL
############################################################################################################################################

#SQMC_Neuro: FILTERING AND LIKLIHOOD ESTIMATION

#INPUTS:

#y=(y_0,...,Y_{T-1}): dx*T vector of data
#dx		    : dimensions of of the hidden process
#theta	 	    : vector of parameters, theta=c(mu,Phi,Psi,C)
#N		    : number of particles
#seed		    : used for SMC (if seed<0, seed is taken at random)
#computeExp	    : computeExp=1 if the algorithm compute E[X_t|y_{0:t}]
#qmc		    : qmc=0 for SMC, qmc=1 for Sobol and qmc=2 for Niederreiter in base 2.
#src		    : src=1 to scramble the sequence,  src=0 otherwise
#ns		    : number of SMC/SQMC runs

#OUTPUTS

#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : matrix of size (T times dx) giving estimate of E[X_t|Y_{0:t}] at time t=0,...,T-1
#EX2		    : matrix of size (T times dx) giving estimate of E[X^2_t|Y_{0:t}] at time t=0,...,T-1

############################################################################################################################################


SQMC_Neuro1=function(y, dx, theta0, N, seed=-1, computeExp=0, qmc=0, src=1, ns=1)
{
	if(qmc>2 || qmc<0)
	{
		print("Error: Not allowed: qmc should be 0,1, or 2")

	}else if( (0<qmc && qmc<3) && (src>1 || src<0))
	{
		print("Error: Not allowed: src should be 0 or 1")
	}else
	{
		T=nrow(y)
		dy=ncol(y)

		nay=dx+dy+dx^2+dx*dy

		theta=c(rep(0,dx),theta0)

		Sx=matrix(theta[(nay+1):(nay+dx^2)],dx,dx)
		Sy=matrix(0,dy,dy)
		mm=length(theta)
		Phi=matrix(theta[(dx+dy+1):(dx+dy+dx^2)],dx,dx,byrow=T)
		beta=matrix(theta[(dx+dy+dx^2+1):(dx+dy+dx^2+dx*dy)],dy,dx, byrow=T)

		vec=diag(1,dx) #variance prior distributiob

		theta2=c(theta[1:dx],log(theta[mm])+theta[(dx+1):(dx+dy)],t(Phi), t(beta), Sx,Sy,c(vec))

		## parameters used to map the particles in (0,1)^d

		PP1=matrix(vec,dx,dx)
		xxx1=theta[1:dx]

		parPsi=rep(0,2*dx)
		for(k in 1:dx)
		{
			parPsi[2*k-1]=xxx1[k]-2*sqrt(PP1[k,k])
			parPsi[2*k]=xxx1[k]+2*sqrt(PP1[k,k])
		}

		filter=SQMC_Neuro(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(T),
		as.double(theta2), as.integer(seed), as.integer(ns), as.integer(N),as.integer(qmc), as.integer(src),
		as.integer(computeExp), as.double(parPsi), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))

		return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))
	}

}












