############################################################################################################################################
#		MULTIVARIATE STOCHASTIC VOLATILITY MODEL
############################################################################################################################################

#SQMC_Univ: FILTERING AND LIKLIHOOD ESTIMATION

#INPUTS:

#y=(y_0,...,Y_{T-1}): dx*T vector of data
#dx		    : dimensions of of the hidden process
#theta	 	    : vector of parameters, theta=c(mu,Phi,Psi,C)
#N		    : number of particles
#seed		    : used for SMC (if seed<0, seed is taken at random)
#computeExp	    : computeExp=1 if the algorithm compute E[X_t|y_{0:t}]
#qmc		    : qmc=0 for SMC, qmc=1 for Sobol and qmc=2 for Niederreiter in base 2
#src		    : src=1 to scramble the sequence,  src=0 otherwise
#ns		    : number of SMC/SQMC runs

#OUTPUTS

#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : matrix of size (T times dx) giving estimate of E[X_t|Y_{0:t}] at time t=0,...,T-1
#EX2		    : matrix of size (T times dx) giving estimate of E[X^2_t|Y_{0:t}] at time t=0,...,T-1


#IMPORTANT NOTE: the prior distribution for X_0 assumed that the latent process is stationnary

############################################################################################################################################


SQMC_SV1=function(y, dx, theta, N, seed=-1, computeExp=0, qmc=0, src=1, ns=1)
{
	if(qmc>2 || qmc<0)
	{
		print("Error: Not allowed: qmc should be 0,1 or 2")

	}else if( (0<qmc && qmc<3) && (src>1 || src<0))
	{
		print("Error: Not allowed: src should be 0 or 1")
	}else
	{
		T=nrow(y)
		dy=ncol(y)

		##precompute some quantities that will be used to evaluate the potential functions

		mu=theta[1:dx]

		Phi=matrix(theta[(dx+1):(dx+dx^2)],dx,dx,byrow=F)

		Psi=matrix(theta[(dx+dx^2+1):(dx+2*dx^2)],dx,dx,byrow=F)

		C=matrix(theta[(dx+2*dx^2+1):(dx+6*dx^2)],2*dx,2*dx,byrow=F)

		C1=C[1:dx,1:dx]
		C2=C[(dx+1):(2*dx),(dx+1):(2*dx)]
		C12=C[1:dx,(dx+1):(2*dx)]

		Sx=sqrt(Psi)%*%C2%*%sqrt(Psi)

		Sy=C1-C12%*%solve(C2)%*%t(C12)
		A=C12%*%solve(C2)%*%solve(sqrt(Psi))
		APhi=A%*%Phi

		theta2=c(mu,t(APhi),t(A), solve(Sy), Sx, t(Phi),log(det(Sy)))

		##parameter of the prior distribution (stationnary distribution of X_t)

		vec=solve(diag(1,dx^2)-Phi%x%Phi)%*%matrix(c(Sx),dx^2)

		theta2=c(theta2,mu,vec)

		## parameters used to map the particles in (0,1)^d

		PP1=matrix(vec,dx,dx)
		xxx1=mu

		parPsi=rep(0,2*dx)
		for(k in 1:dx)
		{
			parPsi[2*k-1]=xxx1[k]-2*sqrt(PP1[k,k])
			parPsi[2*k]=xxx1[k]+2*sqrt(PP1[k,k])
		}


		filter=SQMC_SV(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(T),
		as.double(theta2), as.integer(seed), as.integer(ns), as.integer(N),as.integer(qmc), as.integer(src),
		as.integer(computeExp), as.double(parPsi), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))

		return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))

	}
}

############################################################################################################################################
	#SQMC_Back_SV: QMC FORWARD-BACKWARD SMOOTHING ALGORITHM

#INPUTS:

#y=(y_0,...,Y_{T-1}): dx*T vector of data
#dx		    : dimensions of of the hidden process
#theta	 	    : vector of parameters, theta=c(mu,Phi,Psi,C)
#N		    : number of particles used in the forward step
#Nb		    : number of particles used in the backward step
#seed		    : used for SMC (if seed<0, seed is taken at random)
#qmc		    : qmc=0 for SMC, qmc=1 for scrambled Sobol and qmc=2 for scrambled Niederreiter in base 2 (forward step)
#qmcb		    : qmc=0 for Monte Carlo, qmc=1 for scrambled Sobol and  qmc=2 for scrambled Niederreiter in base 2  (backward step)
#Marg		    : Marg=1 if target matrginal smoothing distributions, 0 otherwise
#ns		    : number of SMC/SQMC runs

#OUTPUTS
#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : matrix of size (T times dx) giving estimate of E[X_t|Y_{0:T-1}] at time t=0,...,T-1
#EX2		    : matrix of size (T times dx) giving estimate of E[X^2_t|Y_{0:T-1}] at time t=0,...,T-1


#IMPORTANT NOTE: the prior distribution for X_0 assumed that the latent process is stationnary

############################################################################################################################################

SQMCBack_SV1=function(y, dx, theta, N, Nb, seed=-1,  qmc=0, qmcb=0, Marg=0, ns=1)
{
    if(qmc>2 || qmc<0)
    {
        print("Error: Not allowed: qmc should be 0,1 or 2")

    }else
    {
        T=nrow(y)
        dy=ncol(y)

        ##precompute some quantities that will be used to evaluate the potential functions

        mu=theta[1:dx]

        Phi=matrix(theta[(dx+1):(dx+dx^2)],dx,dx,byrow=F)

        Psi=matrix(theta[(dx+dx^2+1):(dx+2*dx^2)],dx,dx,byrow=F)

        C=matrix(theta[(dx+2*dx^2+1):(dx+6*dx^2)],2*dx,2*dx,byrow=F)

        C1=C[1:dx,1:dx]
        C2=C[(dx+1):(2*dx),(dx+1):(2*dx)]
        C12=C[1:dx,(dx+1):(2*dx)]

        Sx=sqrt(Psi)%*%C2%*%sqrt(Psi)

        Sy=C1-C12%*%solve(C2)%*%t(C12)
        A=C12%*%solve(C2)%*%solve(sqrt(Psi))
        APhi=A%*%Phi

        theta2=c(mu,t(APhi),t(A), solve(Sy), Sx, t(Phi),log(det(Sy)))

        ##parameter of the prior distribution (stationnary distribution of X_t)

        vec=solve(diag(1,dx^2)-Phi%x%Phi)%*%matrix(c(Sx),dx^2)

        theta2=c(theta2,mu,vec)

        ## parameters used to map the particles in (0,1)^d

        PP1=matrix(vec,dx,dx)
        xxx1=mu

        parPsi=rep(0,2*dx)
        for(k in 1:dx)
        {
            parPsi[2*k-1]=xxx1[k]-2*sqrt(PP1[k,k])
            parPsi[2*k]=xxx1[k]+2*sqrt(PP1[k,k])
        }

        filter=SQMCBack_SV(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(T),as.double(theta2), as.integer(seed), as.integer(ns), as.integer(N), as.integer(Nb), as.integer(qmc),
                           as.integer(qmcb),  as.integer(Marg), as.double(parPsi), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))

       return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))

    }
}

############################################################################################################################################
#SQMC2F_SV: QMC TWO FILTER SMOOTHING ALGORITHM

#INPUTS:

#y=(y_0,...,Y_{T-1}): dx*T vector of data
#dx		    : dimensions of of the hidden process
#theta	 	    : vector of parameters, theta=c(mu,Phi,Psi,C)
#N		    : number of particles used in the forward step (for the two filters)
#Nb		    : number of particles used in the backward step
#tstar		    : starting date of the information filter
#seed		    : used for SMC (if seed<0, seed is taken at random)
#qmc		    : qmc=0 for SMC, qmc=1 for scrambled Sobol and for scrambled  qmc=2 for Niederreiter in base 2 (forward and backward step)
#ns		    : number of SMC/SQMC runs

#OUTPUTS
#L		    : matrix of size (T times ns) giving ns estimations of the log-likelihood at time t=0,...,T-1
#EX		    : matrix of size (T times dx) giving estimate of E[X_t|Y_{0:T-1}] at time t=0,...,T-1
#EX2		    : matrix of size (T times dx) giving estimate of E[X^2_t|Y_{0:T-1}] at time t=0,...,T-1


#IMPORTANT NOTE: the prior distribution for X_0 assumed that the latent process is stationnary

############################################################################################################################################



SQMC2F_SV1=function(y, dx, theta, N, Nb, tstar, seed=-1, qmc=0, ns=1)
{
    if(qmc>2 || qmc<0)
    {
        print("Error: Not allowed: qmc should be 0,1 or 2")

    }else
    {
        T=nrow(y)
        dy=ncol(y)

        ##precompute some quantities that will be used to evaluate the potential functions

        mu=theta[1:dx]

        Phi=matrix(theta[(dx+1):(dx+dx^2)],dx,dx,byrow=F)

        Psi=matrix(theta[(dx+dx^2+1):(dx+2*dx^2)],dx,dx,byrow=F)

        C=matrix(theta[(dx+2*dx^2+1):(dx+6*dx^2)],2*dx,2*dx,byrow=F)

        C1=C[1:dx,1:dx]
        C2=C[(dx+1):(2*dx),(dx+1):(2*dx)]
        C12=C[1:dx,(dx+1):(2*dx)]

        Sx=sqrt(Psi)%*%C2%*%sqrt(Psi)

        Sy=C1-C12%*%solve(C2)%*%t(C12)
        A=C12%*%solve(C2)%*%solve(sqrt(Psi))
        APhi=A%*%Phi

        theta2=c(mu,t(APhi),t(A), solve(Sy), Sx, t(Phi),log(det(Sy)))

        ##parameters of the prior distribution (stationnary distribution of X_t)

        vec=solve(diag(1,dx^2)-Phi%x%Phi)%*%matrix(c(Sx),dx^2)

        theta2=c(theta2,mu,vec)

        ##parameters for the information filter)

        Sigma0=matrix(vec,dx,dx)

        Phi2F=Sigma0%*%Phi%*%solve(Sigma0)

        Sx2F=Sigma0-Sigma0%*%Phi%*%solve(Sigma0)%*%Phi%*%Sigma0

        theta2F=c(mu,t(APhi),t(A), solve(Sy), Sx2F, t(Phi2F),log(det(Sy)),mu,vec)

        ## parameters used to map the particles in (0,1)^d

        PP1=matrix(vec,dx,dx)
        xxx1=mu

        parPsi=rep(0,2*dx)
        for(k in 1:dx)
        {
            parPsi[2*k-1]=xxx1[k]-2*sqrt(PP1[k,k])
            parPsi[2*k]=xxx1[k]+2*sqrt(PP1[k,k])
        }

        filter=SQMC2F_SV(as.double(c(t(y))), as.integer(dy), as.integer(dx), as.integer(tstar), as.integer(T), as.double(theta2), as.double(theta2F), as.integer(seed), as.integer(ns), as.integer(N),as.integer(Nb), as.integer(qmc), as.double(parPsi), lik=double(ns*T), expx=double(dx*T), expx2=double(dx*T))


        return(list(L=matrix(filter$lik,T,ns), EX=matrix(filter$expx,T,dx,byrow=T), EX2=matrix(filter$expx2,T,dx,byrow=T)))


    }
}










