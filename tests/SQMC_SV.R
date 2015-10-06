###############################################################

#MULTIVARIATE STOCHASTIC VOLATILITY MODEL

###############################################################

rm(list=ls())

library(RcppSQMC)


###############################################################
## parameters
###############################################################

dx=2
x0=mu=rep(-9,dx)
Phi=diag(0.90,dx)
Psi=diag(0.1,dx)
C1=matrix(0.6,dx,dx)+diag(0.4,dx,dx)
C2=matrix(0.8,dx,dx)+diag(0.2,dx,dx)
C12=matrix(-0.1,dx,dx)+diag(-0.2,dx,dx)


C=rbind(cbind(C1,C12), cbind(t(C12),C2))
theta=c(mu,Phi,Psi,C)

###############################################################
#data
###############################################################

#Generate data
T=100

Y=SVM_sim(T,mu, Phi, Psi, C, x0)
y=Y$Y
X=Y$X

#Data used in the paper "Sequential quasi-Monte Carlo"

##1D-model
#y=apply(read.table("Data/SVModel/y_SV1.txt"),2,as.numeric)
#X=apply(read.table("Data/SVModel/X_SV1.txt"),2,as.numeric)
#theta=as.matrix(as.numeric(read.table("Data/SVModel/theta_SV1.txt")$x))

##2D-model
#y=apply(read.table("Data/SVModel/y_SV2.txt"),2,as.numeric)
#X=apply(read.table("Data/SVModel/X_SV2.txt"),2,as.numeric)
#theta=as.matrix(as.numeric(read.table("Data/SVModel/theta_SV2.txt")$x))

##4D-model
#y=apply(read.table("Data/SVModel/y_SV4.txt"),2,as.numeric)
#X=apply(read.table("Data/SVModel/X_SV4.txt"),2,as.numeric)
#theta=as.matrix(as.numeric(read.table("Data/SVModel/theta_SV4.txt")$x))

T=nrow(y)
dx=ncol(y)


###############################################################
## FILTERING
###############################################################

N=1000

timer =proc.time()[3]
iter=SQMC_SV1(y,dx,theta, N, seed=-1, computeExp=1, qmc=1, src=1, ns=2)
proc.time()[3]-timer


print(mean(iter$L[T,]))				#mean estimate of the log-likelihood
print(iter$EX[T,])				#mean estimate of E[x_T|y_{1:T}]
print(sqrt(iter$EX2[T,]-iter$EX[T,]^2))		#std of the ns estimates of E[x_T|y_{1:T}]
print(X[T,])

N=2^7
Nb=2^7

##forward-backward algorithm
timer =proc.time()[3]
iter=SQMCBack_SV1(y,dx,theta, N, Nb, seed=-1,  qmc=1,  qmcb=1, Marg=1, ns=2)
proc.time()[3]-timer

t=100
print(mean(iter$L[T,]))				#mean estimate of the log-likelihood
print(iter$EX[t,])				#mean estimate of E[x_t|y_{1:T}]
print(sqrt(iter$EX2[t,]-iter$EX[t,]^2))		#std of the ns estimates of E[x_t|y_{1:T}]
print(X[t,])


##two filter algorithm (Need C12=matrix(0,dx,dx))
timer =proc.time()[3]
iter=SQMC2F_SV1(y, dx, theta, N, Nb, floor(T/2), seed=-1, qmc=1,  ns=2)
proc.time()[3]-timer

###############################################################
## PMMH-SQMC
###############################################################
#Real data example

yySP=read.csv("Data/SVModel/SP500.csv", sep=",", header=TRUE)
yyN=read.csv("Data/SVModel/Nasdac.csv", sep=",", header=TRUE)

ySP=as.numeric(as.character(yySP$Close))
yN=as.numeric(as.character(yyN$Close))

ydata=matrix(0,length(yN)-1,2)


for(i in 2:length(yN))
{
    ydata[i-1,1]=log(ySP[i])-log(ySP[i-1])
    ydata[i-1,2]=log(yN[i])-log(yN[i-1])
}

##mean corrected data

y=matrix(0,length(yN)-1,2)

y[,1]=ydata[,1]-mean(ydata[,1])
y[,2]=ydata[,2]-mean(ydata[,2])


dx=ncol(y)
T=nrow(y)

###############################################################
#Prior distribution (assume no leverage effect)

parPrior=c(10*exp(-10), 10*exp(-3)) #parameters of the gamma priors

priorSV=function(t0,dx, dy,par)	    #prior distribution
{
    nn=(length(t0))

    t=c(t0[1:(nn-1)],rep(0, 2*dx),t0[nn])

    if(sum(t[(dx+1):(2*dx)]^2<1)==dx & sum(t[(2*dx+1):(3*dx)]>0)==dx & sum(t0[(nn-1):nn]^2<1)==dx)
    {
        d=0

        for(i in 1:dx)
        {
            d=d-(1+par[1])*log(t[2*dx+i])-par[2]/t[2*dx+i]
        }

        return(d)
    }else
    {
        return(-Inf)
    }
}


#put the vector of parameters t0 in the "right" format to be used in "SQMC_SV"

parFilterSV=function(t0,dx,dy)
{
    t=c(t0[1:(length(t0)-1)],rep(0, 2*dx),t0[length(t0)])
    mu=t[1:dx]

    Phi=diag(t[(dx+1):(2*dx)],dx)

    Psi=diag(t[(2*dx+1):(3*dx)],dx)

    C=diag(1,2*dx)

    d=1
    for(i in 1:(2*dx-1))
    {
        for(j in (i+1):(2*dx))
        {
            C[i,j]=C[j,i]=t[3*dx+d]
            d=d+1
        }
    }
    return(c(mu,Phi,Psi,C))
}

###############################################################
#initial value of the MArkov chain

theta0=as.numeric(read.table("Data/SVModel/theta0.txt")$x)


###############################################################
#parameter of PMMH-SQMC

m=200				#length of the Markov chain
N=50				#number of particles used in SQMC


#covariance matrix of the gaussian proposal

c=read.table("Data/SVModel/cov.txt")
c=matrix(as.numeric(as.matrix(c)),8,8)

###############################################################
#PMMH-SQMC

iter=Marginal_PMCMC(y, dx, m, c,  theta0, priorSV, parPrior, N, qmc=1, src=1, parFilterSV, SQMC_SV1)


print(iter$P)			#acceptance rate
apply(iter$Y,1,mean)		#posterior mean


#effectiveSize(mcmc(t(iter$Y)))	#compute effective sample size
