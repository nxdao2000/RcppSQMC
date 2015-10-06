###############################################################

#	NEURO DECODING MODEL

###############################################################

rm(list=ls())

library(RcppSQMC)

###############################################################
## parameters
###############################################################

dx=4
dy=10
delta=0.03
alpha=2.5+rnorm(dy,0,1)
beta=matrix(runif(dy*dx),dy,dx)
Phi=matrix(0,dx,dx)
Phi[1:(dx/2),1:(dx/2)]=diag(1,dx/2,dx/2)
Phi[(dx/2+1):dx,(dx/2+1):dx]=diag(1,dx/2,dx/2)
Phi[(dx/2+1):dx,1:(dx/2)]=delta*matrix(1,dx/2,dx/2)
sigma2=0.019
x0=rep(0,dx)


Sx=diag(c(rep(sigma2,dx/2),rep(0,dx/2)))
theta=c(alpha,t(Phi),t(beta),t(Sx),delta)


###############################################################
## data
###############################################################

#Generate data

T=23

yy=Neuro_sim(T,alpha,beta,delta,Phi,Sx,x0)
y=yy$Y
X=yy$X


#Data used in the paper "Sequential quasi-Monte Carlo"

#y=apply(read.table("Data/NeuroModel/y_Neuro.txt"),2,as.numeric)
#X=apply(read.table("Data/NeuroModel/X_Neuro.txt"),2,as.numeric)
#theta=as.matrix(as.numeric(read.table("Data/NeuroModel/theta_Neuro.txt")$x))

T=nrow(y)
dx=ncol(X)

###############################################################
## FILTERING
###############################################################

N=300

timer =proc.time()[3]
iter=SQMC_Neuro1(y,dx,theta, N, seed=-1, computeExp=1, qmc=1, src=1, ns=2)
proc.time()[3]-timer


print(mean(iter$L[T,]))				#mean estimate of the log-likelihood
print(iter$EX[T,])				#mean estimate of E[x_T|y_{1:T}]
print(sqrt(iter$EX2[T,]-iter$EX[T,]^2))		#std of the ns estimates of E[x_T|y_{1:T}]
print(X[T,])








