###############################################################

#NON LINEAR AND NON STATIONNARY UNIVARIARE MODEL

###############################################################

rm(list=ls())

library(RcppSQMC)

#################################################################
#parameters
###############################################################

a=20
b1=0.5
b2=25
b3=8
b4=1.2
sigma=sqrt(10)
x0=0.1
parPrior=c(0,sqrt(2))	#parameters of the gaussian prior

theta=c(b1,b2,b3,b4,sigma,a, parPrior)

#################################################################
#Data
###############################################################

T=100

yy=NL_sim(T,theta,x0)

y=yy$Y
X=yy$X


#Data used in the paper "Sequential quasi-Monte Carlo"

#y=as.matrix(as.numeric(read.table("Data/ToyModel/y_Toy.txt")$V1))
#X=as.matrix(as.numeric(read.table("Data/ToyModel/X_Toy.txt")$V1))
#theta=as.matrix(as.numeric(read.table("Data/ToyModel/theta_Toy.txt")$x))

T=length(y)


#################################################################
## FILTERING
###############################################################

N=2^7

timer =proc.time()[3]
iter=SQMC_Univ1(y,theta, N, seed=-1, computeExp=1, qmc=1, src=1, ns=2)
proc.time()[3]-timer

print(mean(iter$L[T,]))				#mean estimate of the log-likelihood
print(iter$EX[T,])				#mean estimate of E[x_T|y_{1:T}]
print(sqrt(iter$EX2[T,]-iter$EX[T,]^2))		#std of the ns estimates of E[x_T|y_{1:T}]
print(X[T,])

###############################################################
## FORWARD FILTERING-BACKWARD SMOOTHING
###############################################################


N=2^7	##number of particles for the forward pass
Nb=2^7	##number of particles for the backward pass

timer =proc.time()[3]
iter=SQMCBack_Univ1(y,  theta, N, Nb, seed=-1, qmc=1,  qmcb=1, Marg=0, ns=2)
proc.time()[3]-timer


t=floor(T/2)
print(mean(iter$L[T,]))				#mean estimate of the log-likelihood
print(iter$EX[t,])				#mean estimate of E[x_t|y_{1:T}]
print(sqrt(iter$EX2[t,]-iter$EX[t,]^2))		#std of the ns estimates of E[x_t|y_{1:T}]
print(X[t,])

























