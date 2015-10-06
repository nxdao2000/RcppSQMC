## HILBERT RESAMPLING
library(RcppSQMC)
library(gtools)
N=10					#number of particles
dx=8					#dimension of the state space

x=runif(dx*N)				#set of partilces (x^1,...,x^N)

W=c(rdirichlet(1, rep(1,N)))		#normalized weights

parPhi=rep(c(0,1),dx)			#parPhi=c(b_1,c_1,...,b_d,c_d):	Used to rescale the particles:
					#x_i--> (x_i-b_i)/(c_i-b_i), i=1,...,d

quasi=sort(runif(N))			#A set of N points in (0,1), sorted by increasing order.


Index=HilbertResampler(as.double(quasi), as.double(parPhi), as.double(x), as.integer(dx), as.integer(N), as.double(W), double(N))

Index					#Output: N integers in {1,...,N}


##GENERATE Owen (1995) nested scrambled (t,s)-sequence

##NOTE: The point are sorted according to the first coordinate



N=2^7					#Number of points in each point sets
dx=2					#dimension of the point sets
qmc=1					#qmc=1 for Sobol' sequence and qmc=2 for Niederreiter
					#sequence in base b=2

ns=1					#number of point sets to generate
#if(FALSE){
points=RQMC(as.integer(qmc), as.integer(N), as.integer(dx), as.integer(ns), double(dx*N*ns))

res=array(points,c(dx,N,ns))		#res[,,k] is a (dx times N) matrix
					#that contains the points of the
					#k-th sample, k=1,...,ns

k=1					#select sample: k\in{1,...,ns}
plot(res[1,,k],res[2,,k]) 		#plot first two dimensions

#}
##GENERATE (t-s)-sequences
if(FALSE){

N=128					#Number of points to generate
dx=2					#dimension of the sequence
qmc=1					#qmc=1 for Sobol' sequence and qmc=2 for Niederreiter
					#sequence in base b=2

points=QMC(as.integer(qmc), as.integer(N), as.integer(dx), double(dx*N))

res=matrix(points,dx,N)			#(dx times N) matrix containing the point set

plot(res[1,],res[2,]) 			#plot first two dimensions

}



















