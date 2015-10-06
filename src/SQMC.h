

#include<Rcpp.h>
#include<stdlib.h>
#include<stdio.h>


#include "Resampler.h"
#include "functions.h"
#include "generate_RQMC.hpp"
#include "SMC.h"



typedef void (*SimTransitionPF_QMC)(double*, double*, int, int, int, int, double*, double*, double*);

typedef void (*SimInitPF_QMC)(double*, int, int, int, double*, double*);

typedef void (*ResamplingPF_QMC)(double*, double*, double*, int, int, double*, double*);

typedef int* (*ResamplingBack_QMC)(double*, int , double*, double*, int, int, int, double*, double*);

//typedef void *DigitalNetGenerator;

//typedef void *Scrambled;

/*******************************************************************/

//SQMC Algorithm


void SQMC(double*, int, int, int, double*, int, int, ParamTransitionPF, ResamplingPF_QMC,  SimInitPF_QMC, SimTransitionPF_QMC,
PotentialPF, int*, double*, double*, double*);

void SQMC_Forward(double*, int, int, int, double*, int, int, ParamTransitionPF, ResamplingBack_QMC,  SimInitPF_QMC, SimTransitionPF_QMC,
PotentialPF, int*, double*, double*, double*, double*, double*, double*);

void SQMC_BIF(double*, int, int, int, double*, int, int, ParamTransitionPF, ResamplingBack_QMC,  SimInitPF_QMC, SimTransitionPF_QMC,
PotentialPF, int*, double*, double*, double*, double*, double*, double*);

/*******************************************************************/

extern "C"
void SQMC_SV(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int qmc, int src, int computeExp, double *parPsi, double *lik, double *expx, double *expx2);
extern "C"
void SQMCBack_SV(double *y, int dy, int dx, int T, double *theta, int seed, int ns, int N, int Nb, int qmc, int qmcB, int Marg, double *parPsi, double *lik, double *expx, double *expx2);
























