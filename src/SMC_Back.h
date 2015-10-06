

#include "SMC.h"
#include "functions.h"



void SMC_ForBack(double*, int, int, int, double*, int, int, gsl_rng*, ParamTransitionPF, ResamplingPF, SimInitPF, TransitionPF, SimTransitionPF, PotentialPF, int*, double*, double*);

void SMC_Back(int, int, int, double*, int, int, gsl_rng*, ResamplingPF, TransitionPF, double*, double*, double*, double*, double*, double*);

void SMC_ForBackMarg(double*, int, int, int, double*, int, int, gsl_rng*, ParamTransitionPF, ResamplingPF, SimInitPF, TransitionPF, SimTransitionPF, PotentialPF, int*, double*, double*);

void SMC_BackMarg(int, int,int, double*, int, int, gsl_rng*, ResamplingPF, TransitionPF, double*, double*, double*, double*, double*, double*);


void SMC_2F(double*, int, int, int, int, double*, double*, int, int,  gsl_rng*, ParamTransitionPF, ResamplingPF, SimInitPF,
 TransitionPF, SimTransitionPF, PotentialPF, ParamTransitionPF, SimInitPF, TransitionPF, SimTransitionPF, PotentialPF,  int* ,
 double*, double*);
