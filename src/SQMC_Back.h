


#include<stdlib.h>  		
#include<stdio.h>	

#include "SQMC.h"


//QMC forward pass

void SQMC_Back(int, int, int, double*, int, int, int, ResamplingBack_QMC, TransitionPF, double*, double*, double*, double*, double*, double*, double*);


//SQMC Forward-Backward algorithm for additive functions

void SQMC_BackMarg(int, int, int, double*, int, int, int, ResamplingBack_QMC, TransitionPF, double*, double*, double*, double*, double*, double*, double*);



void SQMC_2F(double*, int, int, int, int, double*, double*, int, int, int, ParamTransitionPF, ResamplingBack_QMC, ResamplingBack_QMC, SimInitPF_QMC,
 TransitionPF, SimTransitionPF_QMC, PotentialPF, ParamTransitionPF,  SimInitPF_QMC, TransitionPF, SimTransitionPF_QMC, 
PotentialPF, int*, double*, double*, double*);























