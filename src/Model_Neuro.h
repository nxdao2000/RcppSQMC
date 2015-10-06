#include "functions.h"
#include "Model_LG.h"



double logfactorial(int);

void potential_Neuro(double*, double*, double*, int,  int, int, int, double*, double*);


double* paramTrans_Neuro(double*, int, int);
double* paramTransQMC_Neuro(double*, int, int);

void simTransition_Neuro(gsl_rng*, double*, int, int, int, int, double*, double*, double*);
void simTransitionQMC_Neuro(double*, double*, int, int, int, int, double*, double*, double*);

