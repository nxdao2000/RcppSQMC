



#include "functions.h"


void simLinearGaussian(gsl_rng*, int, int, int, int, double*, double*);
void simLinearGaussianQMC(double*, int, int, int, int, double*, double*);

void simPrior_LG(gsl_rng*, int, int, int, double*, double*);
void simPriorQMC_LG(double*, int, int, int, double*, double*);

double* paramTrans_LG(double*, int, int);
double* paramTransQMC_LG(double*, int, int);

void simTransition_LG(gsl_rng*, double*, int, int, int, int, double*, double*, double*);
void simTransitionQMC_LG(double*, double*, int, int, int, int, double*, double*, double*);



