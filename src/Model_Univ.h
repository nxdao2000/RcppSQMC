


#include "functions.h"


void simPrior_Univ(gsl_rng*, int, int, int, double*, double*);
void simPriorQMC_Univ(double*, int, int, int, double*, double*);

double *paramTrans_Univ(double*, int, int);

void simTransition_Univ(gsl_rng*, double*, int, int, int, int, double*, double*, double*);
void simTransitionQMC_Univ(double*, double*, int, int, int, int, double*, double*, double*);

void potential_Univ(double*, double*, double*, int,  int, int, int, double*, double*);

void transition_Univ(double*, double*, int, int, int, double*, double*);
