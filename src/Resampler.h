#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_heapsort.h> 

#include "functions.h"


void ResamplePF(gsl_rng*, double*, int, int,  int , double*, double*);

int* ResamplePF_index(gsl_rng*, double*, int, int,  int , double*, double*);

void ResampleBack(double, double*, int, int, double*, int, double*);


void ResampleMargBackQMC(double*, double*, int, int, int, double*,  double*, double*);

void ResampleMargBackMC(gsl_rng*, double*, int, int, int, double*,  double*, double*);
