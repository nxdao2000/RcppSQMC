#include <Rcpp.h>

#include <stdlib.h>
#include <stdlib.h>


#include "HilbertCode.hpp"
#include "functions.h"

extern "C"
void HilbertResampler(double *quasi, double *parPhi, double *x, int dx, int N, double *W, double* xh);
void quasi_Resample(double*, int*, double*, int *, double*, double*);
int* ForHilbertResampler(double *quasi, int dimQuasi, double *parPhi, double *x, int dx, int N, int Nb, double *W, double* xh);




void quasi_ResampleFor(double *quasi, int dimQuasi, int dx,  int N, int Nb, double *W, int *J1, double *x, double *xh);


int* BackHilbertResampler(double *quasi, int dimQuasi, double *parPhi, double *x, int dx, int N, int Nb, double *W, double* xh);





void quasi_ResampleBack(double *quasi, int dimQuasi, int dx,  int N, int Nb, double *W, int *J1, int *J2, double *x, double *xh);
