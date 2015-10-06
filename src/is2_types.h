/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Digital nets in base 2                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Copyright (C) 2002 by Thomas Kollig (kollig@informatik.uni-kl.de)       //
//                       Alexander Keller (keller@informatik.uni-kl.de)    //
//                                                                         //
// All rights reserved. You may not distribute this software, in whole or  //
// in part, especially not as part of any commercial product, without the  //
// express consent of the authors.                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef IS2_TYPES_INCLUDED
#define IS2_TYPES_INCLUDED

#define _INCLUDE_LONGLONG
#include<Rcpp.h>
#include<RcppGSL.h>
#include <stdlib.h>  //malloc
#include <gsl/gsl_rng.h>

#include <string.h>
#include <iostream>
#include <math.h>
#include <inttypes.h>


using namespace std;

typedef uint32_t      bitvector;
typedef unsigned long long longbitvector;
typedef unsigned int       uint;
typedef unsigned long long ulonglong;

typedef bitvector bvrand();

#define BV_BYTES   (sizeof(bitvector))
#define BV_BITS    (8*BV_BYTES)
#define BV_MAXBIT  (((bitvector)1 << (BV_BITS-1)))
#define BV_MAX     ((bitvector)-1)
#define BV_RES     (1.0/((double)BV_MAX+1.0))
#define LBV_BYTES  (sizeof(longbitvector))
#define LBV_BITS   (8*LBV_BYTES)
#define LBV_MAXBIT (((longbitvector)1 << (LBV_BITS-1)))
#define LBV_MAX    ((longbitvector)-1)
#define LBV_RES    (1.0/((double)LBV_MAX+1.0))

typedef enum _dgt
{
  dgt_Sobol,
  dgt_SpecialNiederreiter,
  dgt_NiederreiterXing,
  dgt_ShiftNet
} dgt;



typedef void (*SimInitPF)(gsl_rng*, int, int,int, double*, double*);

typedef double* (*ParamTransitionPF)(double*, int, int);

typedef void (*TransitionPF)(double*, double*, int, int, int, double*, double*);

typedef void (*SimTransitionPF)(gsl_rng*, double*, int, int, int, int, double*, double*, double*);

typedef void (*PotentialPF)(double*, double*, double*, int, int, int, int, double*, double*);

typedef void (*ResamplingPF)(gsl_rng*, double*, int, int, int, double*, double*);


#endif
