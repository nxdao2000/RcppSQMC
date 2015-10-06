/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: cap discrepancy of radical inverse point sets                  //
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
#include<Rcpp.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "RadicalInversesBase2.h"
#include "Measurements.h"
using namespace std;
uint rng32a()
{
  return (uint)(drand48()*((double)((uint)-1)+1.0));
}

typedef enum _nt
{
  nt_Sobol,
  nt_Hammersley,
  nt_LarcherPillichshammer
} nt;

class Options
{
public:
  
  nt   ntype;         // type of net
  bool rds;           // random digit scrambling
  int  n;             // number of points in net
  bool sequence;      // calculate discrepancy up to n
  
  Options(int argc, char **argv)
    {
      ntype = nt_Sobol;
      rds = false;
      n = 16;
      sequence = false;
      for (int i = 1; i < argc; ++i)
	if (!strcmp("-nt", argv[i]))
	  if (!strcmp("Sobol", argv[++i]))
	    ntype = nt_Sobol;
	  else if (!strcmp("Hammersley", argv[i]))
	    ntype = nt_Hammersley;
	  else if (!strcmp("LarcherPillichshammer", argv[i]))
	    ntype = nt_LarcherPillichshammer;
	  else
	    {
	      cerr << "Unknown generator type: " << argv[i] << endl;
	      cerr << "Possible types: Sobol, Hammersley, LarcherPillichshammer" << endl;
	      exit(-1);
	    }
	else if (!strcmp("-rds", argv[i]))
	  rds = true;
	else if (!strcmp("-n", argv[i]))
	  {
	    n = atoi(argv[++i]);
	    if (n < 1)
	      {
		cerr << "The net should have at least one point" << endl;
		exit(-1);
	      }
	  }
	else if (!strcmp("-s", argv[i]))
	  sequence = true;
	else if (!strcmp("-help", argv[i]))
	  {
	    cout << "Usage: " << argv[0] << " [-nt nettype] [-rds] [-n points] [-s]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-nt nettype] [-rds] [-n points] [-s]" << endl;
	    exit(-1);
	  }
    }
};

void InitSpherePoints(double *sp, const double *p, const uint n)
{
  double phi, costheta, sintheta;
  int    i;
  
  for (i = 0; i < n; ++i)
    {
      costheta = 2.0 * p[2*i] - 1.0;
      sintheta = sqrt(1.0 - costheta*costheta);
      phi = 2.0 * M_PI * p[2*i+1];
      sp[3*i] = sintheta * cos(phi);
      sp[3*i+1] = sintheta * sin(phi);
      sp[3*i+2] = costheta;
    }
}


