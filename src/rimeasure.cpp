/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: measurements of radical inverses point sets                    //
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
#include "Measurements.h"

#include <stdlib.h>
#include <iostream>

#include "RadicalInversesBase2.h"

using namespace std;
typedef enum _measure
{
  measure_L2,
  measure_Star,
  measure_MinDist,
  measure_MinDistTorus
} measure;

typedef enum _nt
{
  nt_Sobol,
  nt_Hammersley,
  nt_LarcherPillichshammer
} nt;

class Options
{
public:
  
  measure dtype;         // type of discrepancy
  nt      ntype;         // type of net
  bool    ds;            // digit scrambling
  int     n;             // number of points in net
  bool    sequence;      // calculate discrepancy up to n
  
  Options(int argc, char **argv)
    {
      dtype = measure_L2;
      ntype = nt_Sobol;
      ds = false;
      n = 16;
      sequence = false;
      for (int i = 1; i < argc; ++i)
	if (!strcmp("-measure", argv[i]))
	  if (!strcmp("L2", argv[++i]))
	    dtype = measure_L2;
	  else if (!strcmp("Star", argv[i]))
	    dtype = measure_Star;
	  else if (!strcmp("MinDist", argv[i]))
	    dtype = measure_MinDist;
	  else if (!strcmp("MinDistTorus", argv[i]))
	    dtype = measure_MinDistTorus;
	  else
	    {
	      cerr << "Unknown measurement type: " << argv[i] << endl;
	      cerr << "Possible types: L2, Star, MinDist, MinDistTorus" << endl;
	      exit(-1);
	    }
	else if (!strcmp("-nt", argv[i]))
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
	else if (!strcmp("-ds", argv[i]))
	  ds = true;
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
	    cout << "Usage: " << argv[0] << " [-measure measurementtype] [-nt nettype] [-ds] [-n points] [-s]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-measure measurementtype] [-nt nettype] [-ds] [-n points] [-s]" << endl;
	    exit(-1);
	  }
    }
};



