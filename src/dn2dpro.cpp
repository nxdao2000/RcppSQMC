/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: 2D projection of a point set                                   //
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
#include "DigitalNetsBase2.h"

#include <stdlib.h>
#include <stdio.h>
#include <fstream>


#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
// [[Rcpp::export]]
uint32_t rng32()
{
  return (uint32_t)(drand48()*((double)BV_MAX+1.0));
}

typedef enum _rt
{
  rt_None,
  rt_DigitScrambled,
  rt_LinearScrambled,
  rt_Scrambled
} rt;

class Options
{
public:

  dgt  gtype;         // type of generator
  rt   rtype;         // randomization type
  int  d;             // net dimension
  int  d1, d2;        // projection dimensions
  int  n;             // number of points in net
  char filename[255]; // output filename

  Options(int argc, char **argv)
    {
      gtype = dgt_Sobol;
      rtype = rt_None;
      d = 2;
      d1 = 0;
      d2 = 1;
      n = 16;
      sprintf(filename, "dn.tex");
      for (int i = 1; i < argc; ++i)
	if (!strcmp("-dgt", argv[i]))
	  if (!strcmp("Sobol", argv[++i]))
	    gtype = dgt_Sobol;
	  else if (!strcmp("SpecialNiederreiter", argv[i]))
	    gtype = dgt_SpecialNiederreiter;
	  else if (!strcmp("NiederreiterXing", argv[i]))
	    {
	      gtype = dgt_NiederreiterXing;
	      d = atoi(argv[++i]);
	      if ((d < 4) || (d > 32))
		{
		  cerr << "Only dimensions 4 to 32 are supported for NiederreiterXing" << endl;
		  exit(-1);
		}
	    }
	  else if (!strcmp("ShiftNet", argv[i]))
	    {
	      gtype = dgt_ShiftNet;
	      d = atoi(argv[++i]);
	      if ((d < 3) || (d > 39))
		{
		  cerr << "Only dimensions 3 to 39 are supported for ShiftNet" << endl;
		  exit(-1);
		}
	    }
	  else
	    {
	      cerr << "Unknown generator type: " << argv[i] << endl;
	      cerr << "Possible types: Sobol, SpecialNiederreiter, NiederreiterXing dimension, ShiftNet dimension" << endl;
	      exit(-1);
	    }
	else if (!strcmp("-rt", argv[i]))
	  if (!strcmp("None", argv[++i]))
	    rtype = rt_None;
	  else if (!strcmp("DigitScrambled", argv[i]))
	    rtype = rt_DigitScrambled;
	  else if (!strcmp("LinearScrambled", argv[i]))
	    rtype = rt_LinearScrambled;
	  else if (!strcmp("Scrambled", argv[i]))
	    rtype = rt_Scrambled;
	  else
	    {
	      cerr << "Unknown randomization type: " << argv[i] << endl;
	      cerr << "Possible types: None, DigitScrambled, LinearScrambled, Scrambled" << endl;
	      exit(-1);
	    }
	else if (!strcmp("-d", argv[i]))
	  {
	    d1 = atoi(argv[++i]);
	    d2 = atoi(argv[++i]);
	    if ((d1 < 0) || (d2 < 0))
	      {
		cerr << "Projection dimensions have to be non-negative" << endl;
		exit(-1);
	      }
	    d = MAX(d, d1+1);
	    d = MAX(d, d2+1);
	  }
	else if (!strcmp("-n", argv[i]))
	  {
	    n = atoi(argv[++i]);
	    if (n < 1)
	      {
		cerr << "The net should have at least one point" << endl;
		exit(-1);
	      }
	    if (n > 16384)
	      {
		cerr << "Nets with maximal 16384 points can be printed due to LaTeX" << endl;
		exit(-1);
	      }
	  }
	else if (!strcmp("-o", argv[i]))
	  sprintf(filename, "%s.tex", argv[++i]);
	else if (!strcmp("-help", argv[i]))
	  {
	    cout << "Usage: " << argv[0] << " [-dgt generatortype] [-rt randomizationtype] [-d dimension1 dimension2] [-n points] [-o filename]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-dgt generatortype] [-rt randomizationtype] [-d dimension1 dimension2] [-n points] [-o filename]" << endl;
	    exit(-1);
	  }
    }
};


