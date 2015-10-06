/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: Measurements                                                   //
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

#include "Measurements.h"

typedef enum _measure
{
  measure_L2,
  measure_Star,
  measure_MinDist,
  measure_MinDistTorus
} measure;

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
  
  measure dtype;         // type of discrepancy
  dgt     gtype;         // type of generator
  rt      rtype;         // randomization type
  int     d;             // dimension
  int     n;             // number of points in net
  bool    sequence;      // calculate discrepancy up to n
  
  Options(int argc, char **argv)
    {
      dtype = measure_L2;
      gtype = dgt_Sobol;
      rtype = rt_None;
      d = 2;
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
	else if (!strcmp("-dgt", argv[i]))
	  if (!strcmp("Sobol", argv[++i]))
	    gtype = dgt_Sobol;
	  else if (!strcmp("SpecialNiederreiter", argv[i]))
	    gtype = dgt_SpecialNiederreiter;
	  else if (!strcmp("NiederreiterXing", argv[i]))
	    gtype = dgt_NiederreiterXing;
	  else if (!strcmp("ShiftNet", argv[i]))
	    gtype = dgt_ShiftNet;
	  else
	    {
	      cerr << "Unknown generator type: " << argv[i] << endl;
	      cerr << "Possible types: Sobol, SpecialNiederreiter, NiederreiterXing, ShiftNet" << endl;
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
	    d = atoi(argv[++i]);
	    if (d < 1)
	      {
		cerr << "Number of dimensions has to be positive" << endl;
		exit(-1);
	      }
	  }
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
	    cout << "Usage: " << argv[0] << " [-measure measurementtype] [-dgt generatortype] [-rt randomizationtype] [-d dimension] [-n points] [-s]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-measure measurementtype] [-dgt generatortype] [-rt randomizationtype] [-d dimension] [-n points] [-s]" << endl;
	    exit(-1);
	  }
    }
};

/*
int main(int argc, char **argv)
{
  DigitalNetGenerator *dng;
  DigitScrambled      *ds;
  LinearScrambled     *ls;
  Scrambled           *os;
  int                 i, j;
  double              *p;

  Options o(argc, argv);
  
  p = AllocPointMemory(o.d, o.n);

  switch (o.rtype)
    {
    case rt_None:
      dng = new DigitalNetGenerator(o.gtype, o.d);
      for (i = 0; i < o.n; ++i)
	{
	  for (j = 0; j < o.d; ++j)
	    p[i*o.d+j] = (*dng)[j];
	  ++(*dng);
	}
      delete dng;
      break;
    case rt_DigitScrambled:
      ds = new DigitScrambled(o.gtype, o.d, rng32);
      ds->Randomize();
      for (i = 0; i < o.n; ++i)
	{
	  for (j = 0; j < o.d; ++j)
	    p[i*o.d+j] = (*ds)[j];
	  ++(*ds);
	}
      delete ds;
      break;
    case rt_LinearScrambled:
      ls = new LinearScrambled(o.gtype, o.d, rng32);
      ls->Randomize();
      for (i = 0; i < o.n; ++i)
	{
	  for (j = 0; j < o.d; ++j)
	    p[i*o.d+j] = (*ls)[j];
	  ++(*ls);
	}
      delete ls;
      break;
    case rt_Scrambled:
      os = new Scrambled(o.gtype, o.d, o.n, rng32);
      os->Randomize();
      for (i = 0; i < o.n; ++i)
	{
	  for (j = 0; j < o.d; ++j)
	    p[i*o.d+j] = (*os)[j];
	  ++(*os);
	}
      delete os;
      break;
    }

  switch (o.dtype)
    {
    case measure_L2:
      if (o.sequence)
	for (i = 1; i <= o.n; ++i)
	  cout << i << " " << L2Discrepancy(o.d, i, p) << endl;
      else
	cout << o.n << ": " << L2Discrepancy(o.d, o.n, p) << endl;
      break;
    case measure_Star:
      if (o.sequence)
	for (i = 1; i <= o.n; ++i)
	  cout << i << ": " << StarDiscrepancy(o.d, i, p) << endl;
      else
	cout << o.n << ": " << StarDiscrepancy(o.d, o.n, p) << endl;
      break;
    case measure_MinDist:
      if (o.sequence)
	for (i = 1; i <= o.n; ++i)
	  cout << i << ": " << MinimumDistance(o.d, i, p) << endl;
      else
	cout << o.n << ": " << MinimumDistance(o.d, o.n, p) << endl;
      break;
    case measure_MinDistTorus:
      if (o.sequence)
	for (i = 1; i <= o.n; ++i)
	  cout << i << ": " << MinimumDistanceTorus(o.d, i, p) << endl;
      else
	cout << o.n << ": " << MinimumDistanceTorus(o.d, o.n, p) << endl;
      break;
    }
  
  FreePointMemory(p, o.d);
  
  return 0;
}
*/
