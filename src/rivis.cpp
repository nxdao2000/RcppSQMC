/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: visualization of radical inverse point sets                    //
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
#include "RadicalInversesBase2.h"
#include "Measurements.h"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

using namespace std;

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
  char filename[255]; // output filename

  Options(int argc, char **argv)
    {
      ntype = nt_Sobol;
      rds = false;
      n = 16;
      sprintf(filename, "ri.tex");
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
	    cout << "Usage: " << argv[0] << " [-nt nettype] [-rds] [-n points] [-o filename]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-nt nettype] [-rds] [-n points] [-o filename]" << endl;
	    exit(-1);
	  }
    }
};



