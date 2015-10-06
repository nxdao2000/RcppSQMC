/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Example: t-parameter calculation                                        //
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

#include <stdlib.h>

#include "DigitalNetsBase2.h"


class Options
{
public:
  
  dgt gtype;         // type of generator
  uint d;            // dimension
  uint m;            // maximum precision
  
  Options(int argc, char **argv)
    {
      gtype = dgt_Sobol;
      d = 2;
      m = 8;
      for (int i = 1; i < argc; ++i)
	if (!strcmp("-dgt", argv[i]))
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
	else if (!strcmp("-d", argv[i]))
	  {
	    d = atoi(argv[++i]);
	    if (d < 1)
	      {
		cerr << "Number of dimensions has to be positive" << endl;
		exit(-1);
	      }
	  }
	else if (!strcmp("-m", argv[i]))
	  {
	    m = atoi(argv[++i]);
	    if (m < 1)
	      {
		cerr << "Precision has to be positive" << endl;
		exit(-1);
	      }
	    else if (m > 8 * BV_BYTES)
	      {
		cerr << "Maximum precision is " << 8*BV_BYTES << endl;
		exit(-1);
	      }
	  }
	else if (!strcmp("-help", argv[i]))
	  {
	    cout << "Usage: " << argv[0] << " [-dgt generatortype] [-d dimension] [-m precision]" << endl;
	    exit(0);
	  }
	else
	  {
	    cerr << "Unknown option: " << argv[i] << endl;
	    cerr << "Usage: " << argv[0] << " [-dgt generatortype] [-d dimension] [-m precision]" << endl;
	    exit(-1);
	  }
    }
};

int main(int argc, char **argv)
{
  DigitalNetGenerator *dng;
  uint                i;

  Options o(argc, argv);

  dng = new DigitalNetGenerator(o.gtype, o.d);

  cout << "t-Parameters depending on precision: " << endl;
  for (i = 1; i <= o.m; ++i)
    cout << i << ": " << dng->GettParameter(o.d, i) << endl;
}
