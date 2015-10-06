/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Digital nets in base 2                                                  //
//                                                                         //
//                                                                         //
// The method "DigitalNetGenerator::GettParameter" has been adapted from   //
// code kindly provided by Wolfgang Schmid. No more table lookup memory is //
// required anymore.                                                       //
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
#include <fstream>
#include <math.h>
#include "DigitalNetsBase2.h"


#include "DB_NX.h"
#include "DB_ShiftNets.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////

void GeneratorMatrices::Initialize(const uint s)
{
  if (data != NULL)
    delete []data;
  dimension = s;
  data = new bitvector[dimension*BV_BITS];
}

void GeneratorMatrices::Initialize(const uint s, const bitvector *_data)
{
  if (data != NULL)
    delete []data;
  dimension = s;
  data = new bitvector[dimension*BV_BITS];
  memcpy(data, _data, dimension*BV_BITS*BV_BYTES);
}

int GeneratorMatrices::Load(const char *Filename)
{
  char      i;
  ifstream  File(Filename, ios::binary);

  if (data != NULL)
    delete []data;
  File.read(&i, 1);
  if (i != (char)BV_BITS)
    {
      cerr << "Error: Can't load generator matrices" << endl
	   << "       Generator matrices have wrong size" << endl;
      return 1;
    }
  File.read(&i, 1);
  dimension = (uint)i;
  data = new bitvector[dimension*BV_BITS];
  File.read((char *)data, dimension*BV_BITS*BV_BYTES);
  File.close();
  return 0;
}

void GeneratorMatrices::Save(const char *Filename) const
{
  char      i;
  ofstream  File(Filename, ios::binary);

  i = (char)BV_BITS;
  File.write(&i, 1);
  i = (char)dimension;
  File.write(&i, 1);
  File.write((char *)data, dimension*BV_BITS*BV_BYTES);
  File.close();
}

void GeneratorMatrices::Print
(ostream &s, const int m = BV_BITS, const bool binary = true) const
{
  bitvector pf, bv;
  uint      i, j, k;

  for (i = 0; i < dimension; ++i)
    {
      for(j = 0; j < m; ++j)
	{
	  bv = GetRow(i, j);
	  if (binary)
	    {
	      for (k = m, pf = BV_MAXBIT; k > 0; --k, pf >>= 1)
		if (bv & pf)
		  s << "1 ";
		else
		  s << "0 ";
	      s << endl;
	    }
	  else
	    s << (bv >> (BV_BITS - m)) << " ";
	}
      s << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

void DigitalNetGenerator::InitSobol(const uint s)
{
  uint       i, j;
  Polynomial p_i, p, q;
  uint       deg;
  bitvector  bv;

  gm.Initialize(s);
  for (bv = BV_MAXBIT, j = 0; bv; bv >>= 1, ++j)
    gm.SetColumn(0, j, bv);
  p_i = 1;
  for (i = 1; i < s; ++i)
    {
      for (p_i += 2; !p_i.IsPrimitive(); p_i += 2);
      deg = p_i.Degree();
      p = 1;
      for (j = 0; j < BV_BITS; ++j)
	{
	  if (j % deg == 0)
	    {
	      p *= p_i;
	      q = ((longbitvector)BV_MAXBIT << deg);
	    }
	  else
  	    q >>= 1;
	  gm.SetRow(i, j, (bitvector)(q / p));
	}
    }
}

void DigitalNetGenerator::InitSpecialNiederreiter(const uint s)
{
  uint       i, j;
  Polynomial p_i, p, q;
  uint       deg;
  bitvector  bv;

  gm.Initialize(s);
  for (bv = BV_MAXBIT, j = 0; bv; bv >>= 1, ++j)
    gm.SetColumn(0, j, bv);
  p_i = 1;
  for (i = 1; i < s; ++i)
    {
      for (p_i += 2; !p_i.IsIrreducible(); p_i += 2);
      deg = p_i.Degree();
      p = 1;
      for (j = 0; j < BV_BITS; ++j)
	{
	  if (j % deg == 0)
	    {
	      p *= p_i;
	      q = ((longbitvector)BV_MAXBIT << deg);
	    }
	  else
  	    q >>= 1;
	  gm.SetRow(i, j, (bitvector)(q / p));
	}
    }
}

void DigitalNetGenerator::InitNiederreiterXing(const uint s)
{
  if ((s < 4) || (s > 32))
    {
      cerr << "Error: NiederreiterXing can't be initialized" << endl
	   << "       Only dimensions 4 to 32 are supported" << endl;
      exit(-1);
    }
  gm.Initialize(s, _NiederreiterXingMatrices + BV_BITS * (s * (s-1) / 2 - 6));
}

bitvector Reverse(bitvector bv)
{
  bitvector bit, res;

  res = 0;
  for (bit = BV_MAXBIT; bv; bit >>= 1, bv >>= 1)
    if (bv & 1)
      res |= bit;
  return res;
}

void DigitalNetGenerator::InitShiftNet(const uint s)
{
  uint i, j, k, m;

  if ((s < 3) || (s > 39))
    {
      cerr << "Error: ShiftNet can't be initialized" << endl
           << "       Only dimensions 3 to 39 are supported" << endl;
      exit(-1);
    }
  gm.Initialize(s);
  if (s >= 23)
    m = 23;
  else
    m = s;
  k = ShiftNet[m][0];
  for (i = 0; i < k; ++i)
    gm.SetRow(0, i, Reverse(ShiftNet[m][i+1]));
  for (; i < BV_BITS; ++i)
    gm.SetRow(0, i, 0);
  for (j = 1; j < s; ++j)
    {
      for (i = 0; i < s; ++i)
	gm.SetColumn(j, i, gm.GetColumn(j-1, (i+1)%s));
      for (; i < BV_BITS; ++i)
	gm.SetColumn(j, i, 0);
    }
}

DigitalNetGenerator::DigitalNetGenerator(const uint s)
{
  gm.Initialize(s);
  state = new bitvector[s];
  qrn = new double[s];
  Reset();
}

DigitalNetGenerator::DigitalNetGenerator(const char *Filename)
{
  gm.Load(Filename);
  state = new bitvector[gm.GetDimension()];
  qrn = new double[gm.GetDimension()];
  Reset();
}

DigitalNetGenerator::DigitalNetGenerator(const dgt type, const uint s)
{
  switch (type)
    {
    case dgt_Sobol:
      InitSobol(s);
      break;
    case dgt_SpecialNiederreiter:
      InitSpecialNiederreiter(s);
      break;
    case dgt_NiederreiterXing:
      InitNiederreiterXing(s);
      break;
    case dgt_ShiftNet:
      InitShiftNet(s);
      break;
    }
  state = new bitvector[s];
  qrn = new double[s];
  Reset();
}

DigitalNetGenerator::~DigitalNetGenerator()
{
  delete []state;
  delete []qrn;
}

uint DigitalNetGenerator::GettParameter(const uint s, const uint m) const
{
  bitvector mask;
  bitvector **vectors;
  uint      i, j, l, u;
  bitvector r, lindep, vek[32];
  uint      dd[32];

  mask = (BV_MAX << (32 - m));
  vectors = new bitvector*[s];
  for (j = 0; j < s; ++j)
    {
      vectors[j] = new bitvector[m];
      for (i = 0; i < m; ++i)
	vectors[j][i] = gm.GetRow(j, i) & mask;
    }

  for (u = 1; u <= m; ++u)
    {
      dd[0] = u;
      for(j = 1; j < s; ++j)
        dd[j] = 0;
      while (j != s-1)
        {
          l = 0;
          lindep = 0;
          for (i = 0; i < s; ++i)
            if (dd[i] > 0)
              {
                lindep ^= vectors[i][dd[i]-1];
                for (j = 0; j < dd[i]-1; ++j)
                  {
                    vek[l] = vectors[i][j];
                    ++l;
                  };
              };
          if (lindep == 0)
	    {
	      for (j = 0; j < s; ++j)
		delete []vectors[j];
	      delete []vectors;
	      return m+1-u;
	    }
          for (r = 1; r < ((bitvector)1<<l); ++r)
            {
              lindep ^= vek[GrayCodeChange(r)];
              if (lindep == 0)
		{
		  for (j = 0; j < s; ++j)
		    delete []vectors[j];
		  delete []vectors;
		  return m+1-u;
		}
            };
	  for (j = 0; j < s-1; ++j)
	    if (dd[j] > 0)
	      {
		++dd[j+1];
		dd[0] = dd[j]-1;
		if (j > 0)
		  dd[j] = 0;
		break;
	      }
        }
    };
  for (j = 0; j < s; ++j)
    delete []vectors[j];
  delete []vectors;
  return 0;
}

void DigitalNetGenerator::UndoGray()
{
  uint      i;
  bitvector *columns1, *columns2;

  columns1 = gm.GetColumns(0);
  columns2 = gm.GetColumns(1);
  for (i = 0;
       i < (gm.GetDimension() * (BV_BITS-1));
       ++i, ++columns1, ++columns2)
    *columns2 ^= *columns1;
}

/////////////////////////////////////////////////////////////////////////////

DigitScrambled::DigitScrambled(const char *Filename, bvrand &rr)
  : DigitalNetGenerator(Filename)
{
  r = &rr;
}

DigitScrambled::DigitScrambled(const dgt type, const uint s, bvrand &rr)
  : DigitalNetGenerator(type, s)
{
  r = &rr;
}

/////////////////////////////////////////////////////////////////////////////

LinearScrambled::LinearScrambled(const char *Filename, bvrand &rr)
  : DigitalNetGenerator(Filename)
{
  r = &rr;
  originalgm = gm;
}

LinearScrambled::LinearScrambled(const dgt type, const uint s, bvrand &rr)
  : DigitalNetGenerator(type, s)
{
  r = &rr;
  originalgm = gm;
}

/////////////////////////////////////////////////////////////////////////////

Scrambled::Scrambled(const char *Filename, const uint nn, bvrand &rr)
  : DigitalNetGenerator(Filename)
{
  n = nn;
  data = new bitvector[gm.GetDimension()*n];
  r = &rr;
}

Scrambled::Scrambled(const dgt type, const uint s, const uint nn, bvrand &rr)
  : DigitalNetGenerator(type, s)
{
  n = nn;
  data = new bitvector[s*n];
  r = &rr;
}

Scrambled::~Scrambled()
{
  delete []data;
}

/////////////////////////////////////////////////////////////////////////////
