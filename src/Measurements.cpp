/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// SAMPLE PACKage                                                          //
//                                                                         //
// Measurements                                                            //
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
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

typedef unsigned int       uint;
typedef unsigned long long ulonglong;

/////////////////////////////////////////////////////////////////////////////

double *AllocPointMemory(const uint s, const uint n)
{
  double *pp;
  uint   j;
  
  pp = new double[(n+1)*s];
  for (j = 0; j < s; ++j)
    pp[j] = 1.0;
  return pp+s;
}

void FreePointMemory(double *p, const uint s)
{
  delete [](p-s);
}

/////////////////////////////////////////////////////////////////////////////

#define NODESIZE  4
#define STACKSIZE 256

class BSPTree
{
private:
  
  static double *p;
  static uint   s;
  static uint   *pi;
  int           split;
  double        x;
  uint          ln, rn;
  BSPTree       *Left, *Right;
  
  BSPTree
  (const uint start, const uint _n, const uint depth, const double *_x)
    {
      double *__x = new double[s];
      uint   h, l, r;

      if (_n < NODESIZE)
	{
	  split = -1;
	  ln = start;
	  rn = _n;
	  delete []__x;
	  return;
	}
      split = depth % s;
      x = _x[split];
      l = start;
      r = start + _n;
      while (l != r)
	if (p[pi[l]*s+split] < x)
	  ++l;
	else
	  {
	    --r;
	    h = pi[l];
	    pi[l] = pi[r];
	    pi[r] = h;
	  }
      ln = l - start;
      rn = start + _n - r;

      for (h = 0; h < s; ++h)
	__x[h] = _x[h];
      __x[split] = x - 1.0 / (double)(2 + depth / s);
      Left = new BSPTree(start, ln, depth+1, __x);
      for (h = 0; h < s; ++h)
	__x[h] = _x[h];
      __x[split] = x + 1.0 / (double)(2 + depth / s);
      Right = new BSPTree(l, rn, depth+1, __x);
      delete []__x;
    }

public:

  BSPTree(const uint _s, const uint _n, double *_p)
    {
      double *_x = new double[_s];
      uint   h, l, r;
      
      p = _p;
      s = _s;
      pi = new uint[_n+1];
      for (h = 0; h <= _n; ++h)
	pi[h] = h;
      if (_n < NODESIZE)
	{
	  split = -1;
	  ln = 0;
	  rn = _n;
	  delete []_x;
	  return;
	}
      split = 0;
      x = 0.5;
      l = 0;
      r = _n;
      while (l != r)
	if (p[pi[l]*s] < x)
	  ++l;
	else
	  {
	    --r;
	    h = pi[l];
	    pi[l] = pi[r];
	    pi[r] = h;
	  }
      ln = l;
      rn = _n - r;
      _x[0] = 0.25;
      for (h = 1; h < s; ++h)
	_x[h] = 0.5;
      Left = new BSPTree(0, ln, 1, _x);
      _x[0] = 0.75;
      for (h = 1; h < s; ++h)
	_x[h] = 0.5;
      Right = new BSPTree(l, rn, 1, _x);
      delete []_x;
    }
    
  ~BSPTree()
    {
      if (pi != NULL)
        {
          delete []pi;
          pi = NULL;
        }
      if (split != -1)
	{
	  if (pi != NULL)
	    {
	      delete []pi;
	      pi = NULL;
	    }
	  delete Left;
	  delete Right;
	}
    }
  
  void PointsInBox(uint &cin, uint &cb, const double *_x)
    {
      static BSPTree *Next[STACKSIZE];
      static BSPTree *Node;
      static uint    isin[STACKSIZE];
      static uint    aktisin;
      static uint    StackPtr;
      static uint    h, j, ak;

      cin = cb = 0;
      Node = this;
      StackPtr = 0;
      aktisin = (1 << s) - 1;
      for (h = 0; h < s; ++h)
	if (_x[h] == 1.0)
	  aktisin -= (1 << h);
      while (true)
	{
	  if (Node->split == -1)
	    for (h = Node->ln; h < Node->ln+Node->rn; ++h)
	      {
		++cin;
		ak = aktisin;
		for (j = 0; j < s; ++j, ak = (ak >> 1))
		  {
		    if (_x[j] < p[pi[h]*s+j])
		      {
			--cin;
			j = s;
		      }
		    else if (_x[j] == p[pi[h]*s+j])
		      {
			--cin;
			++cb;
			for (++j; j < s; ++j, ak = (ak >> 1))
			  {
			    if (_x[j] < p[pi[h]*s+j])
			      {
				--cb;
				j = s;
			      }
			  }
		      }
		  }
	      }
	  else if (_x[Node->split] < Node->x)
	    {
	      isin[StackPtr] = aktisin;
	      Next[StackPtr++] = Node->Left;
	    }
	  else if (_x[Node->split] == Node->x)
	    {
	      isin[StackPtr] = aktisin;
	      Next[StackPtr++] = Node->Left;
	      isin[StackPtr] = aktisin;
	      Next[StackPtr++] = Node->Right;
	    }
	  else
	    {
	      isin[StackPtr] = aktisin;
	      Next[StackPtr++] = Node->Right;
	      isin[StackPtr] = aktisin & ~(1 << Node->split);
	      if (isin[StackPtr] == 0)
		cin += Node->ln;
	      else
		Next[StackPtr++] = Node->Left;
	    }
	  if (StackPtr == 0)
	    return;
	  Node = Next[--StackPtr];
	  aktisin = isin[StackPtr];
	}
    }
};

double *BSPTree::p;
uint   BSPTree::s;
uint   *BSPTree::pi = NULL;

double StarDiscrepancy(const uint s, const uint n, double *_p)
{
  ulonglong i, ns, k;
  uint      j;
  uint      c, b;
  double    dis = 0.0, disn1, disn2;
  double    *p = _p-s;
  double    *x  = new double[s];
  BSPTree   bsp(s, n, _p);
  
  ns = n+1;
  for (j = 1; j < s; ++j)
    ns *= (n+1);
  for (i = 0; i < ns; ++i)
    {
      for (j = 0, k = 1; j < s; ++j, k *= (n+1))
	x[j] = p[((i/k)%(n+1))*s+j];
      bsp.PointsInBox(c, b, x);
      disn1 = 1.0;
      for (j = 0; j < s; ++j)
	disn1 *= x[j];
      disn2 = disn1;
      disn1 -= (double)c / (double)n;
      disn1 = fabs(disn1);
      if (dis < disn1)
	dis = disn1;
      disn2 -= (double)(c + b) / (double)n;
      disn2 = fabs(disn2);
      if (dis < disn2)
	dis = disn2;
    }
  delete []x;
  
  return dis;
}

/////////////////////////////////////////////////////////////////////////////

#define THRESHOLD 6

double L2DiscrepancyT2(const uint s, const uint n, double *p)
{
  double T2 = 0.0;
  double d;
  uint   i, k;
  
  for (i = 0; i < n; ++i)
    {
      d = 1.0;
      for (k = 0; k < s; ++k)
	d *= (1.0 - p[i*s+k]*p[i*s+k]);
      T2 -= d;
    }
  return T2 / (pow((double)2, (double)(s-1)) * (double)n);
}

double L2DiscrepancyT3
(const uint m, const uint n, uint *A, uint *B, double *Aw, double *Bw,
 const uint sm, const double a, const double b, const uint s, double *p)
{
  uint   i, j, k;
  double hd1, hd2;

  if ((m == 0) || (n == 0))
    return 0.0;
  else if (sm == 0)
    {
      hd1 = hd2 = 0.0;
      for (i = 0; i < m; ++i)
	hd1 += Aw[i];
      for (i = 0; i < n; ++i)
	hd2 += Bw[i];
      return hd1 * hd2;
    }
  else if ((m < THRESHOLD * sm) || (n < THRESHOLD * sm) || (sm < 2))
    {
      hd2 = 0.0;
      for (i = 0; i < m; ++i)
	for (j = 0; j < n; ++j)
	  {
	    hd1 = (1.0 - MAX(p[A[i]*s], p[B[j]*s]));
	    for (k = 1; k < sm; ++k)
	      hd1 *= (1.0 - MAX(p[A[i]*s+k], p[B[j]*s+k]));
	    hd2 += Aw[i] * Bw[j] * hd1;
	  }
      return hd2;
    }
  else
    {
      double theta = 0.5 * (a+b);
      uint   *A2, *B2;
      double *Aw2, *Bw2;
      uint   al, ar, bl, br;
      uint   smax = sm-1;
      uint   hui;
      
      al = 0;
      ar = m;
      Aw2 = new double[m];
      memcpy(Aw2, Aw, m*sizeof(double));
      while (al != ar)
	{
	  if (p[A[al]*s+smax] < theta)
	    ++al;
	  else
	    {
	      --ar;
	      hui = A[al];
	      A[al] = A[ar];
	      A[ar] = hui;
	      hd1 = Aw[al];
	      Aw[al] = Aw[ar];
	      Aw[ar] = hd1;
	      hd1 = Aw2[al];
	      Aw2[al] = Aw2[ar];
	      Aw2[ar] = hd1 * (1.0 - p[hui*s+smax]);
	    }
	}
      A2 = new uint[m];
      memcpy(A2, A, m*sizeof(uint));
      bl = 0;
      br = n;
      Bw2 = new double[n];
      memcpy(Bw2, Bw, n*sizeof(double));
      while (bl != br)
	{
	  if (p[B[bl]*s+smax] < theta)
	    ++bl;
	  else
	    {
	      --br;
	      hui = B[bl];
	      B[bl] = B[br];
	      B[br] = hui;
	      hd1 = Bw[bl];
	      Bw[bl] = Bw[br];
	      Bw[br] = hd1;
	      hd1 = Bw2[bl];
	      Bw2[bl] = Bw2[br];
	      Bw2[br] = hd1 * (1.0 - p[hui*s+smax]);
	    }
	}
      if ((al == 0) || (ar == 0) || (bl == 0) || (br == 0))
        {
          delete []Aw2;
          delete []Bw2;
          delete []A2;
          hd2 = 0.0;
          for (i = 0; i < m; ++i)
            for (j = 0; j < n; ++j)
              {
                hd1 = (1.0 - MAX(p[A[i]*s], p[B[j]*s]));
                for (k = 1; k < sm; ++k)
                  hd1 *= (1.0 - MAX(p[A[i]*s+k], p[B[j]*s+k]));
                hd2 += Aw[i] * Bw[j] * hd1;
              }
          return hd2;
        }
      B2 = new uint[n];
      memcpy(B2, B, n*sizeof(uint));
      hd2 =
	L2DiscrepancyT3(al, n-bl, A2, B2+bl, Aw2, Bw2+bl, smax, 0, 1, s, p)
	+ L2DiscrepancyT3(m-al, bl, A2+al, B2, Aw2+al, Bw2, smax, 0, 1, s, p);
      delete []Aw2;
      delete []Bw2;
      delete []A2;
      delete []B2;
      hd1 =
	L2DiscrepancyT3(al, bl, A, B, Aw, Bw, sm, a, theta, s, p)
	+ L2DiscrepancyT3(m-al, n-bl, A+al, B+bl,
			  Aw+al, Bw+bl, sm, theta, b, s, p);
      return hd1 + hd2;
    }
}

double L2Discrepancy(const uint s, const uint n, double *p)
{
  double T3 = 0.0;
  uint   *Aindex = new uint[n];
  uint   *Bindex = new uint[n];
  double *Aweights = new double[n];
  double *Bweights = new double[n];
  uint   i;
  
  for (i = 0; i < n; ++i)
    {
      Aindex[i] = Bindex[i] = i;
      Aweights[i] = Bweights[i] = 1.0 / (double)n;
    }
  T3 =
    L2DiscrepancyT3(n, n, Aindex, Bindex, Aweights, Bweights,
		    s, 0.0, 1.0, s, p);
  delete []Aindex;
  delete []Bindex;
  delete []Aweights;
  delete []Bweights;
  return sqrt(pow((double)3, -(double)s) + L2DiscrepancyT2(s, n, p) + T3);
}

/////////////////////////////////////////////////////////////////////////////

inline double EuclidianDistance
(const uint s, const double *p1, const double *p2)
{
  double dis = 0.0;
  double d;
  uint   i;
  
  for (i = 0; i < s; ++i)
    {
      d = p1[i] - p2[i];
      dis += d * d;
    }
  return sqrt(dis);
}

inline double EuclidianTorusDistance
(const uint s, const double *p1, const double *p2)
{
  double dis = 0.0;
  double d;
  uint   i;
  
  for (i = 0; i < s; ++i)
    {
      if ((d = p1[i] - p2[i]) > 0.5)
	d = 1.0 - d;
      dis += d * d;
    }
  return sqrt(dis);
}

/////////////////////////////////////////////////////////////////////////////

double MinimumDistance(const uint s, const uint n, double *p)
{
  double dis = sqrt((double)s);
  double d;
  uint   i, j;
  
  for (i = 0; i < s*n-s; i += s)
    for (j = i+s; j < s*n; j += s)
      if ((d = EuclidianDistance(s, p+i, p+j)) < dis)
	dis = d;
  return dis;
}

/////////////////////////////////////////////////////////////////////////////

double MinimumDistanceTorus(const uint s, const uint n, double *p)
{
  double dis = sqrt((double)s);
  double d;
  uint   i, j;
  
  for (i = 0; i < s*n-s; i += s)
    for (j = i+s; j < s*n; j += s)
      if ((d = EuclidianTorusDistance(s, p+i, p+j)) < dis)
	dis = d;
  return dis;
}

/////////////////////////////////////////////////////////////////////////////

#define _epsilon 0.000001

inline void sub(double *res, const double *p1, const double *p2)
{
  res[0] = p1[0] - p2[0];
  res[1] = p1[1] - p2[1];
  res[2] = p1[2] - p2[2];
}

inline void scale(double *p, const double alpha)
{
  p[0] *= alpha;
  p[1] *= alpha;
  p[2] *= alpha;
}

inline double dot(const double *p1, const double *p2)
{
  return p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
}

inline void cross(double *res, const double *p1, const double *p2)
{
  res[0] = p1[1]*p2[2] - p1[2]*p2[1];
  res[1] = p1[2]*p2[0] - p1[0]*p2[2];
  res[2] = p1[0]*p2[1] - p1[1]*p2[0];
}

inline bool normalize(double *p)
{
  double length;

  if ((length = sqrt(dot(p, p))) < _epsilon)
    return true;
  else
    {
      scale(p, 1.0 / length);
      return false;
    }
}

/////////////////////////////////////////////////////////////////////////////

double CapDiscrepancy(const uint n, double *p)
{
  double dis = 1.0;
  double cup, cdown, cplane;
  double ij[3], ik[3], h[3];
  double c[3], no[3], costheta;
  double d;
  double cap, dish;
  uint   i, j, k, l;
  
  for (i = 0; i < 3*n-3; i += 3)
    for (j = i+3; j < 3*n; j += 3)
      {
	sub(ij, p+j, p+i);
	scale(ij, 0.5);
	sub(c, p+j, ij);
	memcpy(no, c, 3*sizeof(double));
	if (normalize(no))
	  continue;
	costheta = dot(p+i, no);
	cup = cdown = 0;
	for (l = 0; l < i; l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      ++cup;
	    else if (d < -_epsilon)
	      ++cdown;
	  }
	for (l += 3; l < j; l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      ++cup;
	    else if (d < -_epsilon)
	      ++cdown;
	  }
	for (l += 3; l < 3*n; l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      ++cup;
	    else if (d < -_epsilon)
	      ++cdown;
	  }
	cplane = n - cup - cdown;
	cap = (double)n * (0.5 - 0.5*costheta);
	if ((dish = fabs(cap - cup)) > dis)
	  dis = dish;
	if ((dish = fabs(cap - cup - cplane)) > dis)
	  dis = dish;
      }
  if (n > 2)
    for (i = 0; i < 3*n-6; i += 3)
      for (j = i+3; j < 3*n-3; j += 3)
	for (k = j+3; k < 3*n; k += 3)
	  {
	    sub(ij, p+j, p+i);
	    sub(ik, p+k, p+i);
	    cross(no, ij, ik);
	    if (normalize(no))
	      continue;
	    memcpy(c, no, 3*sizeof(double));
	    costheta = dot(p+i, no);
	    scale(c, costheta);
	    cup = cdown = 0;
	    for (l = 0; l < i; l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  ++cup;
		else if (d < -_epsilon)
		  ++cdown;
	      }
	    for (l += 3; l < j; l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  ++cup;
		else if (d < -_epsilon)
		  ++cdown;
	      }
	    for (l += 3; l < k; l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  ++cup;
		else if (d < -_epsilon)
		  ++cdown;
	      }
	    for (l += 3; l < 3*n; l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  ++cup;
		else if (d < -_epsilon)
		  ++cdown;
	      }
	    cplane = n - cup - cdown;
	    cap = (double)n * (0.5 - 0.5*costheta);
	    if ((dish = fabs(cap - cup)) > dis)
	      dis = dish;
	    if ((dish = fabs(cap - cup - cplane)) > dis)
	      dis = dish;
	  }
  return dis / (double)n;
}

/////////////////////////////////////////////////////////////////////////////

double MaximumAngle(const uint n, double *p)
{
  double ang = 0.0;
  uint   pin;
  double ij[3], ik[3], h[3];
  double c[3], no[3];
  double d;
  uint   i, j, k, l;
  
  if (n == 1)
    return M_PI;
  for (i = 0; i < 3*n-3; i += 3)
    for (j = i+3; j < 3*n; j += 3)
      {
	sub(ij, p+j, p+i);
	scale(ij, 0.5);
	sub(c, p+j, ij);
	memcpy(no, c, 3*sizeof(double));
	if (normalize(no))
	  continue;
	pin = 3;
	for (l = 0; pin && (l < i); l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      pin &= 2;
	    else if (d < -_epsilon)
	      pin &= 1;
	  }
	for (l += 3; pin && (l < j); l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      pin &= 2;
	    else if (d < -_epsilon)
	      pin &= 1;
	  }
	for (l += 3; pin && (l < 3*n); l += 3)
	  {
	    sub(h, p+l, p+i);
	    d = dot(no, h);
	    if (d > _epsilon)
	      pin &= 2;
	    else if (d < -_epsilon)
	      pin &= 1;
	  }
	if (pin & 1)
	  if ((d = acos(dot(no, p+i))) > ang)
	    ang = d;
	if (pin & 2)
	  if ((d = M_PI-acos(dot(no, p+i))) > ang)
	    ang = d;
      }
  if (n > 2)
    for (i = 0; i < 3*n-6; i += 3)
      for (j = i+3; j < 3*n-3; j += 3)
	for (k = j+3; k < 3*n; k += 3)
	  {
	    sub(ij, p+j, p+i);
	    sub(ik, p+k, p+i);
	    cross(no, ij, ik);
	    if (normalize(no))
	      continue;
	    memcpy(c, no, 3*sizeof(double));
	    scale(c, dot(p+i, no));
	    pin = 3;
	    for (l = 0; pin && (l < i); l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  pin &= 2;
		else if (d < -_epsilon)
		  pin &= 1;
	      }
	    for (l += 3; pin && (l < j); l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  pin &= 2;
		else if (d < -_epsilon)
		  pin &= 1;
	      }
	    for (l += 3; pin && (l < k); l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  pin &= 2;
		else if (d < -_epsilon)
		  pin &= 1;
	      }
	    for (l += 3; pin && (l < 3*n); l += 3)
	      {
		sub(h, p+l, p+i);
		d = dot(no, h);
		if (d > _epsilon)
		  pin &= 2;
		else if (d < -_epsilon)
		  pin &= 1;
	      }
	    if (pin & 1)
	      if ((d = acos(dot(no, p+i))) > ang)
		ang = d;
	    if (pin & 2)
	      if ((d = M_PI-acos(dot(no, p+i))) > ang)
		ang = d;
	  }
  return ang;
}

/////////////////////////////////////////////////////////////////////////////

double MinimumAngle(const uint n, double *p)
{
  double ang = -1.0;
  double d;
  uint   i, j;
  
  for (i = 0; i < 3*n-3; i += 3)
    for (j = i+3; j < 3*n; j += 3)
      if ((d = dot(p+i, p+j)) > ang)
	ang = d;
  return acos(ang);
}

/////////////////////////////////////////////////////////////////////////////
