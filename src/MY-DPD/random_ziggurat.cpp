/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Ziggurat algorithm for generating Gaussian random number
// see https://www.seehuhn.de/pages/ziggurat, Jochen Voss, 2014-06-11

#include <cmath>
#include "random_ziggurat.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RandomZiggurat::RandomZiggurat(LAMMPS *lmp, int seed) : Pointers(lmp),
u(NULL), ktab(NULL), wtab(NULL), ytab(NULL)
{
  if (seed <= 0 || seed > 900000000)
  {
    error->one(FLERR,"Invalid seed for Marsaglia random # generator");
  }

  u = new double[97+1];

  int ij = (seed-1) / 30082;
  int kl = (seed-1) - 30082*ij;
  int i = (ij/177) % 177 + 2;
  int j = ij % 177 + 2;
  int k = (kl/169) % 178 + 1;
  int l = kl % 169;

  for (int ii = 1; ii < 98; ++ii)
  {
    double s = 0.0;
    double t = 0.5;

    for (int jj = 1; jj < 25; ++jj)
    {
      int m = ((i*j)%179) * k % 179;

      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;

      if ((l*m) % 64 > 31) s += t;

      t *= 0.5;
    }

    u[ii] = s;
  }

  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;

  uniform();

  // for Xorshift

  state = seed;

  // for Ziggurat

  ktab = new uint32_t[ZIGG_N];
  wtab = new double[ZIGG_N];
  ytab = new double[ZIGG_N];

  double xtab[ZIGG_N+1];

  xtab[ZIGG_N-1] = ZIGG_R;
  ytab[ZIGG_N-1] = exp(-0.5*ZIGG_R*ZIGG_R);

  for (int ii = ZIGG_N-1; ii > 1; --ii)
  {
    xtab[ii-1] = sqrt(-2.0*log(
      ZIGG_S/xtab[ii] + exp(-0.5*xtab[ii]*xtab[ii])));
    ytab[ii-1] = exp(-0.5*xtab[ii-1]*xtab[ii-1]);
  }

  xtab[0] = 0.0;
  ytab[0] = 1.0;

  xtab[ZIGG_N] = ZIGG_S / ytab[ZIGG_N-1];

  uint32_t pow2_24 = pow(2,24);

  for (int ii = 0; ii < ZIGG_N; ++ii)
  {
    ktab[ii] = pow2_24 * xtab[ii]/xtab[ii+1];
    wtab[ii] = xtab[ii+1] / double(pow2_24);
  }

  invR = 1.0 / ZIGG_R;
}

/* ---------------------------------------------------------------------- */

RandomZiggurat::~RandomZiggurat()
{
  delete [] u;
  delete [] ktab;
  delete [] wtab;
  delete [] ytab;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RandomZiggurat::uniform()
{
  double uni = u[i97] - u[j97];

  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97 == 0) i97 = 97;
  j97--;
  if (j97 == 0) j97 = 97;
  c -= cd;
  if (c < 0.0) c += cm;
  uni -= c;
  if (uni < 0.0) uni += 1.0;
  return uni;
}

/* ----------------------------------------------------------------------
   32-bit RN
------------------------------------------------------------------------- */

uint32_t RandomZiggurat::xorshift()
{
  uint32_t x = state;

  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;

  state = x;

  return x;
}

/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

double RandomZiggurat::gaussian()
{
  uint32_t rand32,i,sign,ix;
  double uni,x;

  while (true)
  {
    rand32 = xorshift();
    i = rand32 & 0x0000007F;  // 7 bit to for segment (0 ... 127 as int)
    sign = rand32 & 0x00000080;  // 1 bit for sign (0 or 128 as int)
    ix = rand32 >> 8;  // 24 bit for x-value

    x = ix*wtab[i];

    if (ix < ktab[i])  break;

    if (i == ZIGG_N-1)
    {
      x = ZIGG_R - log(1.0-uniform()) * invR;
      if (-2.0*log(1.0-uniform()) > (x-ZIGG_R)*(x-ZIGG_R)) break;
    }
    else
    {
      uni = uniform();
      if ((1-uni)*ytab[i+1] + uni*ytab[i] < exp(-0.5*x*x)) break;
    }
  }

  return sign ? x : -x;
}
