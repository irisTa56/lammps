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

// Marsaglia random number generator
// see RANMAR in F James, Comp Phys Comm, 60, 329 (1990)

// Ziggurat algorithm for generating Gaussian random number
// see https://www.seehuhn.de/pages/ziggurat
// (By Jochen Voss, last updated 2014-06-11)

#include <cmath>
#include "random_mars_zigg.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RanMarsZigg::RanMarsZigg(LAMMPS *lmp, int seed) : Pointers(lmp),
u(NULL), xtab(NULL), ytab(NULL)
{
  // for Marsaglia

  int ij,kl,i,j,k,l,ii,jj,m;
  double s,t;

  if (seed <= 0 || seed > 900000000)
  {
    error->one(FLERR,"Invalid seed for Marsaglia random # generator");
  }

  u = new double[97+1];

  ij = (seed-1)/30082;
  kl = (seed-1) - 30082*ij;
  i = (ij/177) % 177 + 2;
  j = ij %177 + 2;
  k = (kl/169) % 178 + 1;
  l = kl % 169;

  for (ii = 1; ii <= 97; ii++)
  {
    s = 0.0;
    t = 0.5;

    for (jj = 1; jj <= 24; jj++)
    {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s = s + t;
      t = 0.5*t;
    }

    u[ii] = s;
  }

  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;

  uniform();

  // for Ziggurat

  xtab = new double[ZIGG_N+1];
  ytab = new double[ZIGG_N];

  xtab[ZIGG_N-1] = ZIGG_R;
  ytab[ZIGG_N-1] = exp(-0.5*ZIGG_R*ZIGG_R);

  for (int i = ZIGG_N-1; i > 1; --i)
  {
    xtab[i-1] = sqrt(-2.0*log(
      ZIGG_S/xtab[i] + exp(-0.5*xtab[i]*xtab[i])));
    ytab[i-1] = exp(-0.5*xtab[i-1]*xtab[i-1]);
  }

  xtab[0] = 0.0;
  ytab[0] = 1.0;

  xtab[ZIGG_N] = ZIGG_S*ytab[ZIGG_N-1];

  invr = 1.0 / ZIGG_R;
}

/* ---------------------------------------------------------------------- */

RanMarsZigg::~RanMarsZigg()
{
  delete [] u;
  delete [] xtab;
  delete [] ytab;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanMarsZigg::uniform()
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
   gaussian RN
------------------------------------------------------------------------- */

double RanMarsZigg::gaussian()
{
  int i;
  double uni,x,y;

  while (true)
  {
    i = ZIGG_N * uniform();

    x = xtab[i+1] * uniform();

    if (x < xtab[i]) break;

    if (i == ZIGG_N-1)
    {
      x = ZIGG_R - log(1.0-uniform()) * invr;

      if (-2.0*log(1.0-uniform()) > (x-ZIGG_R)*(x-ZIGG_R)) break;
    }
    else
    {
      uni = uniform();
      y = (1-uni)*ytab[i+1] + uni*ytab[i];

      if ((1-uni)*ytab[i+1] + uni*ytab[i] < exp(-0.5*x*x)) break;
    }
  }

  return uniform() < 0.5 ? x : -x;
}
