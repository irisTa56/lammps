/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_RANZIGGURAT_H
#define LMP_RANZIGGURAT_H

#include "pointers.h"

namespace LAMMPS_NS {

class RanZiggurat : protected Pointers {

  // this program works only if # of segments is 128)
  const int ZIGG_N = 128;

  // position of right-most step
  const double ZIGG_R = 3.442619855896652;

  // area (all segments have the same area)
  const double ZIGG_S = 9.91256303533647e-03;

 public:
  RanZiggurat(class LAMMPS *, int);
  ~RanZiggurat();
  double uniform();
  uint32_t xorshift();
  double gaussian();

 private:

  // for Marsaglia
  double *u;
  int i97,j97;
  double c,cd,cm;

  // for Xorshift
  uint32_t state;

  // for Ziggurat
  uint32_t *ktab;
  double *wtab,*ytab;
  double invR;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Marsaglia random # generator

The initial seed for this random number generator must be a positive
integer less than or equal to 900 million.

*/
