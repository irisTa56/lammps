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

#ifndef LMP_RANMARSZIGG_H
#define LMP_RANMARSZIGG_H

#include "pointers.h"

namespace LAMMPS_NS {

class RanMarsZigg : protected Pointers {

  // number of segments
  const int ZIGG_N = 128;

  // position of right-most step
  const double ZIGG_R = 3.442619855896652;

  // area (all segments have the same area)
  const double ZIGG_S = 9.91256303533647e-03;

 public:
  RanMarsZigg(class LAMMPS *, int);
  ~RanMarsZigg();
  double uniform();
  double gaussian();

 private:

  // for Marsaglia
  double *u;
  int i97,j97;
  double c,cd,cm;

  // for Ziggurat
  double *x_table,*y_table;
  double invr;

};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Marsaglia random # generator

The initial seed for this random number generator must be a positive
integer less than or equal to 900 million.

*/
