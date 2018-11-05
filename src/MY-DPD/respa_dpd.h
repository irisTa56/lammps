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

#ifdef INTEGRATE_CLASS

IntegrateStyle(respa/dpd,RespaDPD)

#else

#ifndef LMP_RESPA_DPD_H
#define LMP_RESPA_DPD_H

#include "respa.h"

namespace LAMMPS_NS {

class RespaDPD : public Respa {
 public:
  RespaDPD(class LAMMPS *, int, char **);
  virtual ~RespaDPD();
 protected:
  int mlevel;  // MAX(level_pair-1,0)
  int *iloop;
  virtual void recurse(int) override;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Respa levels must be >= 1

Self-explanatory.

E: Cannot set both respa pair and inner/middle/outer

In the rRESPA integrator, you must compute pairwise potentials either
all together (pair), or in pieces (inner/middle/outer).  You can't do
both.

E: Must set both respa inner and outer

Cannot use just the inner or outer option with respa without using the
other.

E: Cannot set respa middle without inner/outer

In the rRESPA integrator, you must define both a inner and outer
setting in order to use a middle setting.

E: Cannot set respa hybrid and any of pair/inner/middle/outer

In the rRESPA integrator, you must compute pairwise potentials either
all together (pair), with different cutoff regions (inner/middle/outer),
or per hybrid sub-style (hybrid).  You cannot mix those.

E: Invalid order of forces within respa levels

For respa, ordering of force computations within respa levels must
obey certain rules.  E.g. bonds cannot be compute less frequently than
angles, pairwise forces cannot be computed less frequently than
kspace, etc.

W: One or more respa levels compute no forces

This is computationally inefficient.

E: Respa inner cutoffs are invalid

The first cutoff must be <= the second cutoff.

E: Respa middle cutoffs are invalid

The first cutoff must be <= the second cutoff.

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

E: Pair style does not support rRESPA inner/middle/outer

You are attempting to use rRESPA options with a pair style that
does not support them.

*/
