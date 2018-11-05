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

#ifndef LMP_COMM_BRICK_DPD_H
#define LMP_COMM_BRICK_DPD_H

#include "comm_brick.h"

namespace LAMMPS_NS {

class CommBrickDPD : public CommBrick {
 public:
  CommBrickDPD(class LAMMPS *, class CommBrick *);

  virtual void init() override;
  virtual void setup() override;
  virtual void forward_comm(int only_dpd = 0) override;
  virtual void reverse_comm() override;
  virtual void borders() override;

 protected:
  int comm_only_dpd;
  double cutdiff;
  int maxneed_dpd[3];
  double cutghost_dpd[3];
  int *size_forward_recv_dpd;
  int *size_reverse_send_dpd,*size_reverse_recv_dpd;
  int *sendnum_dpd,*recvnum_dpd;
  double *slablo_dpd,*slabhi_dpd;

  virtual void allocate_swap(int) override;
  virtual void free_swap() override;

  void presetup(double);
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot change to comm_style brick from tiled layout

Self-explanatory.

*/
