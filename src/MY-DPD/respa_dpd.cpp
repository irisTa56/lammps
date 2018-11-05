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

/* ----------------------------------------------------------------------
   Contributing authors: Mark Stevens (SNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "respa_dpd.h"
#include "neighbor.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix_respa.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "comm_brick_dpd.h"
#include "pair_hybrid_overlay_respa.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RespaDPD::RespaDPD(LAMMPS *lmp, int narg, char **arg) :
  Respa(lmp, narg, arg), iloop(NULL)
{
  iloop = new int[nlevels];

  /* Respa settings */

  if (strcmp(force->pair_style,"hybrid/overlay/respa"))
  {
    error->all(
      FLERR,"Illegal pair_style: respa/dpd accepts only hybrid/overlay/respa");
  }

  PairHybridOverlayRespa *hybrid
    = dynamic_cast<PairHybridOverlayRespa *>(force->pair);
  nhybrid_styles = hybrid->nstyles;

  hybrid_level = new int[nhybrid_styles];
  hybrid_compute = new int[nhybrid_styles];

  for (int i = 0; i < nhybrid_styles; ++i)
  {
    // DPD forces are computed at the innermost Respa loop
    hybrid_level[i]
      = strstr((hybrid->keywords)[i], "dpd") ? 0 : level_pair;
  }

  mlevel = MAX(level_pair-1,0);
  level_pair = -1;

  /* Comm settings */

  if (comm->style != 0)
  {
    error->all(
      FLERR,"Illegal comm_style: respa/dpd accepts only brick");
  }

  CommBrick *oldcomm = dynamic_cast<CommBrick *>(comm);
  comm = new CommBrickDPD(lmp,oldcomm);
  delete oldcomm;
}

/* ---------------------------------------------------------------------- */

RespaDPD::~RespaDPD()
{
  delete [] iloop;
}

/* ---------------------------------------------------------------------- */

void RespaDPD::recurse(int ilevel)
{
  copy_flevel_f(ilevel);

  for (iloop[ilevel] = 0; iloop[ilevel] < loop[ilevel]; iloop[ilevel]++)
  {
    timer->stamp();
    modify->initial_integrate_respa(vflag,ilevel,iloop[ilevel]);
    if (modify->n_post_integrate_respa)
      modify->post_integrate_respa(ilevel,iloop[ilevel]);
    timer->stamp(Timer::MODIFY);

    // at outermost level, check on rebuilding neighbor list
    // at innermost level, communicate
    // at middle levels, do nothing

    if (ilevel == nlevels-1) {

      if (neighbor->decide()) {
        if (modify->n_pre_exchange) {
          timer->stamp();
          modify->pre_exchange();
          timer->stamp(Timer::MODIFY);
        }
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
          domain->reset_box();
          comm->setup();
          if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 &&
            update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (modify->n_pre_neighbor) {
          modify->pre_neighbor();
          timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_post_neighbor) {
          modify->post_neighbor();
          timer->stamp(Timer::MODIFY);
        }

      } else if (ilevel == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
      }

    }
    else if (ilevel == 0)
    {
      timer->stamp();

      if (iloop[mlevel] == loop[mlevel]-1)
      {
        comm->forward_comm();
      }
      else
      {
        comm->forward_comm(1);
      }

      timer->stamp(Timer::COMM);
    }

    // rRESPA recursion thru all levels
    // this used to be before neigh list build,
    // which prevented per-atom energy/stress being tallied correctly
    // b/c atoms migrated to new procs between short/long force calls
    // now they migrate at very start of rRESPA timestep, before all forces

    if (ilevel) recurse(ilevel-1);

    // force computations
    // important that ordering is same as Verlet
    // so that any order dependencies are the same
    // when potentials are invoked at same level

    force_clear(newton[ilevel]);
    if (modify->n_pre_force_respa) {
      timer->stamp();
      modify->pre_force_respa(vflag,ilevel,iloop[ilevel]);
      timer->stamp(Timer::MODIFY);
    }

    timer->stamp();
    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_pair == ilevel && pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_inner == ilevel && pair_compute_flag) {
      force->pair->compute_inner();
      timer->stamp(Timer::PAIR);
    }
    if (level_middle == ilevel && pair_compute_flag) {
      force->pair->compute_middle();
      timer->stamp(Timer::PAIR);
    }
    if (level_outer == ilevel && pair_compute_flag) {
      force->pair->compute_outer(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_bond == ilevel && force->bond) {
      force->bond->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_angle == ilevel && force->angle) {
      force->angle->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_dihedral == ilevel && force->dihedral) {
      force->dihedral->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_improper == ilevel && force->improper) {
      force->improper->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_kspace == ilevel && kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(Timer::KSPACE);
    }

    if (modify->n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);
    }

    if (newton[ilevel]) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
    }
    timer->stamp();
    if (modify->n_post_force_respa)
      modify->post_force_respa(vflag,ilevel,iloop[ilevel]);
    modify->final_integrate_respa(ilevel,iloop[ilevel]);
    timer->stamp(Timer::MODIFY);
  }

  copy_f_flevel(ilevel);
}
