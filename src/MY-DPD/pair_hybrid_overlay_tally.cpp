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

#include <cstdlib>
#include <cstring>
#include <cctype>
#include "pair_hybrid_overlay_tally.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairHybridOverlayTally::settings(int narg, char **arg)
{
  PairHybrid::settings(narg,arg);
  compute_tally_ids.resize(nstyles);
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style and its sub-styles
------------------------------------------------------------------------- */

void PairHybridOverlayTally::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal pair_modify command");

  // if 1st keyword is pair, apply other keywords to one sub-style

  if (strcmp(arg[0],"pair") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal pair_modify command");
    int m;
    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[1],keywords[m]) == 0) break;
    if (m == nstyles) error->all(FLERR,"Unknown pair_modify hybrid sub-style");
    int iarg = 2;

    if (multiple[m]) {
      if (narg < 3) error->all(FLERR,"Illegal pair_modify command");
      int multiflag = force->inumeric(FLERR,arg[2]);
      for (m = 0; m < nstyles; m++)
        if (strcmp(arg[1],keywords[m]) == 0 && multiflag == multiple[m]) break;
      if (m == nstyles)
        error->all(FLERR,"Unknown pair_modify hybrid sub-style");
      iarg = 3;
    }

    // if 2nd keyword (after pair) is special:
    // invoke modify_special() for the sub-style

    if (iarg < narg && strcmp(arg[iarg],"special") == 0) {
      if (narg < iarg+5)
        error->all(FLERR,"Illegal pair_modify special command");
      modify_special(m,narg-iarg,&arg[iarg+1]);
      iarg += 5;
    }

    // if 2nd keyword (after pair) is compute/tally:
    // set flag to register USER-TALLY computes accordingly

    if (iarg < narg && strcmp(arg[iarg],"compute/tally") == 0) {
      if (narg < iarg+2)
        error->all(FLERR,"Illegal pair_modify compute/tally command");
      // Only this line differs from PairHybrid::modify_params()
      compute_tally_ids[m].insert(std::string(arg[iarg+1]));
      iarg += 2;
    }

    // apply the remaining keywords to the base pair style itself and the
    // sub-style except for "pair" and "special".
    // the former is important for some keywords like "tail" or "compute"

    if (narg-iarg > 0) {
      Pair::modify_params(narg-iarg,&arg[iarg]);
      styles[m]->modify_params(narg-iarg,&arg[iarg]);
    }

  // apply all keywords to pair hybrid itself and every sub-style

  } else {
    Pair::modify_params(narg,arg);
    for (int m = 0; m < nstyles; m++) styles[m]->modify_params(narg,arg);
  }

  // reset global compute_flag since there may have been changes
  // to any of the substyles
  compute_flag = 0;
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->compute_flag) compute_flag = 1;
}

/* ---------------------------------------------------------------------- */

void PairHybridOverlayTally::add_tally_callback(Compute *ptr)
{
  bool added = false;
  for (int m = 0; m < nstyles; m++) {
    auto tallyid = std::string(ptr->id)+std::string(keywords[m]);
    if (added_tallys.find(tallyid) == added_tallys.end()) {
      auto id_set = compute_tally_ids[m];
      if (id_set.find(std::string(ptr->id)) != id_set.end()) {
        styles[m]->add_tally_callback(ptr);
        added_tallys.insert(tallyid);
        added = true;
        if (comm->me == 0) {
          if (screen)
            fprintf(
              screen,"Compute/tally %s is added to pair %s\n",
              ptr->id,keywords[m]);
          if (logfile)
            fprintf(
              logfile,"Compute/tally %s is added to pair %s\n",
              ptr->id,keywords[m]);
        }
        break;
      }
    } else added = true;
  }
  if (!added)
    error->warning(
      FLERR,("Orphan compute/tally: "+std::string(ptr->id)).c_str());
}

/* ---------------------------------------------------------------------- */

void PairHybridOverlayTally::del_tally_callback(Compute *ptr)
{
  for (int m = 0; m < nstyles; m++) {
    auto tallyid = std::string(ptr->id)+std::string(keywords[m]);
    if (added_tallys.find(tallyid) != added_tallys.end()) {
      styles[m]->del_tally_callback(ptr);
      added_tallys.erase(tallyid);
      break;
    }
  }
}
