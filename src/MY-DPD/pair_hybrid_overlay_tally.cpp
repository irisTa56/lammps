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

    // search only compute/tally

    for (int i = iarg; i < narg; i++) {
      if (strcmp(arg[i],"compute/tally") == 0 && i+1 < narg) {
        compute_tally_ids[m].insert(std::string(arg[i+1]));
        arg[i+1] = (char*)"yes";  // to avoid error in PairHybrid::modify_params()
        break;
      }
    }
  }

  PairHybrid::modify_params(narg,arg);
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
