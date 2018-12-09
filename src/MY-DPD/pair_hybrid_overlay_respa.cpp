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
#include <regex>
#include <string>
#include "pair_hybrid_overlay_respa.h"
#include "atom.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

/* NOTE:
  compute/tally is added to a substyle only if the compute id contains
  the substyle name; for dpd/trans/tstat, for example, id should be
  something like stress_dpd_trans_tstat.
*/
void PairHybridOverlayRespa::add_tally_callback(Compute *ptr)
{
  for (int m = 0; m < nstyles; m++)
  {
    const char* style_name = std::regex_replace(
      std::string(keywords[m]), std::regex("/"), "_").c_str();

    if (compute_tally[m] && strstr(ptr->id, style_name))
    {
      styles[m]->add_tally_callback(ptr);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairHybridOverlayRespa::del_tally_callback(Compute *ptr)
{
  for (int m = 0; m < nstyles; m++)
  {
    const char* style_name = std::regex_replace(
      std::string(keywords[m]), std::regex("/"), "_").c_str();

    if (compute_tally[m] && strstr(ptr->id, style_name))
    {
      styles[m]->del_tally_callback(ptr);
    }
  }
}
