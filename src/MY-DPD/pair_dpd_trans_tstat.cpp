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

#include <cmath>
#include "pair_dpd_trans_tstat.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "neigh_list.h"
#include "comm.h"
#include "random_ziggurat.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDTransTstat::PairDPDTransTstat(LAMMPS *lmp) : PairDPDTrans(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

void PairDPDTransTstat::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,rinv,ex,ey,ez;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz,randx,randy,randz;
  double dotv,dotr,wd,wr,fd_trans,fr_trans,fpair,fx,fy,fz,factor_dpd;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // adjust sigma if target T is changing

  if (t_start != t_stop) {
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    temperature = t_start + delta * (t_stop-t_start);
    double boltz = force->boltz;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++) {
        sigma[i][j] = sigma[j][i] = sqrt(2.0*boltz*temperature*gamma[i][j]);
        sigma_trans[i][j] = sigma_trans[j][i] = sqrt(2.0*boltz*temperature*gamma_trans[i][j]);
      }
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        rinv = 1.0/r;
        ex = delx * rinv;
        ey = dely * rinv;
        ez = delz * rinv;

        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dotv = ex*delvx + ey*delvy + ez*delvz;

        wr = 1.0 - r*invcut[itype][jtype];
        wd = wr * wr;

        randx = random->gaussian();
        randy = random->gaussian();
        randz = random->gaussian();
        dotr = ex*randx + ey*randy + ez*randz;

        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        /*=== standard part ===*/

        fpair = 0.0;
        fpair -= gamma[itype][jtype]*wd*dotv;
        fpair += sigma[itype][jtype]*wr*dotr*dtinvsqrt;

        fx = ex*fpair;
        fy = ey*fpair;
        fz = ez*fpair;

        /*=== transverse part ===*/

        fd_trans = gamma_trans[itype][jtype]*wd;
        fx -= fd_trans*(delvx-dotv*ex);
        fy -= fd_trans*(delvy-dotv*ey);
        fz -= fd_trans*(delvz-dotv*ez);

        fr_trans = sigma_trans[itype][jtype]*wr*dtinvsqrt;
        fx += fr_trans*(randx-dotr*ex);
        fy += fr_trans*(randy-dotr*ey);
        fz += fr_trans*(randz-dotr*ez);

        fx *= factor_dpd;
        fy *= factor_dpd;
        fz *= factor_dpd;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDTransTstat::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  t_start = force->numeric(FLERR,arg[0]);
  t_stop = force->numeric(FLERR,arg[1]);
  cut_global = force->numeric(FLERR,arg[2]);
  seed = force->inumeric(FLERR,arg[3]);

  temperature = t_start;

  // initialize RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanZiggurat(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
          invcut[i][j] = 1.0/cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDTransTstat::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a0_one = 0.0;
  double gamma_one = force->numeric(FLERR,arg[2]);
  double gamma_trans_one = force->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      gamma_trans[i][j] = gamma_trans_one;
      cut[i][j] = cut_one;
      invcut[i][j] = 1.0/cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDTransTstat::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&gamma_trans[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&invcut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDTransTstat::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&gamma[i][j],sizeof(double),1,fp);
          fread(&gamma_trans[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          fread(&invcut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma_trans[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&invcut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDTransTstat::write_restart_settings(FILE *fp)
{
  fwrite(&t_start,sizeof(double),1,fp);
  fwrite(&t_stop,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDTransTstat::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&t_start,sizeof(double),1,fp);
    fread(&t_stop,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&t_start,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_stop,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  temperature = t_start;

  // initialize RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanZiggurat(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDPDTransTstat::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,gamma[i][i],gamma_trans[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDPDTransTstat::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",
              i,j,gamma[i][j],gamma_trans[i][j],cut[i][j]);
}
