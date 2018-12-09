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
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "comm_brick_dpd.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "math_extra.h"
#include "error.h"
#include "memory.h"
#include "pair_hybrid_overlay_respa.h"

using namespace LAMMPS_NS;

#define BUFMIN 1000
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the CommBrick initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for CommBrick is run and thus creating a shallow copy of "oldcomm".
//           The call to CommBrick::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommBrickDPD::CommBrickDPD(LAMMPS *, CommBrick *oldcomm) :
  CommBrick(*oldcomm), size_forward_recv_dpd(NULL),
  size_reverse_send_dpd(NULL), size_reverse_recv_dpd(NULL),
  sendnum_dpd(NULL), recvnum_dpd(NULL)
{
  style = oldcomm->style;
  layout = oldcomm->layout;
  CommBrick::copy_arrays(oldcomm);
  init_buffers();

  comm_only_dpd = 0;

  ghost_velocity = 1;
}

/* ---------------------------------------------------------------------- */

void CommBrickDPD::init()
{
  if (mode != Comm::SINGLE) {
    error->all(FLERR,"CommBrickDPD supports only single-mode");
  }

  if (triclinic) {
    error->all(FLERR,"CommBrickDPD does not support triclinic");
  }

  CommBrick::init();

  PairHybridOverlayRespa *hybrid
    = dynamic_cast<PairHybridOverlayRespa *>(force->pair);

  double maxcut = 0.0;
  double maxcut_dpd = 0.0;

  int n = atom->ntypes;

  for (int i = 0; i < hybrid->nstyles; ++i)
  {
    bool is_dpd = strstr((hybrid->keywords)[i], "dpd");

    Pair *style_i = (hybrid->styles)[i];

    for (int ii = 1; ii <= n; ++ii)
    {
      for (int jj = 1; jj <= n; ++jj)
      {
        double cut = sqrt(style_i->cutsq[ii][jj]);

        if (is_dpd)
        {
          if (maxcut_dpd < cut) maxcut_dpd = cut;
        }
        else
        {
          if (maxcut < cut) maxcut = cut;
        }
      }
    }
  }

  cutdiff = maxcut - maxcut_dpd;

  if (comm->me == 0)
  {
    if (screen)
    {
      fprintf(screen,"Cutoff difference = %f\n",cutdiff);
    }
    if (logfile)
    {
      fprintf(logfile,"Cutoff difference = %f\n",cutdiff);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CommBrickDPD::presetup(double cut)
{
  double *prd,*sublo,*subhi;

  prd = domain->prd;
  sublo = domain->sublo;
  subhi = domain->subhi;

  cutghost[0] = cutghost[1] = cutghost[2] = cut;

  int *periodicity = domain->periodicity;
  int left,right;

  if (layout == Comm::LAYOUT_UNIFORM)
  {
    maxneed[0] = static_cast<int> (cutghost[0] * procgrid[0] / prd[0]) + 1;
    maxneed[1] = static_cast<int> (cutghost[1] * procgrid[1] / prd[1]) + 1;
    maxneed[2] = static_cast<int> (cutghost[2] * procgrid[2] / prd[2]) + 1;
    if (domain->dimension == 2) maxneed[2] = 0;
    if (!periodicity[0]) maxneed[0] = MIN(maxneed[0],procgrid[0]-1);
    if (!periodicity[1]) maxneed[1] = MIN(maxneed[1],procgrid[1]-1);
    if (!periodicity[2]) maxneed[2] = MIN(maxneed[2],procgrid[2]-1);

    if (!periodicity[0]) {
      recvneed[0][0] = MIN(maxneed[0],myloc[0]);
      recvneed[0][1] = MIN(maxneed[0],procgrid[0]-myloc[0]-1);
      left = myloc[0] - 1;
      if (left < 0) left = procgrid[0] - 1;
      sendneed[0][0] = MIN(maxneed[0],procgrid[0]-left-1);
      right = myloc[0] + 1;
      if (right == procgrid[0]) right = 0;
      sendneed[0][1] = MIN(maxneed[0],right);
    } else recvneed[0][0] = recvneed[0][1] =
             sendneed[0][0] = sendneed[0][1] = maxneed[0];

    if (!periodicity[1]) {
      recvneed[1][0] = MIN(maxneed[1],myloc[1]);
      recvneed[1][1] = MIN(maxneed[1],procgrid[1]-myloc[1]-1);
      left = myloc[1] - 1;
      if (left < 0) left = procgrid[1] - 1;
      sendneed[1][0] = MIN(maxneed[1],procgrid[1]-left-1);
      right = myloc[1] + 1;
      if (right == procgrid[1]) right = 0;
      sendneed[1][1] = MIN(maxneed[1],right);
    } else recvneed[1][0] = recvneed[1][1] =
             sendneed[1][0] = sendneed[1][1] = maxneed[1];

    if (!periodicity[2]) {
      recvneed[2][0] = MIN(maxneed[2],myloc[2]);
      recvneed[2][1] = MIN(maxneed[2],procgrid[2]-myloc[2]-1);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = MIN(maxneed[2],procgrid[2]-left-1);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = MIN(maxneed[2],right);
    } else recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = maxneed[2];
  }
  else
  {
    recvneed[0][0] = updown(0,0,myloc[0],prd[0],periodicity[0],xsplit);
    recvneed[0][1] = updown(0,1,myloc[0],prd[0],periodicity[0],xsplit);
    left = myloc[0] - 1;
    if (left < 0) left = procgrid[0] - 1;
    sendneed[0][0] = updown(0,1,left,prd[0],periodicity[0],xsplit);
    right = myloc[0] + 1;
    if (right == procgrid[0]) right = 0;
    sendneed[0][1] = updown(0,0,right,prd[0],periodicity[0],xsplit);

    recvneed[1][0] = updown(1,0,myloc[1],prd[1],periodicity[1],ysplit);
    recvneed[1][1] = updown(1,1,myloc[1],prd[1],periodicity[1],ysplit);
    left = myloc[1] - 1;
    if (left < 0) left = procgrid[1] - 1;
    sendneed[1][0] = updown(1,1,left,prd[1],periodicity[1],ysplit);
    right = myloc[1] + 1;
    if (right == procgrid[1]) right = 0;
    sendneed[1][1] = updown(1,0,right,prd[1],periodicity[1],ysplit);

    if (domain->dimension == 3) {
      recvneed[2][0] = updown(2,0,myloc[2],prd[2],periodicity[2],zsplit);
      recvneed[2][1] = updown(2,1,myloc[2],prd[2],periodicity[2],zsplit);
      left = myloc[2] - 1;
      if (left < 0) left = procgrid[2] - 1;
      sendneed[2][0] = updown(2,1,left,prd[2],periodicity[2],zsplit);
      right = myloc[2] + 1;
      if (right == procgrid[2]) right = 0;
      sendneed[2][1] = updown(2,0,right,prd[2],periodicity[2],zsplit);
    } else recvneed[2][0] = recvneed[2][1] =
             sendneed[2][0] = sendneed[2][1] = 0;

    int all[6];
    MPI_Allreduce(&recvneed[0][0],all,6,MPI_INT,MPI_MAX,world);
    maxneed[0] = MAX(all[0],all[1]);
    maxneed[1] = MAX(all[2],all[3]);
    maxneed[2] = MAX(all[4],all[5]);
  }
}

/* ---------------------------------------------------------------------- */

void CommBrickDPD::setup()
{
  presetup(MAX(neighbor->cutneighmax-cutdiff,cutghostuser));

  maxneed_dpd[0] = maxneed[0];
  maxneed_dpd[1] = maxneed[1];
  maxneed_dpd[2] = maxneed[2];

  cutghost_dpd[0] = cutghost[0];
  cutghost_dpd[1] = cutghost[1];
  cutghost_dpd[2] = cutghost[2];

  presetup(MAX(neighbor->cutneighmax,cutghostuser));

  if (comm->me == 0)
  {
    if (screen)
    {
      fprintf(
        screen,"Max # of swaps: 2 * (%d * %d * %d)\n",
        maxneed[0],maxneed[1],maxneed[2]);
      fprintf(
        screen,"Max # of swaps for DPD: 2 * (%d * %d * %d)\n",
        maxneed_dpd[0],maxneed_dpd[1],maxneed_dpd[2]);
    }
    if (logfile)
    {
      fprintf(
        logfile,"Max # of swaps: 2 * (%d * %d * %d)\n",
        maxneed[0],maxneed[1],maxneed[2]);
      fprintf(
        logfile,"Max # of swaps for DPD: 2 * (%d * %d * %d)\n",
        maxneed_dpd[0],maxneed_dpd[1],maxneed_dpd[2]);
    }
  }

  double *sublo,*subhi;

  sublo = domain->sublo;
  subhi = domain->subhi;

  // allocate comm memory

  nswap = 2 * (maxneed[0]+maxneed[1]+maxneed[2]);
  if (nswap > maxswap) grow_swap(nswap);

  int iswap = 0;
  for (int dim = 0; dim < 3; dim++) {
    for (int ineed = 0; ineed < 2*maxneed[dim]; ineed++) {
      pbc_flag[iswap] = 0;
      pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] =
        pbc[iswap][3] = pbc[iswap][4] = pbc[iswap][5] = 0;

      if (ineed % 2 == 0) {
        sendproc[iswap] = procneigh[dim][0];
        recvproc[iswap] = procneigh[dim][1];

        if (ineed < 2) slablo[iswap] = -BIG;
        else slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);

        slabhi[iswap] = sublo[dim] + cutghost[dim];

        if (myloc[dim] == 0) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = 1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = 1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = 1;
          }
        }

      } else {
        sendproc[iswap] = procneigh[dim][1];
        recvproc[iswap] = procneigh[dim][0];

        slablo[iswap] = subhi[dim] - cutghost[dim];

        if (ineed < 2) slabhi[iswap] = BIG;
        else slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);

        if (myloc[dim] == procgrid[dim]-1) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = -1;
          if (triclinic) {
            if (dim == 1) pbc[iswap][5] = -1;
            else if (dim == 2) pbc[iswap][4] = pbc[iswap][3] = -1;
          }
        }
      }

      iswap++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void CommBrickDPD::forward_comm(int only_dpd)
{
  comm_only_dpd = only_dpd;

  MPI_Request request;
  AtomVec *avec = atom->avec;

  int iswap = 0;

  for (int dim = 0; dim != 3; ++dim)
  {
    int twoneed = comm_only_dpd ? 2*maxneed_dpd[dim] : 2*maxneed[dim];

    for (int ineed = 0; ineed != twoneed; ++ineed)
    {
      int size_forward_recv_ = comm_only_dpd ?
        size_forward_recv_dpd[iswap] : size_forward_recv[iswap];
      int sendnum_ = comm_only_dpd ?
        sendnum_dpd[iswap] : sendnum[iswap];
      int recvnum_ = comm_only_dpd ?
        recvnum_dpd[iswap] : recvnum[iswap];

      if (sendproc[iswap] != me)
      {
        if (size_forward_recv_)
        {
          MPI_Irecv(buf_recv,size_forward_recv_,
            MPI_DOUBLE,recvproc[iswap],0,world,&request);
        }

        int n = avec->pack_comm_vel(
          sendnum_,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);

        if (n)
        {
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        }

        if (size_forward_recv_)
        {
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }

        avec->unpack_comm_vel(recvnum_,firstrecv[iswap],buf_recv);
      }
      else
      {
        avec->pack_comm_vel(
          sendnum_,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);
        avec->unpack_comm_vel(recvnum_,firstrecv[iswap],buf_send);
      }

      ++iswap;
    }

    iswap += 2*maxneed[dim] - twoneed;
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommBrickDPD::reverse_comm()
{
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  int iswap = nswap - 1;

  for (int dim = 2; dim != -1; --dim)
  {
    int twoneed = comm_only_dpd ? 2*maxneed_dpd[dim] : 2*maxneed[dim];

    iswap -= 2*maxneed[dim] - twoneed;

    for (int ineed = 0; ineed != twoneed; ++ineed)
    {
      int size_reverse_recv_ = comm_only_dpd ?
        size_reverse_recv_dpd[iswap] : size_reverse_recv[iswap];
      int size_reverse_send_ = comm_only_dpd ?
        size_reverse_send_dpd[iswap] : size_reverse_send[iswap];
      int sendnum_ = comm_only_dpd ?
        sendnum_dpd[iswap] : sendnum[iswap];
      int recvnum_ = comm_only_dpd ?
        recvnum_dpd[iswap] : recvnum[iswap];

      if (sendproc[iswap] != me)
      {
        if (comm_f_only)
        {
          if (size_reverse_recv_)
          {
            MPI_Irecv(
              buf_recv,size_reverse_recv_,
              MPI_DOUBLE,sendproc[iswap],0,world,&request);
          }

          if (size_reverse_send_)
          {
            buf = f[firstrecv[iswap]];
            MPI_Send(
              buf,size_reverse_send_,MPI_DOUBLE,recvproc[iswap],0,world);
          }

          if (size_reverse_recv_)
          {
            MPI_Wait(&request,MPI_STATUS_IGNORE);
          }
        }
        else
        {
          if (size_reverse_recv_)
          {
            MPI_Irecv(
              buf_recv,size_reverse_recv_,
              MPI_DOUBLE,sendproc[iswap],0,world,&request);
          }

          int n = avec->pack_reverse(recvnum_,firstrecv[iswap],buf_send);

          if (n)
          {
            MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
          }

          if (size_reverse_recv_)
          {
            MPI_Wait(&request,MPI_STATUS_IGNORE);
          }
        }

        avec->unpack_reverse(sendnum_,sendlist[iswap],buf_recv);
      }
      else
      {
        if (comm_f_only)
        {
          if (sendnum_)
          {
            avec->unpack_reverse(
              sendnum_,sendlist[iswap],f[firstrecv[iswap]]);
          }
        }
        else
        {
          avec->pack_reverse(recvnum_,firstrecv[iswap],buf_send);
          avec->unpack_reverse(sendnum_,sendlist[iswap],buf_send);
        }
      }

      --iswap;
    }
  }

  comm_only_dpd = 0;
}

/* ---------------------------------------------------------------------- */

void CommBrickDPD::borders()
{
  double **x;
  double *buf;

  MPI_Request request;
  AtomVec *avec = atom->avec;

  // max size in atoms of single borders send/recv
  smax = rmax = 0;

  std::vector<int> list_dpd;
  std::vector<int> list_pair;

  list_dpd.reserve(BUFMIN);
  list_pair.reserve(BUFMIN);

  int nfirst,nlast,nsend,nrecv,nsend_dpd,nrecv_dpd;

  int iswap = 0;

  for (int dim = 0; dim != 3; ++dim)
  {
    nlast = 0;
    int twoneed = 2*maxneed[dim];

    for (int ineed = 0; ineed != twoneed; ++ineed)
    {
      // whether DPD needs this swap
      bool dpd_needs = ineed < 2*maxneed_dpd[dim];

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in sendlist for use in future timesteps

      x = atom->x;

      double lo = slablo[iswap];
      double hi = slabhi[iswap];
      double lo_dpd = lo + cutdiff;
      double hi_dpd = hi - cutdiff;

      if (ineed % 2 == 0)
      {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      // compute # of atoms to be sent

      nsend = 0;
      nsend_dpd = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      bool sendflag = ineed/2 >= sendneed[dim][ineed%2] ? 0 : 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus only atoms in bordergroup
      // can only limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag)
      {
        list_dpd.clear();
        list_pair.clear();

        if (!bordergroup || ineed >= 2)
        {
          if (dpd_needs)
          {
            if (ineed % 2 == 0)
            {
              for (int i = nfirst; i != nlast; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] <= hi_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }
            }
            else
            {
              for (int i = nfirst; i != nlast; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] >= lo_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }
            }
          }
          else
          {
            for (int i = nfirst; i != nlast; ++i)
            {
              if (x[i][dim] >= lo && x[i][dim] <= hi)
              {
                list_pair.push_back(i);
              }
            }
          }
        }
        else
        {
          int ngroup = atom->nfirst;

          if (dpd_needs)
          {
            if (ineed % 2 == 0)
            {
              for (int i = 0; i != ngroup; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] <= hi_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }

              for (int i = atom->nlocal; i != nlast; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] <= hi_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }
            }
            else
            {
              for (int i = 0; i != ngroup; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] >= lo_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }

              for (int i = atom->nlocal; i != nlast; ++i)
              {
                if (x[i][dim] >= lo && x[i][dim] <= hi)
                {
                  if (x[i][dim] >= lo_dpd) list_dpd.push_back(i);
                  else list_pair.push_back(i);
                }
              }
            }
          }
          else
          {
            for (int i = 0; i != ngroup; ++i)
            {
              if (x[i][dim] >= lo && x[i][dim] <= hi)
              {
                list_pair.push_back(i);
              }
            }

            for (int i = atom->nlocal; i != nlast; ++i)
            {
              if (x[i][dim] >= lo && x[i][dim] <= hi)
              {
                list_pair.push_back(i);
              }
            }
          }
        }

        nsend_dpd = list_dpd.size();
        nsend = nsend_dpd + list_pair.size();

        if (nsend > maxsendlist[iswap]) grow_list(iswap,nsend);

        for (int i = 0; i != nsend_dpd; ++i)
        {
          sendlist[iswap][i] = list_dpd[i];
        }

        int nsend_pair = nsend - nsend_dpd;

        for (int i = 0; i != nsend_pair; ++i)
        {
          sendlist[iswap][i+nsend_dpd] = list_pair[i];
        }
      }

      // pack up list of border atoms

      if (nsend*size_border > maxsend) grow_send(nsend*size_border,0);

      int n = avec->pack_border_vel(
        nsend,sendlist[iswap],buf_send,pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me)
      {
        MPI_Sendrecv(
          &nsend,1,MPI_INT,sendproc[iswap],0,
          &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);

        if (dpd_needs)
        {
          MPI_Sendrecv(
            &nsend_dpd,1,MPI_INT,sendproc[iswap],0,
            &nrecv_dpd,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
        }
        else
        {
          nrecv_dpd = 0;
        }

        if (nrecv*size_border > maxrecv) grow_recv(nrecv*size_border);

        if (nrecv)
        {
          MPI_Irecv(
            buf_recv,nrecv*size_border,
            MPI_DOUBLE,recvproc[iswap],0,world,&request);
        }

        if (n)
        {
          MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        }

        if (nrecv)
        {
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }

        buf = buf_recv;
      }
      else
      {
        nrecv = nsend;
        nrecv_dpd = nsend_dpd;
        buf = buf_send;
      }

      // unpack buffer

      avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      sendnum_dpd[iswap] = nsend_dpd;
      recvnum[iswap] = nrecv;
      recvnum_dpd[iswap] = nrecv_dpd;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_forward_recv_dpd[iswap] = nrecv_dpd*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_send_dpd[iswap] = nrecv_dpd*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      size_reverse_recv_dpd[iswap] = nsend_dpd*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;

      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void CommBrickDPD::allocate_swap(int n)
{
  CommBrick::allocate_swap(n);

  memory->create(sendnum_dpd,n,"comm:sendnum");
  memory->create(recvnum_dpd,n,"comm:recvnum");
  memory->create(size_forward_recv_dpd,n,"comm:size");
  memory->create(size_reverse_send_dpd,n,"comm:size");
  memory->create(size_reverse_recv_dpd,n,"comm:size");
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void CommBrickDPD::free_swap()
{
  CommBrick::free_swap();

  memory->destroy(sendnum_dpd);
  memory->destroy(recvnum_dpd);
  memory->destroy(size_forward_recv_dpd);
  memory->destroy(size_reverse_send_dpd);
  memory->destroy(size_reverse_recv_dpd);
}
