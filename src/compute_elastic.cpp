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

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "compute_elastic.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeElastic::ComputeElastic(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  v2ptr(NULL), id_press(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute elastic command");
  if (igroup) error->all(FLERR,"Compute elastic must use group all");

  scalar_flag = vector_flag = 0;
  array_flag = 1;
  size_array_rows = 6;
  size_array_cols = 6;
  extarray = 0;
  timeflag = 1;

  // store pressure ID used by elastic computation
  // insure it is valid for pressure computation

  if (strcmp(arg[3],"NULL") == 0) id_press = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[3]);

    int icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute elastic pressure ID");
    if (modify->compute[icompute]->pressflag == 0)
      error->all(FLERR,"Compute elastic pressure ID does not "
                 "compute pressure");
  }

  // process optional args

  if (narg == 4) {
    pressureflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    fixflag = 1;
  } else {
    pressureflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"pressure") == 0) pressureflag = 1;
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
	fixflag = 1;
      } else error->all(FLERR,"Illegal compute elastic command");
      iarg++;
    }
  }

  // error check

  if (kspaceflag)
    error->all(FLERR,"Compute elastic does not support "
	       "kspace integration");
  
  if (pressureflag && id_press == NULL)
    error->all(FLERR,"Compute elastic requires pressure ID "
               "to include pressure tensor");

  array = new double*[6];
  for (int i = 0; i < 6; i++)
    array[i] = new double[6];
  nvirial = 0;
  v2ptr = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeElastic::~ComputeElastic()
{
  delete [] id_press;
  for (int i = 0; i < 6; i++)
    delete [] array[i];
  delete [] array;
  delete [] v2ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeElastic::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (pressureflag) {
    int icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute elastic pressure ID");
    pressure = modify->compute[icompute];
  }

  // detect contributions to virial
  // v2ptr points to all virial2[6][6] contributions

  delete [] v2ptr;
  nvirial = 0;
  v2ptr = NULL;

  if (pairflag && force->pair) nvirial++;
  if (bondflag && atom->molecular && force->bond) nvirial++;
  if (angleflag && atom->molecular && force->angle) nvirial++;
  if (dihedralflag && atom->molecular && force->dihedral) nvirial++;
  if (improperflag && atom->molecular && force->improper) nvirial++;
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->thermo_virial) nvirial++;

  if (nvirial) {
    v2ptr = new double**[nvirial];
    nvirial = 0;
    if (pairflag && force->pair) v2ptr[nvirial++] = force->pair->virial2;
    if (bondflag && force->bond) v2ptr[nvirial++] = force->bond->virial2;
    if (angleflag && force->angle) v2ptr[nvirial++] = force->angle->virial2;
    if (dihedralflag && force->dihedral)
      v2ptr[nvirial++] = force->dihedral->virial2;
    if (improperflag && force->improper)
      v2ptr[nvirial++] = force->improper->virial2;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->thermo_virial)
          v2ptr[nvirial++] = modify->fix[i]->virial2;
  }
}


/* ----------------------------------------------------------------------
   compute Born contribution to stiffness tensor in NVT
   assume pressure tensor has already been computed
------------------------------------------------------------------------- */

void ComputeElastic::compute_array()
{
  // invoke pressure if it hasn't been already

  double *pressure_tensor;
  if (pressureflag) {
    if (pressure->invoked_vector != update->ntimestep)
      pressure->compute_vector();
    pressure_tensor = pressure->vector;
  }

  // compute stiffness tensor in Voigt notation
  
  double e[6][6];
  int i,j,k,l,m,n,p,delta;
  for (i = 0; i < 6; i++) {
    if (i == 0) {
      k = 0;
      l = 0;
    }
    else if (i == 1) {
      k = 1;
      l = 1;
    }
    else if (i == 2) {
      k = 2;
      l = 2;
    }
    else if (i == 3) {
      k = 1;
      l = 2;
    }
    else if (i == 4) {
      k = 0;
      l = 2;
    }
    else {
      k = 0;
      l = 1;
    }
    for (j = 0; j < 6; j++) {
      if (j == 0) {
        m = 0;
        n = 0;
      }
      else if (j == 1) {
        m = 1;
        n = 1;
      }
      else if (j == 2) {
        m = 2;
        n = 2;
      }
      else if (j == 3) {
        m = 1;
        n = 2;
      }
      else if (j == 4) {
        m = 0;
        n = 2;
      }
      else {
        m = 0;
        n = 1;
      }

      if (l == n) delta = 1;
      else delta = 0;

      if (k == 0) {
        if (m == 0) p = 0;
	else if (m == 1) p = 3;
	else p = 4;
      }
      else if (k == 1) {
        if (m == 0) p = 3;
	else if (m == 1) p = 1;
	else p = 5;
      }
      else {
        if (m == 0) p = 4;
	else if (m == 1) p = 5;
	else p = 3;
      }

      inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
      e[i][j] = virial2[i][j]*inv_volume*nktv2p 
	      + delta*pressure_tensor[p];
    }

  }
}

/* ---------------------------------------------------------------------- */

void ComputeElastic::virial_compute()
{
  int i,j,k;
  double v2[6][6],**v2component;

  for (i = 0; i < 6; i++) 
    for (j = 0; j < 6; j++)
	  v2[i][j] = 0.0;

  // sum contributions to virial from forces and fixes

  for (k = 0; k < nvirial; k++) {
    v2component = v2ptr[k];
    for (i = 0; i < 6; i++)
      for (j = 0; j < 6; j++) v2[i][j] += v2component[i][j];
  }

  // sum virial across procs

  for (i = 0; i < 6; i++)
    MPI_Allreduce(v2[i],virial2[i],6,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void ComputeElastic::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_press;
  int n = strlen(id_new) + 1;
  id_press = new char[n];
  strcpy(id_press,id_new);
}
