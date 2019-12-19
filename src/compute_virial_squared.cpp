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
#include "compute_virial_squared.h"
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

ComputeVirialSquared::ComputeVirialSquared(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute virial/squared command");
  if (igroup) error->all(FLERR,"Compute virial/squared must use group all");

  array_flag = 1;
  size_array_rows = 6;
  size_array_cols = 6;
  extarray = 0;
  pressflag = 0;
  timeflag = 1;

  // store virial ID used by elastic computation
  // insure it is valid for pressure computation

  if (strcmp(arg[3],"NULL") == 0) id_virial = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    id_virial = new char[n];
    strcpy(id_virial,arg[3]);

    int icompute = modify->find_compute(id_virial);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute virial/squared virial ID");
    if (modify->compute[icompute]->pressflag == 0)
      error->all(FLERR,"Compute virial/squared virial ID does not "
                 "compute pressure");
  }
  // process optional args

  if (narg != 4) error->all(FLERR,"Illegal compute virial/squared command");
  
  array = new double*[6];
  for (int i = 0; i < 6; i++)
    array[i] = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeVirialSquared::~ComputeVirialSquared()
{
  for (int i = 0; i < 6; i++)
    delete [] array[i];
  delete [] array;
}

/* ---------------------------------------------------------------------- */

void ComputeVirialSquared::init()
{
  // set virial compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  int icompute = modify->find_compute(id_virial);
  if (icompute < 0)
    error->all(FLERR,"Could not find compute virial/squared virial ID");
  virial = modify->compute[icompute];
}

/* ----------------------------------------------------------------------
   compute squared virial tensor
------------------------------------------------------------------------- */

void ComputeVirialSquared::compute_array()
{
  invoked_array = update->ntimestep;

  // invoke virial if it hasn't been already

  double *virial_tensor;
  if (virial->invoked_vector != update->ntimestep)
    virial->compute_vector();
  virial_tensor = virial->vector;

  // compute squared virial tensor

  int i,j,k,l,m,n,p,q;
  double virial_squared[3][3][3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++) {
	        p = virial_index(i,j);
	        q = virial_index(k,l);
	        virial_squared[i][j][k][l] = virial_tensor[p] * virial_tensor[q];
	      }

  // transform to Voigt notation
  
  for (i = 0; i < 6; i++) {
    voigt_index(i,k,l);
    for (j = 0; j < 6; j++) {
      voigt_index(j,m,n);
      array[i][j] = virial_squared[k][l][m][n];
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVirialSquared::voigt_index(int i, int &j, int &k)
{
  if (i == 0) {
    j = 0;
    k = 0;
  }
  else if (i == 1) {
    j = 1;
    k = 1;
  }
  else if (i == 2) {
    j = 2;
    k = 2;
  }
  else if (i == 3) {
    j = 1;
    k = 2;
  }
  else if (i == 4) {
    j = 0;
    k = 2;
  }
  else {
    j = 0;
    k = 1;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeVirialSquared::virial_index(int i, int j)
{
  if (i == 0) {
    if (j == 0) return 0;
    else if (j == 1) return 3;
    else return 4;
  }
  else if (i == 1) {
    if (j == 0) return 3;
    else if (j == 1) return 1;
    else return 5;
  }
  else {
    if (j == 0) return 4;
    else if (j == 1) return 5;
    else return 2;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVirialSquared::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_virial;
  int n = strlen(id_new) + 1;
  id_virial = new char[n];
  strcpy(id_virial,id_new);
}
