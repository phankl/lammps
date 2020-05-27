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

#include "fix_electric.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "neighbor.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;


FixElectric::FixElectric(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if (narg != 7) error->all(FLERR,"Illegal fix electric command");

  k = force->numeric(FLERR,arg[3]);
  evec[0] = force->numeric(FLERR,arg[4]);
  evec[1] = force->numeric(FLERR,arg[5]);
  evec[2] = force->numeric(FLERR,arg[6]);

  printf("Parameters: %e %e %e %e\n",k,evec[0],evec[1],evec[2]);

  eflag = 0;
  ee = 0.0;
}

/* ---------------------------------------------------------------------- */

FixElectric::~FixElectric()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

int FixElectric::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectric::init()
{
}

/* ---------------------------------------------------------------------- */

void FixElectric::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixElectric::post_force(int /*vflag*/)
{
  // apply electric force to each bond

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nbondlist = neighbor->nbondlist;
  int newton_bond = force->newton_bond;

  eflag = 0;
  ee = 0.0;

  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];

    double *x1 = x[i1];
    double *x2 = x[i2];

    double delta[3];
    sub3(x2,x1,delta);
    normalize3(delta,delta);
    
    double de = dot3(delta,evec);
    
    double f1[3];
    double f2[3];
    copy3(evec,f1);
    copy3(delta,f2);
    scale3(-2.0*k*de,f1);
    scaleadd3(k*de*de,f2,f1,f1);

    if (mask[i1] & groupbit) add3(f1,f[i1],f[i1]);
    if (mask[i2] & groupbit) sub3(f[i2],f1,f[i2]);  
  }
}

