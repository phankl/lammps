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
  if (narg != 5) error->all(FLERR,"Illegal fix electric command");

  n = force->numeric(arg[0],0,NULL);
  k = force->numeric(arg[1],0,NULL);
  evec[0] = force->numeric(arg[2],0,NULL);
  evec[1] = force->numeric(arg[3],0,NULL);
  evec[2] = force->numeric(arg[4],0,NULL);

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
  // apply electric force to each particle

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
    
    double de = dot3(delta,evec);
    double dd = dot3(delta,delta);

    double f1[3] = {0.0,0.0,0.0};
    scaleadd3(dd,evec,f1,f1);
    scaleadd3(-de,delta,f1,f1);

    double a = k * n * pow(dd,-1.0-n/2.0) * pow(de,n-1);
    scale3(a,f1);

    if (newton_bond || i1 < nlocal) add3(f1,f[i1],f[i1]);
    if (newton_bond || i2 < nlocal) sub3(f[i2],f1,f[i2]);
  }
}
