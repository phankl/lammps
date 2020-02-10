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
  Fix(lmp, narg, arg),
{
  if (narg != 5) error->all(FLERR,"Illegal fix electric command");

  int num = sscanf(arg,"%d %lg %lg %lg %lg",
                    &n,&k,&evec[0],&evec[1],&evec[2]);
  if (num != 5) 
    error->one(FLERR,"Incorrect format of fix electric parameters");

  eflag = 0;
  ee = 0.0;
}

