/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(virial/squared,ComputeVirialSquared)

#else

#ifndef LMP_COMPUTE_VIRIAL_SQUARED_H
#define LMP_COMPUTE_VIRIAL_SQUARED_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVirialSquared : public Compute {
 public:
  ComputeVirialSquared(class LAMMPS *, int, char **);
  ~ComputeVirialSquared();
  void init();
  void compute_array();
  void reset_extra_compute_fix(const char *);

 protected:
  Compute *virial;
  char *id_virial;

  void voigt_index(int, int &, int &);
  int virial_index(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute virial must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

E: Must use 'kspace_modify virial/scalar no' for tensor components with kspace_style msm

Otherwise MSM will compute only a scalar virial.  See the kspace_modify
command for details on this setting.

*/
