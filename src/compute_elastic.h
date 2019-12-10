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

ComputeStyle(elastic,ComputeElastic)

#else

#ifndef LMP_COMPUTE_ELASTIC_H
#define LMP_COMPUTE_ELASTIC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeElastic : public Compute {
 public:
  ComputeElastic(class LAMMPS *, int, char **);
  ~ComputeElastic();
  void init();
  void compute_array();
  void virial_compute();
  void reset_extra_compute_fix(const char *);

 protected:
  double boltz,nktv2p,inv_volume;
  int nvirial,dimension;
  double virial2[6][6];
  double ***v2ptr;
  Compute *pressure;
  char *id_press;
  int pressureflag;
  int pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int fixflag,kspaceflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute elastic must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Could not find compute elastic pressure ID

The compute ID for calculating pressure does not exist.

E: Compute elastic pressure ID does not compute pressure

The compute ID assigned to an elastic computation must compute
pressure.

E: Compute elastic requires pressure ID to include stress term

The keflag cannot be used unless a temperature compute is provided.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
