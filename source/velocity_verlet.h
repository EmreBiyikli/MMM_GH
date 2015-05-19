// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// velocity_verlet.h
// *****************************************************************************
// Consists of VelocityVerlet class that handles time integration. In addition
// to the construction of the class, it needs to be initialized by the
// Init function. Time integration is performed by the Run function.

#ifndef MMM_V14_6_VELOCITY_VERLET_H_
#define MMM_V14_6_VELOCITY_VERLET_H_

#include <vector>

#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class VelocityVerlet : public Mediator {
 public:
  explicit VelocityVerlet() : Mediator() {}
  ~VelocityVerlet() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // VelocityVerlet(VelocityVerlet&) = delete;
  // VelocityVerlet& operator=(const VelocityVerlet&) = delete;

  // Initializes the class by initializing the input member variables.
  void Init(double timestep, int total_iteration);
  // Performs the time integration.
  void Run();

  // Accessor and mutator functions:
  // timestep_
  double timestep() const { return timestep_; }
  void set_timestep(double timestep) { timestep_ = timestep; }
  // total_iteration_
  int total_iteration() const { return total_iteration_; }
  void set_total_iteration(int total_iteration) {
    total_iteration_ = total_iteration;
  }

 private:
  // Applies displacement loading.
  void ApplyDisplacementLoad();
  // Checks whether to terminate the time integration.
  void CheckStopCriteria();
  // Displaces ghost atoms.
  void DisplaceGhostAtom();
  // Displaces rep atoms.
  void DisplaceRepAtom();
  // Performs an iteration.
  void Iterate();
  // Updates velocities of atoms.
  void UpdateVelocity();

  // Accessor and mutator functions:
  // timestep_half_
  double timestep_half() const { return timestep_half_; }
  void set_timestep_half(double timestep_half) {
    timestep_half_ = timestep_half;
  }
  // timestep_square_half_
  double timestep_square_half() const { return timestep_square_half_; }
  void set_timestep_square_half(double timestep_square_half) {
    timestep_square_half_ = timestep_square_half;
  }

  double timestep_;               // timestep
  double timestep_half_;          // half of timestep
  double timestep_square_half_;   // half of square of timestep
  int total_iteration_;           // total number of iterations
};

#endif  // MMM_V14_6_VELOCITY_VERLET_H_

