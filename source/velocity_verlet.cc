// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// velocity_verlet.cc
// *****************************************************************************

#include "velocity_verlet.h"

#include <mpi.h>

#include <cmath>

#include "atom_group.h"
#include "matrix.h"
#include "mesh.h"
#include "mmm.h"
#include "model.h"
#include "neighbor.h"
#include "output.h"
#include "potential.h"
#include "temperature.h"

// Public
// *****************************************************************************

// Init
// *****************************************************************************
void VelocityVerlet::Init(double timestep, int total_iteration) {
  set_timestep(timestep);
  set_total_iteration(total_iteration);
  set_timestep_half(timestep / 2);
  set_timestep_square_half(pow(timestep, 2) / 2);
}

// Run
// *****************************************************************************
void VelocityVerlet::Run() {
  output()->OutputToLog("");
  output()->OutputToLog("velocity verlet run");
  output()->OutputToLog("dashed line");
  mmm()->set_current_iteration(-1);
  mmm()->set_is_iteration_loop_continue(true);
  // Loop
  while (mmm()->is_iteration_loop_continue()) {
    mmm()->increase_current_iteration();
    neighbor()->BuildNeighbor(atom_group()->position_.get(),
                              mmm()->current_iteration());
    if (mmm()->current_iteration() == 0) potential()->Compute(true, true);
    temperature()->Compute();
    mmm()->set_total_energy(mmm()->potential_energy() +
                            mmm()->kinetic_energy());
    if (mmm()->current_iteration() == 0) {
      if (temperature()->is_on()) {
        temperature()->Set();  // Set to initial temperature
      }
    } else {
      Iterate();
    }
    temperature()->Compute();
    mmm()->set_total_energy(mmm()->potential_energy() +
                                mmm()->kinetic_energy());
    if (mmm()->current_iteration() % output()->screen_update_frequency() == 0 ||
        !mmm()->is_iteration_loop_continue()) {
      if (mmm()->current_iteration() == 0) {
        string temp = "iteration total_energy kinetic_energy potential_";
        temp += "energy temperature";
        output()->OutputToLog(temp);
        output()->OutputToLog("dashed line");
      }
      output()->OutputToLog(output()->ToString(mmm()->current_iteration()) +
          " " + output()->ToString(mmm()->total_energy()) +
          " " + output()->ToString(mmm()->kinetic_energy()) +
          " " + output()->ToString(mmm()->potential_energy()) +
          " " + output()->ToString(temperature()->current_temperature()));
    }
    output()->OutputAllToFile();
    mmm()->CheckLostAtom();
  }
}

// Private
// *****************************************************************************

// ApplyDisplacementLoad
// *****************************************************************************
void VelocityVerlet::ApplyDisplacementLoad() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i) ||
        !atom_group()->is_load_.get_1d_element(i)) 
        continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->load_type_.get_element(i, j) == 1) {
        atom_group()->displacement_.set_element(i, j,
            atom_group()->load_value_.get_element(i, j));
      } else if (atom_group()->load_type_.get_element(i, j) == 2) {
        atom_group()->displacement_.add_element(i, j,
            atom_group()->load_value_.get_element(i, j));
      }
    }
  }
}

// CheckStopCriteria
// *****************************************************************************
void VelocityVerlet::CheckStopCriteria() {
  if (mmm()->current_iteration() == total_iteration()) {
    mmm()->set_is_iteration_loop_continue(false);
  }
}

// DisplaceGhostAtom
// *****************************************************************************
void VelocityVerlet::DisplaceGhostAtom() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      atom_group()->displacement_.set_element(i, j, 0);
    }
    int element = model()->atom_element_.get_element(i, 0);
    for (int j = 0; j < mesh()->element_node_num(); j++) {
      int node = mesh()->element_node_at(element, j);
      double interpolation_weight = model()->shape_function_.get_element(i, j);
      for (int k = 0; k < mmm()->dimension(); k++) {
        atom_group()->displacement_.add_element(i, k, interpolation_weight * 
            atom_group()->displacement_.get_element(node, k));
      }
    }
  }
}

// DisplaceRepAtom
// *****************************************************************************
void VelocityVerlet::DisplaceRepAtom() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      atom_group()->displacement_.set_element(i, j,
          atom_group()->velocity_.get_element(i, j) * timestep() +
          (atom_group()->applied_force_.get_element(i, j) /
          model()->atom_mass_.get_1d_element(i)) * timestep_square_half());
    }
  }
}

// Iterate
// *****************************************************************************
void VelocityVerlet::Iterate() {
  DisplaceRepAtom();
  ApplyDisplacementLoad();
  if (!model()->is_full_atomistic()) DisplaceGhostAtom();
  atom_group()->position_.add(atom_group()->displacement_.get());
  atom_group()->old_applied_force_.set(atom_group()->applied_force_.get());
  potential()->Compute(true, true);
  if (temperature()->is_on()) {  // Apply thermostat (if applied on force)
    temperature()->ApplyThermostat("after_force");
  }
  UpdateVelocity();
  if (temperature()->is_on()) {  // Apply thermostat (if applied on velocity)
    temperature()->ApplyThermostat("after_velocity");
  }
  CheckStopCriteria();
}

// UpdateVelocity
// *****************************************************************************
void VelocityVerlet::UpdateVelocity() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      atom_group()->velocity_.add_element(i, j,
          ((atom_group()->applied_force_.get_element(i, j) +
              atom_group()->old_applied_force_.get_element(i, j)) /
              model()->atom_mass_.get_1d_element(i)) * timestep_half());   
    }
  }
}

