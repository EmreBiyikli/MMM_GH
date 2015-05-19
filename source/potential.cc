// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// potential.cc
// *****************************************************************************
// classes are placed after the functions

#include "potential.h"

#include <mpi.h>

#include <string>
#include <vector>

#include "atom_group.h"
#include "matrix.h"
#include "mesh.h"
#include "mmm.h"
#include "model.h"
#include "neighbor.h"
#include "time.h"

using std::string;
using std::vector;

// Public
// *****************************************************************************

// Compute
// *****************************************************************************
void Potential::Compute(bool is_energy_compute, bool is_force_compute) {
  set_is_energy_compute(is_energy_compute);
  set_is_force_compute(is_force_compute);
  ResetEnergyForce();
  double start_time = MPI_Wtime();
  LoopAtom();
  double finish_time = MPI_Wtime();
  time()->add_parallel_time(finish_time - start_time);
  if (mmm()->is_MPI()) Communicate();
  if (!model()->is_full_atomistic()) {
    if (is_energy_compute) MMMEnergy();
    if (is_force_compute) MMMForce();
  }
  if (is_energy_compute) ComputeTotalPotentialEnergy();
  if (is_force_compute) {
    atom_group()->applied_force_.set(atom_group()->force_.get());
    ApplyLoad();
  }
}

// Init
// *****************************************************************************
void Potential::Init(const vector<double>& parameter, double cutoff_radius,
                     string potential) {
  set_cutoff_radius(cutoff_radius);
  set_potential(potential);
  set_parameter(parameter);
  if (potential.compare("lennard_jones") == 0) {
    pair_potential = &pair_lennard_jones;
    set_equilibrium_spacing(pow(2.0, 1.0/6.0) * parameter[1]);
  } else if (potential.compare("morse") == 0) {
    pair_potential = &pair_morse;
    set_equilibrium_spacing(parameter[2]);
  } else if (potential.compare("spring") == 0) {
    pair_potential = &pair_spring;
    set_equilibrium_spacing(parameter[0]);
  } 
  pair_potential->Init(mmm()->dimension(), parameter);
  set_cutoff_radius_square(pow(cutoff_radius, 2));
}

// Private
// *****************************************************************************

// ApplyForceLoad
// *****************************************************************************
void Potential::ApplyLoad() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!atom_group()->is_load_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->load_type_.get_element(i, j) == 1) {
        atom_group()->applied_force_.set_element(i, j, 0.0);
      } else if (atom_group()->load_type_.get_element(i, j) == 3) {
        atom_group()->applied_force_.set_element(i, j,
            atom_group()->load_value_.get_element(i, j));
      } else if (atom_group()->load_type_.get_element(i, j) == 4) {
        atom_group()->applied_force_.add_element(i, j,
            atom_group()->load_value_.get_element(i, j));
      }
    }
  }
}

// Communicate
// *****************************************************************************
void Potential::Communicate() {
  MPI_Barrier(MPI_COMM_WORLD);
  if (is_energy_compute()) {
    MPI_Allreduce(MPI_IN_PLACE, atom_group()->energy_.address(),
                  atom_group()->atom_num(), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
  if (is_force_compute()) {
    MPI_Allreduce(MPI_IN_PLACE, atom_group()->force_.address(),
                  atom_group()->degree_of_freedom_num(), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// ComputeDistance
// *****************************************************************************
void Potential::ComputeDistance() {      
  set_distance_vector(0, atom_group()->position_.get_element(first_atom(), 0) -
                      atom_group()->position_.get_element(second_atom(), 0));
  set_distance_vector(1, atom_group()->position_.get_element(first_atom(), 1) -
                      atom_group()->position_.get_element(second_atom(), 1));
  set_distance_vector(2, atom_group()->position_.get_element(first_atom(), 2) -
                      atom_group()->position_.get_element(second_atom(), 2));
  set_distance_magnitude_square(pow(distance_vector(0), 2) +
                                pow(distance_vector(1), 2) +
                                pow(distance_vector(2), 2));
}

// ComputeEnergy
// *****************************************************************************
void Potential::ComputeEnergy() {
  double energy = pair_potential->Energy(distance_magnitude()) / 2.0;
  atom_group()->energy_.add_1d_element(first_atom(), energy);
  atom_group()->energy_.add_1d_element(second_atom(), energy);
}

// ComputeForce
// *****************************************************************************
void Potential::ComputeForce() {
  vector<double> force(3, 0);
  pair_potential->Force(distance_magnitude(), distance_vector(), &force);
  if (!model()->is_full_atomistic()) {
    double force_weight = (model()->force_weight_.get_1d_element(first_atom()) +
        model()->force_weight_.get_1d_element(second_atom())) / 2;
    force[0] *= force_weight;
    force[1] *= force_weight;
    force[2] *= force_weight;
  }
  for (int k = 0; k < mmm()->dimension(); k++) {
    atom_group()->force_.add_element(first_atom(), k, force[k]);
    atom_group()->force_.add_element(second_atom(), k, -force[k]);
  }
}

// ComputeTotalPotentialEnergy
// *****************************************************************************
void Potential::ComputeTotalPotentialEnergy() {
  mmm()->set_potential_energy(0.0);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    mmm()->add_potential_energy(atom_group()->energy_.get_1d_element(i));
  }
}

// LoopAtom
// *****************************************************************************
void Potential::LoopAtom() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    set_first_atom(i);
    if (atom_group()->processing_element_.get_1d_element(first_atom()) !=
        mmm()->processing_element()) {
      continue;
    }
    for (int j = 0; j < neighbor()->neighbor_num(first_atom()); j++) {
      set_second_atom(neighbor()->neighbor_list_at(first_atom(), j));
      if (second_atom() < first_atom()) continue;
      if (model()->atom_type_.get_1d_element(first_atom()) == 5 &&
          model()->atom_type_.get_1d_element(second_atom()) == 5) {
        continue;
      }
      vector<double> distance_vector(3, 0);
      ComputeDistance();
      if (distance_magnitude_square() > cutoff_radius_square()) continue;
      set_distance_magnitude(sqrt(distance_magnitude_square()));
      if (is_energy_compute()) ComputeEnergy();
      if (is_force_compute()) ComputeForce();
    }
  }
}

// MMMEnergy
// *****************************************************************************
// Sets energy of nsmp atoms to the energy of corresponding ssmp atom.
void Potential::MMMEnergy() {
  for (int i = 0; i < mesh()->element_num(); i++) {
    if (model()->element_atom_num_.get_1d_element(i) == 0) continue;    
    int psmp = model()->element_center_atom_.get_1d_element(i);
    double element_energy = atom_group()->energy_.get_1d_element(psmp);
    for (int j = 0; j < model()->element_atom_num_.get_1d_element(i); j++) {
      int atom = model()->element_atom_.get_element(i, j);
      if (model()->atom_type_.get_1d_element(atom) == 5) {
        atom_group()->energy_.set_1d_element(atom, element_energy);
      }
    }
  }
}

// MMMForce
// *****************************************************************************
// Adds force on ghost atoms to irep atoms.
void Potential::MMMForce() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (model()->is_rep_.get_1d_element(i)) continue;
    int element = model()->atom_element_.get_element(i, 0);
    for (int j = 0; j < mesh()->element_node_num(); j++) {
      int node = mesh()->element_node_at(element, j);
      double interpolation_weight = model()->shape_function_.get_element(i, j);      
      for (int k = 0; k < mmm()->dimension(); k++) {
        atom_group()->force_.add_element(node, k, interpolation_weight *
            atom_group()->force_.get_element(i, k));
      }
    }
  }
}

// ResetEnergyForce
// *****************************************************************************
void Potential::ResetEnergyForce() {
  if (is_energy_compute()) atom_group()->energy_.reset();
  if (is_force_compute()) atom_group()->force_.reset();
}

