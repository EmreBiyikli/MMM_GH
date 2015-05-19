// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// temperature.cc
// *****************************************************************************

#include "temperature.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "atom_group.h"
#include "matrix.h"
#include "model.h"
#include "mmm.h"
#include "velocity_verlet.h"

using std::string;
using std::vector;

// Public
// *****************************************************************************

// Compute
// *****************************************************************************
// Note that this function ignored loads on atoms.
void Temperature::Compute() {
  // Kinetic energy
  double degree_of_freedom_num = 0;
  mmm()->set_kinetic_energy(0);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->is_load_.get_1d_element(i) &&
          (atom_group()->load_type_.get_element(i, j) == 1 ||
           atom_group()->load_type_.get_element(i, j) == 3)) {
        continue;
      }
      degree_of_freedom_num++;
      mmm()->add_kinetic_energy(0.5 * model()->atom_mass_.get_1d_element(i) *
                                    atom_group()->velocity_.get_element(i, j) *
                                    atom_group()->velocity_.get_element(i, j));
    }
  }
  // Temperature
  set_current_temperature(2 * mmm()->kinetic_energy() / (degree_of_freedom_num *
                              kBoltzmannConst));
}

// ApplyThermostat
// *****************************************************************************
// State indicates when the thermostat is called.
void Temperature::ApplyThermostat(string state) {
  if (state.compare("after_force") == 0) {
    if (style().compare("berendsen") == 0) {
      ApplyBerendsen();
    } else if (style().compare("langevin") == 0) {
      ApplyLangevin();
    }
  } else if (state.compare("after_velocity") == 0) {
    if (style().compare("velocity_rescaling") == 0) {
      ApplyVelocityRescaling();
    }
  }
}

// InitBerendsen
// *****************************************************************************
void Temperature::InitBerendsen(double dissipation_coefficient,
                                double initial_temperature,
                                double target_temperature) {
  set_dissipation_coefficient(dissipation_coefficient);
  set_initial_temperature(initial_temperature);
  set_target_temperature(target_temperature);
  set_style("berendsen");
  set_is_on(true);
}

// InitLangevin
// *****************************************************************************
void Temperature::InitLangevin(double dissipation_coefficient,
                               double initial_temperature,
                               double target_temperature) {
  set_dissipation_coefficient(dissipation_coefficient);
  set_initial_temperature(initial_temperature);
  set_target_temperature(target_temperature);
  set_style("langevin");
  set_is_on(true);
}

// InitVelocityRescaling
// *****************************************************************************
void Temperature::InitVelocityRescaling(double initial_temperature,
                                        double target_temperature,
                                        int update_frequency) {
  set_initial_temperature(initial_temperature);
  set_target_temperature(target_temperature);
  set_update_frequency(update_frequency);
  set_style("velocity_rescaling");
  set_is_on(true);
}

// Set
// *****************************************************************************
// Set temperature by setting atom velocities to number drawn from a normal
// (Gaussian) distribution with 0 mean and sqrt(Boltzmann_constant *
// initial_temperature / mass).
// Note that this function ignored loads on atoms.
void Temperature::Set() {
  double double_RAND_MAX_plus_one = static_cast<double>(RAND_MAX + 1.0);
  double mean = 0;
  double standard_deviation = sqrt(kBoltzmannConst * initial_temperature() /
                                       mmm()->mass());
  double velocity_sum_count = 0;
  vector<double> velocity_sum(3, 0);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->is_load_.get_1d_element(i) &&
          (atom_group()->load_type_.get_element(i, j) == 1 ||
            atom_group()->load_type_.get_element(i, j) == 3)) {
        continue;
      }
      // This gives equal distribution in (0,1].
      double random_num_1 = static_cast<double>(rand() + 1.0) /
                                double_RAND_MAX_plus_one;
      double random_num_2 = static_cast<double>(rand() + 1.0) /
                                double_RAND_MAX_plus_one;
      // Below, sqrt and cos are for conversion from uniform to normal
      // distribution.
      double random_num = mean + standard_deviation *
          sqrt(-2 * log(random_num_1)) * cos(2 * kPi * random_num_2);
      atom_group()->velocity_.set_element(i, j, random_num);
      velocity_sum[j] += random_num;
      velocity_sum_count++;
    }
  }
  // Assure a velocity of 0 mean.
  velocity_sum_count /= mmm()->dimension();
  for (int i = 0; i < mmm()->dimension(); i++) {
    velocity_sum[i] /= velocity_sum_count;
  }
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;    
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->is_load_.get_1d_element(i) &&
          (atom_group()->load_type_.get_element(i, j) == 1 ||
            atom_group()->load_type_.get_element(i, j) == 3)) {
        continue;
      }
      atom_group()->velocity_.add_element(i, j, -velocity_sum[j]);
    }
  }
}

// Private
// *****************************************************************************

// ApplyBerendsen
// *****************************************************************************
// It is applied on the force.
// For further details, see Dr. To's lecture notes.
void Temperature::ApplyBerendsen() {
  double scaling_factor = 1.0 - target_temperature() / current_temperature();
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->is_load_.get_1d_element(i) &&
          (atom_group()->load_type_.get_element(i, j) == 1 ||
           atom_group()->load_type_.get_element(i, j) == 3)) {
        continue;
      }
      atom_group()->force_.add_element(i, j,
          -dissipation_coefficient() * model()->atom_mass_.get_1d_element(i) *
          atom_group()->velocity_.get_element(i, j) * scaling_factor);
      atom_group()->applied_force_.add_element(i, j,
          -dissipation_coefficient() * model()->atom_mass_.get_1d_element(i) *
          atom_group()->velocity_.get_element(i, j) * scaling_factor);
    }
  }
}

// ApplyLangevin
// *****************************************************************************
// It is applied on the force.
// For further details, see corresponding file in LAMMPS.
// The first number in Fr_const is 24 in LAMMPS = 2 from theory * 12 from
// uniform distribution. However, for us, 12 worked better in some cases. I am
// not sure of the reason.
void Temperature::ApplyLangevin() {
  double double_RAND_MAX = static_cast<double>(RAND_MAX);
  double Fr_const = sqrt(24 * kBoltzmannConst * target_temperature() / 
      (velocity_verlet()->timestep() * dissipation_coefficient()));
  int seed = static_cast<int>(fmod(MPI_Wtime() / MPI_Wtick(), 1e6));
  srand(seed);  // It is very important to have high resolution here. 
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (atom_group()->is_load_.get_1d_element(i) &&
          (atom_group()->load_type_.get_element(i, j) == 1 ||
           atom_group()->load_type_.get_element(i, j) == 3)) {
        continue;
      }
      double Ff = -1.0 * model()->atom_mass_.get_1d_element(i) *
          atom_group()->velocity_.get_element(i, j) / dissipation_coefficient();
      double random_num = static_cast<double>(rand()) / double_RAND_MAX - 0.5;
      double Fr = random_num * sqrt(model()->atom_mass_.get_1d_element(i)) *
                                        Fr_const;
      atom_group()->force_.set_element(i, j, Ff + Fr);
      atom_group()->applied_force_.set_element(i, j, Ff + Fr);
    }
  }
}

// ApplyVelocityRescaling
// *****************************************************************************
// It is applied on the velocity.
void Temperature::ApplyVelocityRescaling() {
  if (mmm()->current_iteration() == 1 ||
      mmm()->current_iteration() % update_frequency() == 0) {
    Compute();
    double scaling_factor = sqrt(target_temperature() / current_temperature());
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      if (!model()->is_rep_.get_1d_element(i)) continue;
      for (int j = 0; j < mmm()->dimension(); j++) {
        if (atom_group()->is_load_.get_1d_element(i) &&
            (atom_group()->load_type_.get_element(i, j) == 1 ||
             atom_group()->load_type_.get_element(i, j) == 3)) {
          continue;
        }
        atom_group()->velocity_.multiply(i, j, scaling_factor);
      }
    }
  }
}

