// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// potential.h
// *****************************************************************************
// Consists of Potential class that handles a pairwise potential. In addition to
// the construction of the class, it needs to be initialized by the Init
// function. Computation of energy and/or force is performed by Compute
// function.
// The file also includes pair potentials of PairBase as a base and some others
// PairLennardJones, PairMorse, and PairSpring polymorph from the base.

#ifndef MMM_V14_6_POTENTIAL_H_
#define MMM_V14_6_POTENTIAL_H_

#include <cmath>
#include <string>
#include <vector>

#include "mediator.h"

using std::string;
using std::vector;

// PairBase
// *****************************************************************************
// See comment at top of file for a complete description.
class PairBase {
 public:
  PairBase() {}
  virtual ~PairBase() {}
  virtual void Init(int, vector<double>) {}
  virtual double Energy(double) { return 0; }
  virtual void Force(double, vector<double>, vector<double>*) {}
  int dimension_;
};

// PairLennardJones
// *****************************************************************************
// See comment at top of file for a complete description.
class PairLennardJones : public PairBase {
 public:
  // using PairBase::PairBase;  // TODO: activate with C++11 compiler
  PairLennardJones() {}
  virtual ~PairLennardJones() {}

  virtual void Init(int dimension, vector<double> parameter) {
    dimension_ = dimension;
    epsilon_ = parameter[0];
    sigma_ = parameter[1];
    A_ = 4 * epsilon_ * pow(sigma_, 12);
    B_ = 4 * epsilon_ * pow(sigma_, 6);
  }
  virtual double Energy(double distance_magnitude) {
    double temp = pow(distance_magnitude, 6);
    return A_ / (temp * temp) - B_ / temp;
  }
  virtual void Force(double distance_magnitude, vector<double> distance_vector,
                     vector<double>* force) {
    double temp = -6 * B_ / pow(distance_magnitude, 8) + 12 * A_ /
        pow(distance_magnitude, 14);
    for (int i = 0; i < dimension_; i++) {
      (*force)[i] = temp * distance_vector[i];
    }
  }

 private:
  double A_, B_, epsilon_, sigma_;
};

// PairMorse
// *****************************************************************************
// See comment at top of file for a complete description.
class PairMorse : public PairBase {
 public:
  // using PairBase::PairBase;  // TODO: activate with C++11 compiler
  PairMorse() {}
  virtual ~PairMorse() {}

  virtual void Init(int dimension, vector<double> parameter) {
    dimension_ = dimension;
    D0_ = parameter[0];
    alpha_ = parameter[1];
    initial_spacing_ = parameter[2];
  }
  virtual double Energy(double distance_magnitude) {
    double temp = exp(-alpha_ * (distance_magnitude - initial_spacing_));
    return D0_ * (temp * temp - 2 * temp);
  }
  virtual void Force(double distance_magnitude, vector<double> distance_vector,
                     vector<double>* force) {
    double temp_1 = exp(-alpha_ * (distance_magnitude - initial_spacing_));
    double temp_2 = 2 * alpha_ * D0_ * (temp_1 * temp_1 - temp_1) /
        distance_magnitude;
    for (int i = 0; i < dimension_; i++) {
      (*force)[i] = temp_2 * distance_vector[i];
    }
  }

 private:
  double alpha_, D0_, initial_spacing_;
};

// PairSpring
// *****************************************************************************
// See comment at top of file for a complete description.
class PairSpring : public PairBase {
 public:
  // using PairBase::PairBase;  // TODO: activate with C++11 compiler
  PairSpring() {}
  virtual ~PairSpring() {}

  virtual void Init(int dimension, vector<double> parameter) {
    dimension_ = dimension;
    initial_spacing_ = parameter[0];
    stiffness_ = parameter[1];
  }
  virtual double Energy(double distance_magnitude) {
    double displacement_difference = distance_magnitude - initial_spacing_;
    return 0.5 * stiffness_ * displacement_difference * displacement_difference;
  }
  virtual void Force(double distance_magnitude, vector<double> distance_vector,
                     vector<double>* force) {
    double distance_difference = distance_magnitude - initial_spacing_;
    double temp = -stiffness_ * distance_difference / distance_magnitude;
    for (int i = 0; i < dimension_; i++) {
      (*force)[i] = temp * distance_vector[i];
    }
  }

 private:
  double initial_spacing_, stiffness_;
};

// Potential
// *****************************************************************************
// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Potential : public Mediator {
 public:
  explicit Potential() : Mediator() {
    distance_vector_.resize(3);
  }
  ~Potential() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Potential(Potential&) = delete;
  // Potential& operator=(const Potential&) = delete;

  // Computes energy and/or force with respect to the input choices.
  void Compute(bool is_energy_compute, bool is_force_compute);
  // Initializes the class by initializing the input member variables.
  void Init(const vector<double>& parameter, double cutoff_radius,
            string potential);

  // Accessor and mutator functions:
  // cutoff_radius_
  double cutoff_radius() const { return cutoff_radius_; }
  void set_cutoff_radius(double cutoff_radius) {
    cutoff_radius_ = cutoff_radius;
  }
  // cutoff_radius_square_
  double cutoff_radius_square() const { return cutoff_radius_square_; }
  void set_cutoff_radius_square(double cutoff_radius_square) {
    cutoff_radius_square_ = cutoff_radius_square;
  }
  // distance_magnitude_
  double distance_magnitude() const { return distance_magnitude_; }
  void set_distance_magnitude(double distance_magnitude) {
    distance_magnitude_ = distance_magnitude;
  }
  // distance_magnitude_square_
  double distance_magnitude_square() const {
    return distance_magnitude_square_;
  }
  void set_distance_magnitude_square(double distance_magnitude_square) {
    distance_magnitude_square_ = distance_magnitude_square;
  }
  // distance_vector_
  const vector<double>& distance_vector() { return distance_vector_; }
  double distance_vector(int dimension) const {
    return distance_vector_[dimension];
  }
  void set_distance_vector(int dimension, double distance) {
    distance_vector_[dimension] = distance;
  }
  // equilibrium_spacing_
  double equilibrium_spacing() const { return equilibrium_spacing_; }
  void set_equilibrium_spacing(double equilibrium_spacing) {
    equilibrium_spacing_ = equilibrium_spacing;
  }
  // first_atom_
  int first_atom() const { return first_atom_; }
  void set_first_atom(int first_atom) { first_atom_ = first_atom; }
  // is_energy_compute_
  bool is_energy_compute() const { return is_energy_compute_; }
  void set_is_energy_compute(bool is_energy_compute) {
    is_energy_compute_ = is_energy_compute;
  }
  // is_force_compute_
  bool is_force_compute() const { return is_force_compute_; }
  void set_is_force_compute(bool is_force_compute) {
    is_force_compute_ = is_force_compute;
  }
  // parameter_
  double parameter(int index) const { return parameter_[index]; }
  void set_parameter(const vector<double>& parameter) {
    parameter_ = parameter;
  }
  // potential_
  string potential() const { return potential_; }
  void set_potential(string potential) { potential_ = potential; }
  // second_atom_
  int second_atom() const { return second_atom_; }
  void set_second_atom(int second_atom) { second_atom_ = second_atom; }

 private:
  // Applies loading.
  void ApplyLoad();
  // Shares information between processing elements.
  void Communicate();
  // Computes the distance between current first_atom and second_atom.
  void ComputeDistance();
  // Computes the energy between current first_atom and second_atom.
  void ComputeEnergy();
  // Computes the force between current first_atom and second_atom.
  void ComputeForce();
  // Computes the total potential energy.
  void ComputeTotalPotentialEnergy();
  // Loops atoms.
  void LoopAtom();
  // Manipulates energy of atoms with respect to the MMM scheme, if MMM is on.
  void MMMEnergy();
  // Manipulates force of atoms with respect to the MMM scheme, if MMM is on.
  void MMMForce();
  // Resets energies and forces of atoms to 0.
  void ResetEnergyForce();

  bool is_energy_compute_;              // whether to compute the energy
  bool is_first_time_;                  // whether it is the first time
  bool is_force_compute_;               // whether to compute the force
  double distance_magnitude_;           // magnitude of the distance
  double distance_magnitude_square_;    // square of the magnitude of distance
  double equilibrium_spacing_;          // equilibrium spacing
  double cutoff_radius_;                // cut-off distance
  double cutoff_radius_square_;         // square of the cut-off distance
  int first_atom_;                      // current first atom
  int second_atom_;                     // current second atom
  string potential_;                    // potential
  vector<double> distance_vector_;      // distance vector
  vector<double> parameter_;            // potential parameters

  PairBase* pair_potential;             // base for pair potentials
  PairLennardJones pair_lennard_jones;  // pair lennard-jones potential
  PairMorse pair_morse;                 // pair morse potential
  PairSpring pair_spring;               // pair spring potential
};

#endif  // MMM_V14_6_POTENTIAL_H_

