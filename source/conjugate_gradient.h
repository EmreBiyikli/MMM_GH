// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// conjugate_gradient.h
// *****************************************************************************
// Consists of ConjugateGradient class that handles energy minimization. In
// addition to the construction of the class, it needs to be initialized by the
// Init function. Minimization is performed when the Run function is called.

#ifndef MMM_V14_6_CONJUGATE_GRADIENT_H_
#define MMM_V14_6_CONJUGATE_GRADIENT_H_

#include <vector>

#include "matrix.h"
#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class ConjugateGradient : public Mediator {
 public:
  explicit ConjugateGradient() : Mediator() {}
  ~ConjugateGradient() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // ConjugateGradient(ConjugateGradient&) = delete;
  // ConjugateGradient& operator=(const ConjugateGradient&) = delete;

  // Initializes the class by initializing the input member variables.
  void Init(double max_displacement, double tolerance, int max_iteration);
  // Performs the minimization.
  void Run();

  // Accessor and mutator functions:
  // max_displacement_
  double max_displacement() const { return max_displacement_; }
  void set_max_displacement(double max_displacement) {
    max_displacement_ = max_displacement;
  }
  // max_iteration_
  int max_iteration() const { return max_iteration_; }
  void set_max_iteration(int max_iteration) { max_iteration_ = max_iteration; }
  // tolerance_
  double tolerance() const { return tolerance_; }
  void set_tolerance(double tolerance) { tolerance_ = tolerance; }

 private:
  // Applies displacement loading.
  void ApplyDisplacementLoad();
  // Checks whether to terminate the line search.
  void CheckLineSearchStopCriteria();
  // Check whether to terminate the minimization.
  void CheckStopCriteria();
  // Computes beta.
  void ComputeBeta();
  // Computes the direction to take a step, force_dot_direction, and
  // max_direction.
  void ComputeDirection();
  // Displaces ghost atoms.
  void DisplaceGhostAtom();
  // Displaces rep atoms.
  void DisplaceRepAtom();
  // Performs an iteration.
  void Iterate();
  // Performs line search.
  void LineSearch();

  // Accessor and mutator functions:
  // alpha_
  double alpha() const { return alpha_; }
  void set_alpha(double alpha) { alpha_ = alpha; }
  void multiply_alpha(double alpha) { alpha_ *= alpha; }
  // beta_
  double beta() const { return beta_; }
  void set_beta(double beta) { beta_ = beta; }
  // energy_check_
  double energy_check() const { return energy_check_; }
  void set_energy_check(double energy_check) { energy_check_ = energy_check; }
  // force_dot_direction_
  double force_dot_direction() const { return force_dot_direction_; }
  void set_force_dot_direction(double force_dot_direction) {
    force_dot_direction_ = force_dot_direction;
  }
  void add_force_dot_direction(double force_dot_direction) {
    force_dot_direction_ += force_dot_direction;
  }
  // is_continue_line_search_loop_
  bool is_continue_line_search_loop() const {
    return is_continue_line_search_loop_;
  }
  void set_is_continue_line_search_loop(bool is_continue_line_search_loop) {
    is_continue_line_search_loop_ = is_continue_line_search_loop;
  }
  // line_search_potential_energy_
  double line_search_potential_energy() const {
    return line_search_potential_energy_;
  }
  void set_line_search_potential_energy(double line_search_potential_energy) {
    line_search_potential_energy_ = line_search_potential_energy;
  }
  // max_direction_
  double max_direction() const { return max_direction_; }
  void set_max_direction(double max_direction) {
    max_direction_ = max_direction;
  }
  // old_potential_energy_
  double old_potential_energy() const { return old_potential_energy_; }
  void set_old_potential_energy(double old_potential_energy) {
    old_potential_energy_ = old_potential_energy;
  }
  // potential_energy_difference_
  double potential_energy_difference() const {
    return potential_energy_difference_;
  }
  void set_potential_energy_difference(double potential_energy_difference) {
    potential_energy_difference_ = potential_energy_difference;
  }

  bool is_continue_line_search_loop_;     // whether to continue line search
  double alpha_;                          // alpha used in line search algorithm
  // beta used in calculation of direction
  double beta_;
  // energy value to check in deciding to terminate minimization
  double energy_check_;
  double force_dot_direction_;            // force multiplied by direction
  double line_search_potential_energy_;   // potential energy of the line search
  double max_direction_;                  // max of direction
  // max allowed displacement for an atom for an iteration
  double max_displacement_;
  // potential energy from previous iteration
  double old_potential_energy_;
  // difference in potential energy from previous iteration to current iteration
  double potential_energy_difference_;
  double tolerance_;                      // tolerance
  int max_iteration_;                     // max number of iterations
  Matrix<double> direction_;              // direction
};

#endif  // MMM_V14_6_CONJUGATE_GRADIENT_H_

