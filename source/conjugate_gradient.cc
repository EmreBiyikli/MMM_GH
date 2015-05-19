// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// conjugate_gradient.cc
// *****************************************************************************

#include "conjugate_gradient.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "add_in.h"
#include "atom_group.h"
#include "matrix.h"
#include "mesh.h"
#include "mmm.h"
#include "model.h"
#include "neighbor.h"
#include "output.h"
#include "potential.h"

using std::abs;
using std::min;
using std::max;
using std::vector;

// Public
// *****************************************************************************

// Init
// *****************************************************************************
void ConjugateGradient::Init(double max_displacement, double tolerance,
                             int max_iteration) {
  set_max_displacement(max_displacement);
  set_tolerance(tolerance);
  set_max_iteration(max_iteration);
  direction_.set_size(atom_group()->atom_num(), 3);
}

// Run
// *****************************************************************************
void ConjugateGradient::Run() {
  output()->OutputToLog("");
  output()->OutputToLog("conjugate gradient run");
  output()->OutputToLog("dashed line");
  output()->OutputToLog("iteration potential_energy energy_check");
  output()->OutputToLog("dashed line");
  set_energy_check(0.0);
  mmm()->set_current_iteration(-1);
  mmm()->set_is_iteration_loop_continue(true);  
  // Loop  
  while (mmm()->is_iteration_loop_continue()) {
    mmm()->increase_current_iteration();    
    neighbor()->BuildNeighbor(atom_group()->position_.get(),
                              mmm()->current_iteration());  
    if (mmm()->current_iteration() == 0) {
      // Init
      potential()->Compute(true, true);
      for (int i = 0; i < atom_group()->atom_num(); i++) {
        if (!model()->is_rep_.get_1d_element(i)) continue;
        for (int j = 0; j < mmm()->dimension(); j++) {
          atom_group()->old_applied_force_.set_element(i, j,
              atom_group()->applied_force_.get_element(i, j));
          direction_.set_element(i, j, 
              atom_group()->applied_force_.get_element(i, j));
        }
      }
    } else {
      // Iterate
      Iterate();
    }
    if (mmm()->current_iteration() % output()->screen_update_frequency() == 0 ||
        !mmm()->is_iteration_loop_continue()) {
       output()->OutputToLog(output()->ToString(mmm()->current_iteration()) +
                                " " +
                                output()->ToString(mmm()->potential_energy()) +
                                " " +
                                output()->ToString(energy_check()));
    }
    output()->OutputAllToFile();
    mmm()->CheckLostAtom();
    // Adds displacement in the static paper crack example.
    add_in()->StaticPaperCrackDisplacement();
  }
}

// Private
// *****************************************************************************

// ApplyDisplacementLoad
// *****************************************************************************
void ConjugateGradient::ApplyDisplacementLoad() {
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

// CheckLineSearchStopCriteria
// *****************************************************************************
// Terminates if either the potential energy difference is smaller than or equal
// to an ideal potential energy difference or alpha is too small.
void ConjugateGradient::CheckLineSearchStopCriteria() {
  double ideal_potential_energy_difference = -0.4 * alpha() *
      force_dot_direction();
  set_potential_energy_difference(mmm()->potential_energy() -
                                      line_search_potential_energy());
  if (potential_energy_difference() <= ideal_potential_energy_difference ||
      alpha() < 0.0001) {
    set_is_continue_line_search_loop(false);
  }
}

// CheckStopCriteria
// *****************************************************************************
// Terminates if normalized potential energy difference is smaller than the
// tolerance or maximum allowed number of iterations is reached.
void ConjugateGradient::CheckStopCriteria() {
  double energy_norm = 0;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    energy_norm += atom_group()->energy_.get_1d_element(i) *
        atom_group()->energy_.get_1d_element(i);
  }
  energy_norm = sqrt(energy_norm);
  set_energy_check(abs(potential_energy_difference() / energy_norm));
  if (mmm()->current_iteration() > 0 &&
      (energy_check() <= tolerance() ||
       mmm()->current_iteration() == max_iteration())) {
    mmm()->set_is_iteration_loop_continue(false);
  }
}

// ComputeBeta
// *****************************************************************************
void ConjugateGradient::ComputeBeta() {
  double temp_1 = 0, temp_2 = 0, temp_3 = 0;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      temp_1 += atom_group()->applied_force_.get_element(i, j) *
                    atom_group()->applied_force_.get_element(i, j);
      temp_2 += atom_group()->applied_force_.get_element(i, j) *
                    atom_group()->old_applied_force_.get_element(i, j);
      temp_3 += atom_group()->old_applied_force_.get_element(i, j) *
                    atom_group()->old_applied_force_.get_element(i, j);      
    }
  }
  set_beta(max((temp_1 - temp_2) / temp_3, 0.0));
}

// ComputeDirection
// *****************************************************************************
void ConjugateGradient::ComputeDirection() {
  set_force_dot_direction(0.0);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      direction_.set_element(i, j,
                             atom_group()->applied_force_.get_element(i, j) +
                                 beta() * direction_.get_element(i, j));
      add_force_dot_direction(atom_group()->applied_force_.get_element(i, j) *
                                  direction_.get_element(i, j));
    }
  }
  if (force_dot_direction() <= 0) {
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      if (!model()->is_rep_.get_1d_element(i)) continue;
      for (int j = 0; j < mmm()->dimension(); j++) {
        direction_.set_element(i, j, atom_group()->applied_force_.get_element(i, j));
      }
    }
  }
  set_max_direction(-1.0e6);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      if (direction_.get_element(i, j) > max_direction()) {
        set_max_direction(direction_.get_element(i, j));
      }
    }
  }
}

// DisplaceGhostAtom
// *****************************************************************************
void ConjugateGradient::DisplaceGhostAtom() {
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
        atom_group()->displacement_.add_element(i, k,
            interpolation_weight *
                atom_group()->displacement_.get_element(node, k));
      }
    }
  }
}

// DisplaceRepAtom
// *****************************************************************************
void ConjugateGradient::DisplaceRepAtom() {
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (!model()->is_rep_.get_1d_element(i)) continue;
    for (int j = 0; j < mmm()->dimension(); j++) {
      atom_group()->displacement_.set_element(i, j,
          alpha() * direction_.get_element(i, j));
    }
  }
}

// Iterate
// *****************************************************************************
void ConjugateGradient::Iterate() {
  ComputeBeta();
  ComputeDirection();
  LineSearch();      
  atom_group()->old_applied_force_.set(atom_group()->applied_force_.get());
  potential()->Compute(false, true);
  CheckStopCriteria();
}

// LineSearch
// *****************************************************************************
// For details of the algorithm, check the book Numerical Optimization by
// Nocedal.
void ConjugateGradient::LineSearch() {
  set_is_continue_line_search_loop(true);
  set_alpha(min(max_displacement() / max_direction(), 1.0));
  set_line_search_potential_energy(mmm()->potential_energy());
  vector<double> line_search_position = atom_group()->position_.get();
  double reduce_ratio = 0.5;
  int step = 0;
  while (is_continue_line_search_loop()) {
    step++;
    if (step != 1) multiply_alpha(reduce_ratio);
    atom_group()->position_.set(line_search_position);
    DisplaceRepAtom();
    ApplyDisplacementLoad();
    if (!model()->is_full_atomistic()) DisplaceGhostAtom();    
    atom_group()->position_.add(atom_group()->displacement_.get());    
    potential()->Compute(true, false);   
    CheckLineSearchStopCriteria();
  }
  set_old_potential_energy(line_search_potential_energy());  
}

