// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// atom_group.cc
// *****************************************************************************

#include "atom_group.h"

#include <vector>

using std::vector;

// Public
// *****************************************************************************

// GetBoundary
// *****************************************************************************
vector<double> AtomGroup::GetBoundary() const {
  // TODO: activate with C++11 compiler
  // vector<double> boundary = {1e6, -1e6, 1e6, -1e6, 1e6, -1e6};
  vector<double> boundary(6, 1e6);
  boundary[1] = boundary[3] = boundary[5] = -1e6;
  for (int i = 0; i < atom_num(); i++) {
    for (int j = 0; j < 3; j++) {
      if (position_.get_element(i, j) < boundary[2 * j]) {
        boundary[2 * j] = position_.get_element(i, j);
      }
      if (position_.get_element(i, j) > boundary[2 * j + 1]) {
        boundary[2 * j + 1] = position_.get_element(i, j);
      }
    }
  }
  return boundary;
}

// Init
// *************************************
void AtomGroup::Init(int atom_num) {
  applied_force_.set_size(atom_num, 3);
  atom_num_ = atom_num;
  degree_of_freedom_num_ = 3 * atom_num_;
  displacement_.set_size(atom_num, 3);
  energy_.set_1d_size(atom_num);
  force_.set_size(atom_num, 3);
  initial_position_.set_size(atom_num, 3);
  is_load_.set_1d_size(atom_num);
  load_type_.set_size(atom_num, 3);
  load_value_.set_size(atom_num, 3);
  old_applied_force_.set_size(atom_num, 3);
  position_.set_size(atom_num, 3);
  processing_element_.set_1d_size(atom_num);
  velocity_.set_size(atom_num, 3);
}

