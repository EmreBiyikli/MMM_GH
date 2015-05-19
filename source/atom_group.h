// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// atom_group.h
// *****************************************************************************
// Consists of AtomGroup class that stores information about a group of atoms.
// In addition to the construction of the class, it needs to be initialized
// by the Init function. Information is stored in Matrix type (see matrix.h).

#ifndef MMM_V14_6_ATOM_GROUP_H_
#define MMM_V14_6_ATOM_GROUP_H_

#include <vector>

#include "matrix.h"

// See comment at top of file for a complete description.
class AtomGroup {
 public:
  AtomGroup() {}
  ~AtomGroup() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // AtomGroup(AtomGroup&) = delete;
  // AtomGroup& operator=(const AtomGroup&) = delete;

  // Returns boundaries of the positions of the atom group.
  vector<double> GetBoundary() const;
  // Initializes the class by the input atom number by means of resizing the
  // member variables.
  void Init(int atom_num);

  // Accessor and mutator functions:
  // atom_num_
  int atom_num() const { return atom_num_; }
  void set_atom_num(int atom_num) { atom_num_ = atom_num; }
  // degree_of_freedom_num_
  int degree_of_freedom_num() const { return degree_of_freedom_num_; }
  void set_degree_of_freedom_num(int degree_of_freedom_num) {
    degree_of_freedom_num_ = degree_of_freedom_num;
  }

  Matrix<bool> is_load_;              // whether there is a loading on atoms
  Matrix<double> displacement_;       // displacement of atoms
  Matrix<double> energy_;             // energy of atoms
  // applied force on atoms (force loading is imposed)
  Matrix<double> applied_force_;
  // computed force on atoms (force loading is not imposed)
  Matrix<double> force_;
  Matrix<double> initial_position_;   // initial positions of atoms
  Matrix<double> load_value_;         // value of loading, if there is any
  // applied force on atoms from previous iteration
  Matrix<double> old_applied_force_;
  Matrix<double> position_;           // position of atoms
  Matrix<double> velocity_;           // velocity of atoms
  // types of loading
  // 1: displacement set, 2: displacement add, 3: force set, 4: force add
  Matrix<int> load_type_;
  // which processing element atoms belong to
  Matrix<int> processing_element_;

 private:
  int atom_num_;                      // atom number of atom group
  // degree of freedom number of atom group
  int degree_of_freedom_num_;
};

#endif  // MMM_V14_6_ATOM_GROUP_H_

