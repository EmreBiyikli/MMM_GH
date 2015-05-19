// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// model.h
// *****************************************************************************
// Consists of Model class that handles the model. In addition to the
// construction of the class, it needs to be initialized by either the
// InitFullAtomistic function or the InitMMM function, depending on the choice
// of model. The actual model is built by the BuildModel function. Besides to
// the model, types of a list of atoms can explicitly be set.

#ifndef MMM_V14_6_Model_H_
#define MMM_V14_6_Model_H_

#include <string>
#include <vector>

#include "matrix.h"
#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Model : public Mediator {
 public:
  // maximum number of element an atom can belong to
  const int kMaxAtomElementNum = 100;
  const int kMaxElementAtomNum = 1000;  // maximum number of atoms in an element
  
  explicit Model() : Mediator() {}
  ~Model() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Model(Model&) = delete;
  // Model& operator=(const Model&) = delete;

  // Build the model.
  void BuildModel();
  // Initializes the class to the full_atomisitc model.
  void InitFullAtomistic();
  // Initializes the class to the MMM model by initializing the scheme.
  void InitMMM(string scheme);
  // Sets input atom_list to the input_type.
  void SetAtomListType(const vector<int>& atom_list, int input_type);

  // Accessor and mutator functions:
  // is_full_atomistic_
  bool is_full_atomistic() const { return is_full_atomistic_; }
  void set_is_full_atomistic(bool is_full_atomistic) {
    is_full_atomistic_ = is_full_atomistic;
  }
  // scheme_
  string scheme() const { return scheme_; }
  void set_scheme(string scheme) { scheme_ = scheme; }
  // style_
  string style() const { return style_; }
  void set_style(string style) { style_ = style; }

  Matrix<bool> is_rep_;               // whether the atoms are rep atoms
  Matrix<double> atom_mass_;          // masses of atoms
  Matrix<double> force_weight_;       // force weights of atoms
  Matrix<double> shape_function_;     // shape functions of atoms
  Matrix<int> atom_element_;          // elements of atoms they belong to
  // number of elements of atoms they belong to
  Matrix<int> atom_element_num_;
  // types of atoms
  // 1: irep, 2: nirep, 3: psmp, 4: ssmp, 5: nsmp
  Matrix<int> atom_type_;
  Matrix<int> atom_type_num_;         // atom numbers of each atom type
  Matrix<int> element_atom_;          // atoms of elements
  Matrix<int> element_atom_num_;      // number of atoms of elements
  // center atoms of elements (-1 if there is no ghost atom in the element)
  Matrix<int> element_center_atom_;

 private:
  // Builds elements of atoms.
  void BuildAtomElement();
  // Builds line elements of atoms.
  void BuildAtomElementLine();
  // Builds tetrahedra elements of atoms.
  void BuildAtomElementTetrahedron();
  // Builds triangle elements of atoms.
  void BuildAtomElementTriangle();
  // Builds atoms of elements.
  void BuildElementAtom();
  // Builds center atoms of elements.
  void BuildElementCenterAtom();
  // Builds force weights.
  void BuildForceWeight();
  // Builds masses.
  void BuildMass();
  // Builds scheme.
  void BuildScheme();
  // Builds shape functions.
  void BuildShapeFunction();
  // Counts atom number of each type.
  void Model::CountAtomTypeNum();

  bool is_full_atomistic_;  // whether model is full atomistic
  // scheme
  // no_ssmp: no secondary sampling atoms
  // all_ssmp: all secondary sampling atoms
  // ssmp_around_irep: secondary sampling atoms in potential cut-off radius of
  //                   interpolating rep atoms
  // ssmp_around_psmp: secondary sampling atoms in potential cut-off radius of
  //                   primary sampling atoms
  string scheme_;
  string style_;            // style: full_atomistic or MMM
};

#endif  // MMM_V14_6_Model_H_

