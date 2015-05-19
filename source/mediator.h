// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// mediator.h
// *****************************************************************************
// Consists of Mediator class that render a group of other classes visible to
// each other. All other classes inherit from Mediator class. See MANUAL.txt for
// details.

#ifndef MMM_V14_6_MEDIATOR_H_
#define MMM_V14_6_MEDIATOR_H_

#include "mmm.h"

// See comment at top of file for a complete description.
class Mediator {
 public:
  // Inline Constructor that initializes member class pointers to the member
  // classes of the input MMM class.
  Mediator() {}
  ~Mediator() {}
  void InitPointer(MMM* mmm) {
    add_in_ = mmm->add_in();
    atom_group_ = mmm->atom_group();
    conjugate_gradient_ = mmm->conjugate_gradient();
    input_ = mmm->input();
    mesh_ = mmm->mesh();
    mmm_ = mmm;
    model_ = mmm->model();
    neighbor_ = mmm->neighbor();
    output_ = mmm->output();
    potential_ = mmm->potential();
    temperature_ = mmm->temperature();
    time_ = mmm->time();
    velocity_verlet_ = mmm->velocity_verlet();
  }

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Mediator(Mediator&) = delete;
  // Mediator& operator=(const Mediator&) = delete;

  // Accessor and mutator functions:
  // add_in_
  AddIn* add_in() { return add_in_; }
  void set_add_in(AddIn* add_in) { add_in_ = add_in; }
  // atom_group_
  AtomGroup* atom_group() { return atom_group_; }
  void set_atom_group(AtomGroup* atom_group) { atom_group_ = atom_group; }
  // conjugate_gradient_
  ConjugateGradient* conjugate_gradient() { return conjugate_gradient_; }
  void set_conjugate_gradient(ConjugateGradient* conjugate_gradient) {
    conjugate_gradient_ = conjugate_gradient;
  }
  // input_
  Input* input() { return input_; }
  void set_input(Input* input) { input_ = input; }
  // mesh_
  Mesh* mesh() { return mesh_; }
  void set_mesh(Mesh* mesh) { mesh_ = mesh; }
  // mmm_
  MMM* mmm() { return mmm_; }
  void set_mmm(MMM* mmm) { mmm_ = mmm; }
  // model_
  Model* model() { return model_; }
  void set_model(Model* model) { model_ = model; }
  // neighbor_
  Neighbor* neighbor() { return neighbor_; }
  void set_neighbor(Neighbor* neighbor) { neighbor_ = neighbor; }
  // output_
  Output* output() { return output_; }
  void set_output(Output* output) { output_ = output; }
  // potential_
  Potential* potential() { return potential_; }
  void set_pair_potential(Potential* potential) { potential_ = potential; }
  // temperature_
  Temperature* temperature() { return temperature_; }
  void set_temperature(Temperature* temperature) { temperature_ = temperature; }
  // time_
  Time* time() { return time_; }
  void set_time(Time* time) { time_ = time; }
  // velocity_verlet_
  VelocityVerlet* velocity_verlet() { return velocity_verlet_; }
  void set_velocity_verlet(VelocityVerlet* velocity_verlet) {
    velocity_verlet_ = velocity_verlet;
  }

 private:
  AddIn* add_in_;                           // add in
  AtomGroup* atom_group_;                   // atom group
  ConjugateGradient* conjugate_gradient_;   // conjugate gradient
  Input* input_;                            // input
  Mesh* mesh_;                              // mesh
  MMM* mmm_;                                // mmm
  Model* model_;                            // model
  Neighbor* neighbor_;                      // neighbor
  Output* output_;                          // output
  Potential* potential_;                    // potential
  Temperature* temperature_;                // temperature
  Time* time_;                              // time
  VelocityVerlet* velocity_verlet_;         // velocity verlet
};

#endif  // MMM_V14_6_MEDIATOR_H_

