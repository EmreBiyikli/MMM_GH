// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// mmm.h
// *****************************************************************************
// Consists of MMM class that handles the simulation. In addition to the
// construction of the class, it needs to be initialized by the Init function.
// Simulation is performed by the Simulate function. The class needs to be
// finalized by the Final function.

#ifndef MMM_V14_6_MMM_H_
#define MMM_V14_6_MMM_H_

const double kPi = 3.14159265;  // PI
// Boltzmann constant in units of (g/mole)(A^2/10fs^2)(1/K)
const double kBoltzmannConst = 0.83144286*1e-4;

#include <string>
#include <vector>

using std::string;
using std::vector;

// See comment at top of file for a complete description.
class MMM {
 public:  // public 1
  MMM();
  ~MMM();

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // MMM(MMM&) = delete;
  // MMM& operator=(const MMM&) = delete;

  // Checks for lost atoms and terminates the simulation if there is any.
  void CheckLostAtom();
  // Finalizes the class.
  void Final();
  // Initializes the class by initializing the simulation name.
  void Init(string simulation_name);
  // Performs simulation.
  void Simulate();

  // Accessor and mutator functions:
  // current_iteration_
  int current_iteration() const { return current_iteration_; }
  void set_current_iteration(int current_iteration) {
    current_iteration_ = current_iteration;
  }
  void increase_current_iteration() { current_iteration_++; }
  // dimension_
  int dimension() const { return dimension_; }
  void set_dimension(int dimension) { dimension_ = dimension; }
  // domain_boundary_
  vector<double> domain_boundary() const { return domain_boundary_; }
  double domain_boundary_at(int dimension) const {
    return domain_boundary_[dimension];
  }
  void set_domain_boundary(vector<double> domain_boundary) {
    domain_boundary_ = domain_boundary;
  }
  // is_dynamic_
  bool is_dynamic() const { return is_dynamic_; }
  void set_is_dynamic(bool is_dynamic) { is_dynamic_ = is_dynamic; }
  // is_iteration_loop_continue_
  bool is_iteration_loop_continue() const {
    return is_iteration_loop_continue_;
  }
  void set_is_iteration_loop_continue(bool is_iteration_loop_continue) {
    is_iteration_loop_continue_ = is_iteration_loop_continue;
  }
  // is_MPI_
  bool is_MPI() const { return is_MPI_; }
  void set_is_MPI(bool is_MPI) { is_MPI_ = is_MPI; }
  // kinetic_energy_
  double kinetic_energy() const { return kinetic_energy_; }
  void set_kinetic_energy(double kinetic_energy) {
    kinetic_energy_ = kinetic_energy;
  }
  void add_kinetic_energy(double kinetic_energy) {
    kinetic_energy_ += kinetic_energy;
  }
  // mass_
  double mass() const { return mass_; }
  void set_mass(double mass) { mass_ = mass; }
  // name_
  string name() const { return name_; }
  void set_name(string name) { name_ = name; }
  // potential_energy_
  double potential_energy() const { return potential_energy_; }
  void set_potential_energy(double potential_energy) {
    potential_energy_ = potential_energy;
  }
  void add_potential_energy(double potential_energy) {
    potential_energy_ += potential_energy;
  }
  // processing_element_
  int* processing_element_address() { return &processing_element_; }
  int processing_element() const { return processing_element_; }
  void set_processing_element(int processing_element) {
    processing_element_ = processing_element;
  }
  // processing_element_num_
  int* processing_element_num_address() { return &processing_element_num_; }
  int processing_element_num() const { return processing_element_num_; }
  void set_processing_element_num(int processing_element_num) {
    processing_element_num_ = processing_element_num;
  }
  // solver_
  string solver() const { return solver_; }
  void set_solver(string solver) { solver_ = solver; }
  // total_energy_
  double total_energy() const { return total_energy_; }
  void set_total_energy(double total_energy) { total_energy_ = total_energy; }

 private:
  // Executes dimension command.
  void MMM::CommandDimension();
  // Executes initial configuration command.
  void CommandInitialConfiguration();
  // Executes load command.
  void CommandLoad();
  // Executes mesh command.
  void CommandMesh();
  // Executes model command.
  void CommandModel();
  // Executes neighbor command.
  void CommandNeighbor();
  // Executes output command.
  void CommandOutput();
  // Executes potential command.
  void CommandPotential();
  // Executes run command.
  void CommandRun();
  // Executes select command.
  void CommandSelect();
  // Executes temperature command.
  void CommandTemperature();
  // Executes type command.
  void CommandType();

  bool is_dynamic_;                   // whether the simulation is dynamic
  bool is_iteration_loop_continue_;   // whether the iteration continues
  bool is_MPI_;                       // whether MPI is on
  double kinetic_energy_;             // total kinetic energy
  double potential_energy_;           // total potential energy
  double total_energy_;               // total energy
  double mass_;                       // mass of an atom
  int current_iteration_;             // current iteration
  int dimension_;                     // dimension
  int processing_element_;            // current processing element
  int processing_element_num_;        // number of processing elements
  string name_;                       // simulation name
  // solver: conjugate gradient or velocity verlet
  string solver_;
  vector<double> domain_boundary_;    // boundaries of the positions of atoms

  class AddIn* add_in_;                           // add in
  class AtomGroup* atom_group_;                   // atom group
  class ConjugateGradient* conjugate_gradient_;   // conjugate gradient
  class Input* input_;                            // input
  class Mesh* mesh_;                              // mesh
  class MMM* mmm_;                                // mmm
  class Model* model_;                            // model
  class Neighbor* neighbor_;                      // neighbor
  class Output* output_;                          // output
  class Potential* potential_;                    // potential
  class Select* select_;                          // select
  class Temperature* temperature_;                // temperature
  class Time* time_;                              // time
  class VelocityVerlet* velocity_verlet_;         // velocity verlet

 public:  // public 2
  // Accessor and mutator functions for sub-classes. See MANUAL.txt for details.
  // add_in_
  AddIn* add_in() { return add_in_; }
  // atom_group_
  AtomGroup* atom_group() { return atom_group_; }
  // conjugate_gradient_
  ConjugateGradient* conjugate_gradient() { return conjugate_gradient_; }
  // input_
  Input* input() { return input_; }
  // mesh_
  Mesh* mesh() { return mesh_; }
  // model_
  Model* model() { return model_; }
  // neighbor_
  Neighbor* neighbor() { return neighbor_; }
  // output_
  Output* output() { return output_; }
  // potential_
  Potential* potential() { return potential_; }
  // select
  Select* select() { return select_; }
  // temperature_
  Temperature* temperature() { return temperature_; }
  // time_
  Time* time() { return time_; }
  // velocity_verlet_
  VelocityVerlet* velocity_verlet() { return velocity_verlet_; }
};

#endif  // MMM_V14_6_MMM_H_

