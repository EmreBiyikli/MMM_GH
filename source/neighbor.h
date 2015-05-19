// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// neighbor.h
// ****************************************************************************
// (independent)
// Consists of Neighbor class that handles the neighbors. In addition to the
// construction of the class, it needs to be initialized by either the InitAuto
// function or the InitEvery function, depending on the choice of update style.
// Neighbor lists are built by the BuildNeighbor function. Lists are built by a
// KDTree search algorithm.

#ifndef MMM_V14_6_NEIGHBOR_H_
#define MMM_V14_6_NEIGHBOR_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

// See comment at top of file for a complete description.
class Neighbor{
 public:
  const int kMaxNeighborNum = 1000;  // maximum number of neighbors of an atom

  Neighbor() {}
  ~Neighbor() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Neighbor(Neighbor&) = delete;
  // Neighbor& operator=(const Neighbor&) = delete;

  // Builds neighbor lists with input position, if it is time to update.
  // Iteration is used to check for update.
  void BuildNeighbor(const vector<double>& position, int iteration);
  // Initializes the class to the Auto update style by initializing the input
  // member variables.
  void InitAuto(double neighbor_cutoff_radius,
                double potential_cutoff_radius);
  // Initializes the class to the Every update style by initializing the input
  // member variables.
  void InitEvery(double neighbor_cutoff_radius,
                 double potential_cutoff_radius, int update_frequency);

  // Accessor and mutator functions:
  // atom_num_
  int atom_num() const { return atom_num_; }
  void set_atom_num(int atom_num) { atom_num_ = atom_num; }
  // build_count_
  int build_count() const { return build_count_; }
  void set_build_count(int build_count) { build_count_ = build_count; }
  void increase_build_count() { build_count_++; }
  // critical_displacement_square_
  double critical_displacement_square() const {
    return critical_displacement_square_;
  }
  void set_critical_displacement_square(double critical_displacement_square) {
    critical_displacement_square_ = critical_displacement_square;
  }
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
  // is_first_update_
  bool is_first_update() const { return is_first_update_; }
  void set_is_first_update(bool is_first_update) {
    is_first_update_ = is_first_update;
  }
  // memorized_position_
  double memorized_position(int atom, int dimension) const {
    return memorized_position_[3 * atom + dimension];
  }
  void set_memorized_position(const vector<double> &memorized_position) {
    memorized_position_ = memorized_position;
  }
  // negihbor_list_
  const vector<int>& neighbor_list() const { return neighbor_list_; }
  int neighbor_list_at(int atom, int neighbor) const {
    return neighbor_list_[kMaxNeighborNum * atom + neighbor];
  }
  void set_neighbor_list(const vector<int>& neighbor_list) {
    neighbor_list_ = neighbor_list;
  }
  // neighbor_num_
  int neighbor_num(int atom) const { return neighbor_num_[atom]; }
  void set_neighbor_num(const vector<int>& neighbor_num) {
    neighbor_num_ = neighbor_num;
  }
  // update_frequency_
  int update_frequency() const { return update_frequency_; }
  void set_update_frequency(int update_frequency) {
    update_frequency_ = update_frequency;
  }
  // update_style_
  string update_style() const { return update_style_; }
  void set_update_style(string update_style) { update_style_ = update_style; }

 private:
  bool is_first_update_;                  // whether it is the first update
  double critical_displacement_square_;   // square of the critical displacement
  double cutoff_radius_;                // cutoff distance
  double cutoff_radius_square_;         // square of the cutoff distance
  int atom_num_;                          // atom number
  int build_count_;                       // count of builds
  int update_frequency_;                  // frequency to updates
  string update_style_;                   // update style
  vector<double> memorized_position_;     // positions saved in memory
  vector<int> neighbor_list_;             // neighbor lists of atoms
  vector<int> neighbor_num_;              // neighbor number of atoms
};

#endif  // MMM_V14_6_NEIGHBOR_H_

