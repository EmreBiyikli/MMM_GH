// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// neighbor.cc
//******************************************************************************

#include "neighbor.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "kdtree.h"

using std::vector;

// Independent functions (by means of unnamed namespace)
// *****************************************************************************

namespace {

// KdTree
//******************************************************************************
// (independent)
// Builds neighbor lists by a KDTree search algorithm.
// Note that this function is external so that it is not styled.
void KdTree(const vector<double>& position, double cutoff_radius,
            int max_neighbor_num, vector<int>* neighbor_list,
            vector<int>* neighbor_num) {
  // init
  int atom_num = position.size() / 3;
  neighbor_list->clear();
  neighbor_num->clear();
  int *id;
  kdres *set;
  kdtree *kd;
  kd = kd_create(3);
  int *ids = new int[atom_num]();
  // insert points
  for (int i = 0; i < atom_num; i++) {
    ids[i] = i;
    kd_insert3(kd, position[3 * i + 0], position[3 * i + 1],
               position[3 * i + 2], &ids[i]);
  }
  // build neighbor
  for (int i = 0; i < atom_num; i++) {
    // set
    set = kd_nearest_range3(kd, position[3 * i + 0], position[3 * i + 1],
                            position[3 * i + 2], cutoff_radius);
    // set size
    int setSize = kd_res_size(set);
    // set size -> neighbor_num[i]
    neighbor_num->push_back(setSize - 1); // -1 is for the atom itself
    if ((*neighbor_num)[i] > max_neighbor_num) {
      printf("nearest neighbor num %d is bigger than the allowed maximum %s\n",
             (*neighbor_num)[i], max_neighbor_num);
      exit(1);
    }
    // set -> neighbor_list
    int count = 0;
    for (int j = 0; j < setSize; j++) {
      id = (int*)kd_res_item_data(set);
      kd_res_next(set);
      if (*id == i) continue; // eliminate the atom itself
      (*neighbor_list).push_back(*id);
      count++;
    }
    neighbor_list->insert(neighbor_list->end(),
                          max_neighbor_num - (*neighbor_num)[i], -1);
    // clear set
    kd_res_free(set);
  }
  // clear tree
  kd_free(kd);
}

}  // namespace

// Public
// *****************************************************************************

// InitAuto
//******************************************************************************
void Neighbor::InitAuto(double neighbor_cutoff_radius,
                        double potential_cutoff_radius) {
  set_cutoff_radius(neighbor_cutoff_radius);
  set_critical_displacement_square(pow((cutoff_radius() -
       potential_cutoff_radius) / 2.0, 2));
  set_update_style("auto");
  set_cutoff_radius_square(cutoff_radius() * cutoff_radius());
  set_build_count(0);
  set_is_first_update(true);
}

// InitEvery
//******************************************************************************
void Neighbor::InitEvery(double neighbor_cutoff_radius,
                         double potential_cutoff_radius,
                         int update_frequency) {
  set_cutoff_radius(neighbor_cutoff_radius);
  set_critical_displacement_square(pow((cutoff_radius() -
       potential_cutoff_radius) / 2.0, 2));
  set_update_frequency(update_frequency);
  set_update_style("every");
  set_cutoff_radius_square(cutoff_radius() * cutoff_radius());
  set_build_count(0);
  set_is_first_update(true);
}

// BuildNeighbor
//******************************************************************************
void Neighbor::BuildNeighbor(const vector<double>& position, int iteration) {
  set_atom_num(position.size() / 3);
  if (is_first_update() && update_style().compare("auto") == 0) {
    set_memorized_position(position);
  }
  bool is_update = false;
  if (iteration == 0) {
    is_update = true;  // First time build
  } else if (update_style().compare("auto") == 0) {  // Auto build
    for (int i = 0; i < atom_num(); i++) {
      double displacement_square = 
        pow(position[3 * i + 0] - memorized_position(i, 0), 2) +
        pow(position[3 * i + 1] - memorized_position(i, 1), 2) +
        pow(position[3 * i + 2] - memorized_position(i, 2), 2);
      if (displacement_square > critical_displacement_square()) {
        is_update = true;
        break;
      }
    }
  } else if (update_style().compare("every") == 0) {  // Every build
    if (iteration % update_frequency() == 0) is_update = true;
  }
  // Build
  if (is_update) {
    increase_build_count();
    vector<int> temp_neighbor_list;
    vector<int> temp_neighbor_num;
    KdTree(position, cutoff_radius(), kMaxNeighborNum, &temp_neighbor_list,
           &temp_neighbor_num);
    set_neighbor_list(temp_neighbor_list);
    set_neighbor_num(temp_neighbor_num);
    if (update_style().compare("auto") == 0) set_memorized_position(position);
  }
  if (is_first_update()) set_is_first_update(false);
}

