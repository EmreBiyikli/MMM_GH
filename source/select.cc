// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// select.cc
// *****************************************************************************

#include "select.h"

#include <mpi.h>

#include <cmath>
#include <string>
#include <vector>

using std::string;
using std::vector;

// Public
// *****************************************************************************

// GetPosition
// *****************************************************************************
void Select::GetPosition(int ID, vector<double>* atom_list_position) {
  for (int i = 0; i < selection(ID).size(); i++) {
    for (int j = 0; j < 3; j++) {
      atom_list_position->push_back(position(selection_at(ID, i), j));
    }
  }
}

// Init
// *****************************************************************************
void Select::Init(const vector<double> &position) {
  set_position(position);
  set_atom_num(position.size() / 3);
  resize_selection(kMaxSelectionNum);
  for (int i = 0; i < atom_num(); i++) {
    push_selection(0, i);
  }
}

// SelectBlock
// *****************************************************************************
void Select::SelectBlock(int ID, vector<double> block) {
  for (int i = 0; i < atom_num(); i++) {
    bool is_inside = true;
    for (int j = 0; j < 3; j++) {
      if (position(i, j) < block[2 * j] || position(i, j) > block[2 * j + 1]) {
        is_inside = false;
      }
    }
    if (is_inside) push_selection(ID, i);
  }
  SelectSortAndUnique(ID);
}

// SelectFile
// *****************************************************************************
void Select::SelectFile(int ID, string input_file) {
  int processing_element, processing_element_num;
  MPI_Comm_rank(MPI_COMM_WORLD, &processing_element);
  MPI_Comm_size(MPI_COMM_WORLD, &processing_element_num);
  if (processing_element == 0) {
    input_file.insert(0, "input/");
    FILE* input_file_ptr = fopen(input_file.c_str(), "r");
    if (input_file_ptr == NULL) {
      printf("cannot open the selection input file: %s\n", input_file.c_str());
      exit(1);
    }
    int atom, atom_num;
    fscanf(input_file_ptr, "%d", &atom_num);
    for (int i = 0; i < atom_num; i++) {
      if (input_file_ptr == NULL) {
        printf("inconsistent atom num in file selection\n");
        exit(1);
      }
      fscanf(input_file_ptr, "%d", &atom);
      push_selection(ID, atom);
    }
    fclose(input_file_ptr);
    SelectSortAndUnique(ID);
  }
  if (processing_element_num > 1) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(selection_address(ID), selection_size(ID), MPI_INT, 0,
              MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// SelectGrid
// *****************************************************************************
// Boundary of the input selection is extended by 0.01 to make sure that the
// domain is covered.
void Select::SelectGrid(int ID, int in_ID, vector<double> grid) {
  vector<double> boundary = GetBoundary(in_ID);
  boundary[0] -= 0.01;
  boundary[1] += 0.01;
  boundary[2] -= 0.01;
  boundary[3] += 0.01;
  boundary[4] -= 0.01;
  boundary[5] += 0.01;
  // Grid point numbers
  vector<double> actual_grid(3, 0);
  vector<double> length(3, 0);
  vector<int> point_num(3, 1);
  int total_point_num = 1;
  for (int i = 0; i < 3; i++) {
    length[i] = boundary[2 * i + 1] - boundary[2 * i];
    point_num[i] = static_cast<int>(ceil(length[i] / grid[i]) + 1);
    actual_grid[i] = length[i] / static_cast<double>(point_num[i] - 1);
    total_point_num *= point_num[i];
  }
  // Grid point positions
  vector< vector<double> > point_position;
  vector<double> one_point_position;
  int count = 0;
  for (int i = 0; i < point_num[0]; i++) {
    double x_position = boundary[0] + i * actual_grid[0];
    for (int j = 0; j < point_num[1]; j++) {
      double y_position = boundary[2] + j * actual_grid[1];
      for (int k = 0; k < point_num[2]; k++) {
        double z_position = boundary[4] + k * actual_grid[2];
        one_point_position.clear();
        one_point_position.push_back(x_position);
        one_point_position.push_back(y_position);
        one_point_position.push_back(z_position);
        point_position.push_back(one_point_position);
        count++;
      }
    }
  }
  // Select atoms closest to grid points
  for (int i = 0; i < total_point_num; i++) {
    double minimum_distance_square = 1e6;
    int closest_atom;
    for (int j = 0; j < selection(in_ID).size(); j++) {
      double distance_square = 0;
      for (int k = 0; k < 3; k++) {
        double difference = position(selection_at(in_ID, j), k) -
            point_position[i][k];
        distance_square += difference*difference;
      }
      if (distance_square < minimum_distance_square) {
        minimum_distance_square = distance_square;
        closest_atom = selection_at(in_ID, j);
      }
    }
    push_selection(ID, closest_atom);
  }
  SelectSortAndUnique(ID);
}

// SelectID
// *****************************************************************************
void Select::SelectID(int ID, vector<string> expression) {
  for (int i = 0; i < expression.size(); i++) {
    int num;
    int first_colon_place = expression[i].find(":");
    int second_colon_place = expression[i].find(":", first_colon_place + 1);
    if (first_colon_place >= 0) {
      int end, incerement;
      string sub_string = expression[i].substr(0, first_colon_place);
      int start = atoi(sub_string.c_str());
      if (second_colon_place >= 0) {
        sub_string = expression[i].substr(first_colon_place + 1,
            second_colon_place - first_colon_place - 1);
        incerement = atoi(sub_string.c_str());
        sub_string = expression[i].substr(second_colon_place + 1);
        end = atoi(sub_string.c_str());
      } else {
        sub_string = expression[i].substr(first_colon_place + 1);
        end = atoi(sub_string.c_str());
        incerement = 1;
      }
      num = start;
      while (num <= end) {
        push_selection(ID, num);
        num += incerement;
      }
    } else {
      num = atoi(expression[i].c_str());
      push_selection(ID, num);
    }
  }
  SelectSortAndUnique(ID);
}

// SelectRadial
// *****************************************************************************
void Select::SelectRadial(double radius, int ID, int in_ID, int center_ID) {
  double radius_square = radius * radius;
  for (int i = 0; i < selection(center_ID).size(); i++) {
    int center_atom = selection_at(center_ID, i);
    for (int j = 0; j < selection(in_ID).size(); j++) {
      int in_atom = selection_at(in_ID, j);
      if (in_atom == center_atom) continue;
      double distance_square = 0;
      for (int k = 0; k < 3; k++) {
        double difference = position(in_atom, k) - position(center_atom, k);
        distance_square += difference * difference;
      }
      if (distance_square < radius_square) {
        push_selection(ID, selection_at(in_ID, j));
      }
    }
  }
  SelectSortAndUnique(ID);}

// SelectSubtract
// *****************************************************************************
void Select::SelectSubtract(int ID, int ID_1, int ID_2) {
  vector<int>::const_iterator iterator;
  for (int i = 0; i < selection(ID_1).size(); i++) {
    iterator = find(selection(ID_2).begin(), selection(ID_2).end(),
                    selection(ID_1)[i]);
    if (iterator == selection(ID_2).end()) {
      push_selection(ID, selection_at(ID_1, i));
    }
  }
  SelectSortAndUnique(ID);
}

// SelectSurface
// *****************************************************************************
void Select::SelectSurface(int ID, double depth, string surface) {
  // Surface flags
  vector<bool> surface_flag(6, false);
  if (surface.find("x-") != string::npos) surface_flag[0] = true;
  if (surface.find("x+") != string::npos) surface_flag[1] = true;
  if (surface.find("y-") != string::npos) surface_flag[2] = true;
  if (surface.find("y+") != string::npos) surface_flag[3] = true;
  if (surface.find("z-") != string::npos) surface_flag[4] = true;
  if (surface.find("z+") != string::npos) surface_flag[5] = true;
  // Boundary
  vector<double> boundary = GetBoundary(0);
  for (int i = 0; i < 6; i++) {
    if (surface_flag[i]) {
      if (i % 2 == 0) {
        boundary[i] += depth;
      } else {
        boundary[i] -= depth;
      }
    }
  }
  // Select
  for (int i = 0; i < atom_num(); i++) {
    bool is_selected = false;
    for (int j = 0; j < 3; j++) {
      if (surface_flag[2 * j] && position(i, j) <= boundary[2 * j]) {
        is_selected = true;
      }
      if (surface_flag[2 * j + 1] && position(i, j) >= boundary[2 * j + 1]) {
        is_selected = true;
      }
    }
    if (is_selected) push_selection(ID, i);
  }
  SelectSortAndUnique(ID);
}

// SelectUnite
// *****************************************************************************
void Select::SelectUnite(int ID, int ID_1, int ID_2) {
  insert_selection(ID, selection(ID_1));
  insert_selection(ID, selection(ID_2));
  SelectSortAndUnique(ID);
}

// SelectSortAndUnique
// *****************************************************************************
void Select::SelectSortAndUnique(int ID) {
  sort_selection(ID);
  unique_selection(ID);
}

// Private
// *****************************************************************************

// GetBoundary
// *****************************************************************************
vector<double> Select::GetBoundary(int ID) {  
  // TODO: activate with C++11 compiler
  // vector<double> boundary = {1e6, -1e6, 1e6, -1e6, 1e6, -1e6};
  vector<double> boundary(6, 1e6);
  boundary[1] = boundary[3] = boundary[5] = -1e6;
  for (int i = 0; i < selection(ID).size(); i++) {
    int atom = selection_at(ID, i);
    for (int j = 0; j < 3; j++) {
      if (position(atom, j) < boundary[2 * j]) {
        boundary[2 * j] = position(atom, j);
      }
      if (position(atom, j) > boundary[2 * j + 1]) {
        boundary[2 * j + 1] = position(atom, j);
      }
    }
  }
  return boundary;
}

