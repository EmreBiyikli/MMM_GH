// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// build_initial_configuration.cc
// *****************************************************************************

#include "build_initial_configuration.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

// Build1DUniform
// *****************************************************************************
// (independent)
// Starts from left most position (taken to be equal to 0) and advances to the
// right.
void Build1DUniform(const vector<double>& input_domain_size, 
                    double initial_spacing, vector<double>* initial_position) {
  double x_increment = initial_spacing;
  double x_position = 0;
  while (x_position < input_domain_size[0]) {
    (*initial_position).push_back(x_position);
    (*initial_position).push_back(0);
    (*initial_position).push_back(0);
    x_position += x_increment;
  }
}

// Build2DHexagonal
// *****************************************************************************
// (independent)
// Starts from left-bottom position (taken to be equal to (0, 0)) and advances
// to the top in the outer loop and to the right in the inner loop.
void Build2DHexagonal(const vector<double>& input_domain_size,
                      double initial_spacing,
                      vector<double>* initial_position) {
  double x_increment = initial_spacing;
  double y_increment = (sqrt(3.0) / 2) * initial_spacing;
  int y_count = 0;
  double y_position = 0.0;
  while (y_position < input_domain_size[1]) {
    y_count++;
    int x_count = 0;
    double x_position = 0.0;
    while (x_position < input_domain_size[0]) {
      x_count++;
      if (x_count == 1 && y_count % 2 == 0) x_position += x_increment / 2;
      (*initial_position).push_back(x_position);
      (*initial_position).push_back(y_position);
      (*initial_position).push_back(0);
      x_position += x_increment;
    }
    y_position += y_increment;
  }
}

// Build3DFCC
// *****************************************************************************
// (independent)
// Starts from left-bottom-front position (taken to be equal to (0, 0, 0)) and
// advances to the back in the outer loop, to the top in the middle loop, and to
// the right in the inner loop.
void Build3DFCC(const vector<double>& input_domain_size, double initial_spacing,
                vector<double>* initial_position) {
  double x_increment = sqrt(2.0) * initial_spacing;
  double y_increment = sqrt(2.0) / 2.0 * initial_spacing;
  double z_increment = sqrt(2.0) / 2.0 * initial_spacing;
  int z_count = 0;
  double z_position = 0.0;
  while (true) {
    z_count++;
    z_position = (z_count - 1) * z_increment;
    if (z_position > input_domain_size[2]) break;
    int y_count = 0;
    while (true) {
      y_count++;
      double y_position = (y_count - 1) * y_increment;
      if (y_position > input_domain_size[1]) break;
      int x_count = 0;
      while (true) {
        x_count++;
        double x_position = (x_count - 1) * x_increment;
        if (y_count % 2 == 0) {
          if (z_count % 2 == 0) {
            x_position -= x_increment / 2.0;
          } else {
            x_position += x_increment / 2.0;
          }
        }
        if (z_count % 2 == 0) x_position += x_increment / 2.0;
        if (x_position > input_domain_size[0]) break;
        (*initial_position).push_back(x_position);
        (*initial_position).push_back(y_position);
        (*initial_position).push_back(z_position);
      }
    }
  }
}

// BuildFromFile
// *****************************************************************************
// (independent)
// Reads atom positions from a file formatted as:
// atom_number
// atom_1_x atom_1_y atom_1_z
// atom_2_x atom_2_y atom_2_z
// ...
// y & z or z dimensions are omitted in 1-D and 2-D cases, respectively.
void BuildFromFile(string input_configuration_file,
                   vector<double>* initial_position) {
  input_configuration_file.insert(0, "input/");
  FILE* input_configuration_file_ptr = fopen(input_configuration_file.c_str(),
                                             "r");
  if (input_configuration_file_ptr == NULL) {
    printf("cannot open the initial configuration file: %s\n",
           input_configuration_file.c_str());
    exit(1);
  }
  int atom_num;
  fscanf(input_configuration_file_ptr, "%d", &atom_num);
  for (int i = 0; i < atom_num; i++) {
    if (input_configuration_file_ptr == NULL) {
      printf("inconsistent atom num in the initial configuration file\n");
      exit(1);
    }
    for (int j = 0; j < 3; j++) {
      double position;
      fscanf(input_configuration_file_ptr, "%lf", &position);
      (*initial_position).push_back(position);
    }
  }
  fclose(input_configuration_file_ptr);
}

