// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// build_initial_configuration.h
// ****************************************************************************
// (independent)
// Consists of 4 standalone functions those, depending on the choice, builds an
// initial configuration. The choices are 1-D uniform, 2-D hexagonal, 3-D face-
// centered cubic (FCC), or N-D file. The last choice reads an N-D configuration
// from the input file.

#ifndef MMM_V14_6_BUILD_INITIAL_CONFIGURATION_H_
#define MMM_V14_6_BUILD_INITIAL_CONFIGURATION_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

// Builds a 1-D uniform configuration. Input domain size and input initial
// spacing indicate the size of the domain and space between atoms,
// respectively. The configuration is returned by initial_position parameter.
void Build1DUniform(const vector<double>& input_domain_size,
                    double initial_spacing, vector<double>* initial_position);
// Builds a 2-D hexagonal configuration. Input domain size and input initial
// spacing indicate the size of the domain and space between atoms,
// respectively. The configuration is returned by initial_position parameter.
void Build2DHexagonal(const vector<double>& input_domain_size,
                      double initial_spacing, vector<double>* initial_position);
// Builds a 3-D face-centered cubic (FCC) configuration. Input domain size and
// input initial spacing indicate the size of the domain and space between
// atoms, respectively. The configuration is returned by initial_position
// parameter.
void Build3DFCC(const vector<double>& input_domain_size, double initial_spacing,
                vector<double>* initial_position);
// Builds an N-D configuration read from the input file. The configuration is
// returned by initial_position parameter.
void BuildFromFile(string input_configuration_file,
                   vector<double>* initial_position);

#endif  // MMM_V14_6_BUILD_INITIAL_CONFIGURATION_H_

