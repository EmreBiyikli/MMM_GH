// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// main.cc
// *****************************************************************************
// Run with the simulation name, it constructs, initializes, simulates, and
// finalizes MMM. It also initializes and finalizes MPI before and after calls
// of MMM.

#include <mpi.h>

#include <string>

#include "mmm.h"

using std::string;

// See comment at top of file for a complete description.
int main(int argc, char **argv) {
  // This should normally be simulation_name = argv[1];
  string simulation_name = "dynamic_paper_2d_wave";
  MPI_Init(&argc, &argv);
  MMM mmm;
  mmm.Init(simulation_name);
  mmm.Simulate();
  mmm.Final();
  MPI_Finalize();
}

