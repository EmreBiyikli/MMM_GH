// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// output.cc
// *****************************************************************************

#include "output.h"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "atom_group.h"
#include "input.h"
#include "matrix.h"
#include "mesh.h"
#include "mmm.h"
#include "model.h"
#include "neighbor.h"
#include "time.h"

using std::cerr;
using std::endl;
using std::ios;
using std::string;
using std::vector;

// Independent functions (by means of unnamed namespace)
// *****************************************************************************

namespace {

// OutputLammpsTrajectory
// *****************************************************************************
// (independent)
// Appends data to the input file in the lammps trajectory file in the format:
// ITEM: TIMESTEP
// #
// ITEM: NUM OF ATOMS
// #
// ITEM: BOX BOUNDS
// #
// ITEM: ATOMS
// atom_1_id atom_1_type atom_1_x atom_1_y atom_1_z
// atom_2_id atom_2_type atom_2_x atom_2_y atom_2_z
// ...
// In case of 1-D and 2-D, missin bounds are set to -1.0 and +1.0 and missing
// positions are set to 0.0.
void OutputLammpsTrajectory(const vector<int>& type,
                            const vector<double>& domain_boundary,
                            const vector<double>& position,
                            int timestep,
                            int dimension,
                            FILE* lammps_trajectory_file_ptr) {
  fprintf(lammps_trajectory_file_ptr, "ITEM: TIMESTEP\n");
  fprintf(lammps_trajectory_file_ptr, "%d\n", timestep);
  fprintf(lammps_trajectory_file_ptr, "ITEM: NUMBER OF ATOMS\n");
  fprintf(lammps_trajectory_file_ptr, "%d\n", position.size() / 3);
  fprintf(lammps_trajectory_file_ptr, "ITEM: BOX BOUNDS\n");
  for (int i = 0; i < 3; i++) {
    if (i < dimension) {
      fprintf(lammps_trajectory_file_ptr, "%.8f %.8f\n", domain_boundary[2 * i],
              domain_boundary[2 * i + 1]);
    } else {
      fprintf(lammps_trajectory_file_ptr, "-1.0 1.0\n");
    }
  }
  fprintf(lammps_trajectory_file_ptr, "ITEM: ATOMS\n");
  for (int i = 0; i < position.size() / 3; i++) {
    fprintf(lammps_trajectory_file_ptr, "%d %d", i, type[i]);
    for (int j = 0; j < 3; j++) {
      if (j < dimension) {
        fprintf(lammps_trajectory_file_ptr, " %.8f", position[3 * i + j]);
      } else {
        fprintf(lammps_trajectory_file_ptr, " 0.0");
      }
    }
    fprintf(lammps_trajectory_file_ptr, "\n");
  }
}

// OutputMesh
// *****************************************************************************
// (independent)
// Appends data to the input file inthe mesh file in the format:
// current_iteration element_number
// element_1_node_1 element_1_node_2 ...
// ...
// The line continues until element_node_number and file until element_number.
void OutputMesh(const vector<int>& element_node,
                int element_node_num,
                int element_num,
                int iteration,
                FILE* mesh_file_ptr) {
  fprintf(mesh_file_ptr, "%d %d\n", iteration, element_num);
  for (int i = 0; i < element_num; i++) {
    for (int j = 0; j < element_node_num; j++) {
      fprintf(mesh_file_ptr, "%d ", element_node[element_node_num * i + j]);
    }
    fprintf(mesh_file_ptr, "\n");
  }
}

}  // namespace

// Public
// *****************************************************************************

// Final
// *****************************************************************************
void Output::Final() {
  OutputToLog("empty line");
  OutputToLog("simulation stats");
  OutputToLog("dashed line");
  OutputToLog("neighbor build count: " + ToString(neighbor()->build_count()));
  OutputToLog("mesh build count: " + ToString(mesh()->build_count()));
  OutputTimeToLog();
  OutputToLog("empty line");
  OutputToLog("-- end of simulation");
  if (mmm()->processing_element() == 0) {
    log_file_object()->close();
    fclose(lammps_trajectory_file_ptr());
    fclose(mesh_file_ptr());
    fclose(selected_variable_file_ptr());
    if (model()->is_full_atomistic()) {
      remove((output_file_namebase() + "_mesh.txt").c_str());
    }
  }
}

// Init
// *****************************************************************************
void Output::Init() {
  set_frame_count(0);
  set_output_file_namebase("output/" + mmm()->name());
  if (mmm()->processing_element() == 0) {
    log_file_object()->open((output_file_namebase() + "_log.txt").c_str(),
                            ios::out);
    set_mesh_file_ptr(fopen((output_file_namebase() + "_mesh.txt").c_str(), 
                            "w"));
    set_lammps_trajectory_file_ptr(fopen((output_file_namebase() + 
                                             ".lammpstrj").c_str(), "w"));
    set_selected_variable_file_ptr(fopen((output_file_namebase() +
                                             ".txt").c_str(), "w"));
  }
  OutputToLog("simulation name: " + mmm()->name());
  OutputToLog("dashed line");
  OutputToLog("number of processing elements: " +
      ToString(mmm()->processing_element_num()));
}

// OutputAllToFile
// *****************************************************************************
void Output::OutputAllToFile() {
  if (mmm()->processing_element() == 0) {
    if (mmm()->current_iteration() % file_update_frequency() == 0 ||
        !mmm()->is_iteration_loop_continue()) {
      add_frame_count();
      OutputLammpsTrajectory(model()->atom_type_.get(), 
                             mmm()->domain_boundary(),
                             atom_group()->position_.get(), frame_count() - 1,
                             mmm()->dimension(), lammps_trajectory_file_ptr());
      OutputSelectedVariable();
      if (!model()->is_full_atomistic()) {
        OutputMesh(mesh()->element_node(), mesh()->element_node_num(),
                   mesh()->element_num(), mmm()->current_iteration(),
                   mesh_file_ptr());
      }
    }
  }
}

// ThrowError
// *****************************************************************************
void Output::ThrowError(string statement) {
  if (mmm()->processing_element() == 0) {
    *log_file_object() << statement << endl;
    cerr << statement << endl;
  }
  exit(1);
}

// OutputInputFileToLog
// *****************************************************************************
void Output::OutputInputFileToLog() {
  if (mmm()->processing_element() == 0) {
    OutputToLog("empty line");
    OutputToLog("input file");
    OutputToLog("dashed line");
    FILE* input_file_ptr = fopen(input()->input_file().c_str(), "r");
    char input_file_line[100];
    while (fgets(input_file_line, 100, input_file_ptr) != NULL) {
      OutputToLog(false, ToString(input_file_line));
    }
    fclose(input_file_ptr);
    OutputToLog("dashed line");
    OutputToLog("empty line");
  }
}

// OutputTimeToLog
// *****************************************************************************
void Output::OutputTimeToLog() {
  // Title
  OutputToLog("empty line");
  OutputToLog("time");
  OutputToLog("dashed line");
  // Simulation time
  double hour, minute, second;
  Time::ConvertTimeInSecondToHumanReadable(time()->simulation_time(), &hour,
                                           &minute, &second);
  OutputToLog("simulation time: " + ToString(hour) + " h " + ToString(minute) +
                  " m " + ToString(second) + " s ");
  time()->ComputeParallelTimeStat();
  // Parallel time
  if (mmm()->processing_element() == 0) {
    time()->set_average_parallel_time(time()->parallel_time_ratio_stat(1) *
                                          time()->simulation_time());
    double hour, minute, second;
    Time::ConvertTimeInSecondToHumanReadable(time()->average_parallel_time(),
                                             &hour, &minute, &second);
    OutputToLog("average parallel time: " + ToString(hour) + " h " +
                    ToString(minute) + " m " + ToString(second) + " s " + "(" +
                    ToString(100 * time()->parallel_time_ratio_stat(1)) + "%)");
    OutputToLog("paralell time distribution (min average max): " +
                    ToString(time()->parallel_time_ratio_stat(0)) + " " +
                    ToString(time()->parallel_time_ratio_stat(1)) + " " +
                    ToString(time()->parallel_time_ratio_stat(2)));
  }
}

// OutputToLog
// *****************************************************************************
void Output::OutputToLog(string statement) {
  TrimOutputStatement(statement);
  *log_file_object() << endl;
}
void Output::OutputToLog(bool is_new_line, string statement) {
  TrimOutputStatement(statement);
  if (is_new_line) *log_file_object() << endl;
}

// Private
// *****************************************************************************

// OutputSelectedVariable
// *****************************************************************************
// Appends data into the current file in the format:
// current_iteration atom_number 0 0 0 0 0 0 0
// atom_1: id(1) type(2) position(3-5) energy(6) force(7-9)
// atom_2: id(1) type(2) position(3-5) energy(6) force(7-9)
// ...
void Output::OutputSelectedVariable() {
  fprintf(selected_variable_file_ptr(), "%d %d 0 0 0 0 0 0 0\n",
          mmm()->current_iteration(), atom_group()->atom_num());
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    fprintf(selected_variable_file_ptr(),
            "%d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
            i,
            model()->atom_type_.get_1d_element(i),
            atom_group()->position_.get_element(i, 0),
            atom_group()->position_.get_element(i, 1),
            atom_group()->position_.get_element(i, 2),
            atom_group()->energy_.get_1d_element(i),
            atom_group()->force_.get_element(i, 0),
            atom_group()->force_.get_element(i, 1),
            atom_group()->force_.get_element(i, 2));
  }
}

// TrimOutputStatement
// *****************************************************************************
// Does not insert a new line at the end.
void Output::TrimOutputStatement(string statement) {
  if (mmm()->processing_element() == 0) {
    if (statement.compare("dashed line") == 0) {
      *log_file_object() << "----------------------------------------";
      *log_file_object() << "----------------------------------------";
    } else if (statement.compare("empty line") == 0) {
      *log_file_object();
    } else {
      // Normal statements
      if(statement.size() <= 80) {
        *log_file_object() << statement;
      } else {
        while (statement.size() > 80) {
          string to_output = statement.substr(0, 80);
          *log_file_object() << to_output << endl;
          statement = statement.substr(80);
        }
        *log_file_object() << statement;
      }
    }
  }
}