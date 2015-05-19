// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// output.h
// *****************************************************************************
// Consists of Output class that is responsible of logging, throwing some
// errors, and writing to files. In addition to the construction of the class,
// it needs to be initialized by the Init function and finalized by the Final
// function.

#ifndef MMM_V14_6_OUTPUT_H_
#define MMM_V14_6_OUTPUT_H_

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mediator.h"

using std::ofstream;
using std::stringstream;

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Output : public Mediator {
 public:
  explicit Output() : Mediator() {}
  ~Output() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Output(Output&) = delete;
  // Output& operator=(const Output&) = delete;

  // Finalizes the class.
  void Final();
  // Initializes the class.
  void Init();
  // Outputs some data to files.
  void OutputAllToFile();
  // Outputs input file to the log file.
  void OutputInputFileToLog();
  // Outputs time information to the log file. In particular, outputs
  // simulation time, average parallel time, and parallel time ratio statistics.
  void OutputTimeToLog();
  // Outputs input statement to the log file and continues to a new line.
  void OutputToLog(string statement);
  // Outputs input statement to the log file and stays in the current line.
  void OutputToLog(bool is_new_line, string statement);
  // Throws an error and exits the program.
  void ThrowError(string statement);
  // Converts the input type to a string.
  template<typename T>
  string ToString(T value) const {
    stringstream string_stream;
    string_stream << value;
    return string_stream.str();
  }

  // Accessor and mutator functions:
  // frame_count_
  int frame_count() const { return frame_count_; }
  void set_frame_count(int frame_count) { frame_count_ = frame_count; }
  void add_frame_count() { frame_count_++; }
  // file_update_frequency_
  int file_update_frequency() const { return file_update_frequency_; }
  void set_file_update_frequency(int file_update_frequency) {
    file_update_frequency_ = file_update_frequency;
  }
  // screen_update_frequency_
  int screen_update_frequency() const { return screen_update_frequency_; }
  void set_screen_update_frequency(int screen_update_frequency) {
    screen_update_frequency_ = screen_update_frequency;
  }

 private:
  // Outputs data of some selected variables to a file. See .cc file for format
  // details.
  void OutputSelectedVariable();
  // Outputs by trimming the statement to the 80 columns of horizontal width.
  void TrimOutputStatement(string statement);

  // Accessor and mutator functions:
  // lammps_trajectory_file_ptr_
  FILE* lammps_trajectory_file_ptr() { return lammps_trajectory_file_ptr_; }
  void set_lammps_trajectory_file_ptr(FILE* lammps_trajectory_file_ptr) { 
    lammps_trajectory_file_ptr_ = lammps_trajectory_file_ptr; 
  }
  // log_file_object_
  ofstream* log_file_object() { return &log_file_object_; }
  // mesh_file_ptr_
  FILE* mesh_file_ptr() { return mesh_file_ptr_; }
  void set_mesh_file_ptr(FILE* mesh_file_ptr) { 
    mesh_file_ptr_ = mesh_file_ptr; 
  }
  // output_file_namebase;_
  string output_file_namebase() const { return output_file_namebase_; }
  void set_output_file_namebase(string output_file_namebase) {
    output_file_namebase_ = output_file_namebase;
  }
  // selected_variable_file_ptr_
  FILE* selected_variable_file_ptr() { return selected_variable_file_ptr_; }
  void set_selected_variable_file_ptr(FILE* selected_variable_file_ptr) { 
    selected_variable_file_ptr_ = selected_variable_file_ptr;
  }

  FILE* lammps_trajectory_file_ptr_;  // lammps trajectory
  FILE* mesh_file_ptr_;               // mesh file
  FILE* selected_variable_file_ptr_;  // selected variables
  int frame_count_;                   // count of frames (time of output of all)
  int file_update_frequency_;         // frequency of file update
  int screen_update_frequency_;       // frequency of screen update
  ofstream log_file_object_;          // log file
  string output_file_namebase_;       // namebase to output files
};

#endif  // MMM_V14_6_OUTPUT_H_

