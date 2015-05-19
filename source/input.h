// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// input.h
// *****************************************************************************
// Consists of Input class that is responsible of reading the simulation input
// file, parsing it, and saving it in the command table. All are performed by
// the ReadInputFile function.

#ifndef MMM_V14_6_INPUT_H_
#define MMM_V14_6_INPUT_H_

#include <string>
#include <vector>

#include "mediator.h"

using std::string;

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Input : public Mediator {
 public:
  explicit Input() : Mediator() {
    valid_char = 
        "qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM1234567890+-._:";
  }
  ~Input() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Input(Input&) = delete;
  // Input& operator=(const Input&) = delete;

  // See comment at top of file for a complete description.
  void ReadInputFile();

  // Accessor and mutator functions:
  // command_line_
  const vector<string>& command_line() const { return command_line_; }
  string command_line_at(int index) const { return command_line_[index]; }
  void set_command_line(vector<string> command_line) {
    command_line_ = command_line;
  }
  // command_table_
  const vector< vector<string> >& command_table() const { return command_table_; }
  const vector<string>& command_table_at(int index) const {
    return command_table_[index];
  }
  void push_command_table(vector<string> line) {
    command_table_.push_back(line);
  }
  // input_file_
  string input_file() const { return input_file_; }
  void set_input_file(string input_file) { input_file_ = input_file; }

 private:
  // Parses the current line_ into parsed_line_.
  void ParseLine();

  // Accessor and mutator functions:
  // line_
  const string* line_address() const { return &line_; }
  string line() const { return line_; }
  string line_sub(int start, int end) const { return line_.substr(start, end); }
  int line_find_first_of(string to_find, int start) const {
    return line_.find_first_of(to_find, start);
  }
  int line_find_first_not_of(string to_find, int start) const {
    return line_.find_first_not_of(to_find, start);
  }
  void set_line(string line) { line_ = line; }
  // parsed_line_
  const vector<string>& parsed_line() const { return parsed_line_; }
  void clear_parsed_line() { parsed_line_.clear(); }
  void push_parsed_line(string parse) { parsed_line_.push_back(parse); }

  string input_file_;                     // input filename
  string line_;                           // current input file line
  // characters those are considered valid when reading the input file
  string valid_char;
  vector<string> command_line_;           // current command line
  vector<string> parsed_line_;            // current parsed line
  vector< vector<string> > command_table_;  // command table
};

#endif  // MMM_V14_6_INPUT_H_

