// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// input.cc
// *****************************************************************************

#include "input.h"

#include <mpi.h>

#include <cstdlib>
#include <string>

#include "mmm.h"
#include "output.h"

using std::string;

// Public
// *****************************************************************************

// ReadInputFile
// *****************************************************************************
// Reads input file into vector<vector<string>> command_table. In that, some
// rules are followed:
// - Line length is limited to 100 characters including everything.
// - Lines that start with pound ("#") or space (" ") are ignored.
// - Numbers in scientific notation cannot be read.
// - Selection ID style must use colons for sequences (e.g., ID 2:2:20).
// - Selection ID style accepts input of multiple blocks (e.g., ID 1 2 4:2:8).
void Input::ReadInputFile() {  
  set_input_file("input/" + mmm()->name() + ".in");
  FILE* input_file_ptr = fopen(input_file().c_str(), "r");
  if (input_file_ptr == NULL) {
    output()->ThrowError("cannot open the input file: " + input_file());
  }
  char char_line[100];
  while (fgets(char_line, 100, input_file_ptr) != NULL) {
    set_line(char_line);
    string first_char = line_sub(0, 1);
    if (first_char.find_first_not_of("\t\v\r\n") != 0 ||
        first_char.compare(" ") == 0 || first_char.compare("#") == 0) {
      continue;
    }
    ParseLine();      
    push_command_table(parsed_line());
  }
  fclose(input_file_ptr);
}

// Private
// *****************************************************************************

// ParseLine
// *****************************************************************************
// Right and left bounds are not included in the block, e.g., for "0123" at the
// start of a line, left bound is -1 and right bound is 4.
void Input::ParseLine() {
  clear_parsed_line();
  bool is_end_of_line = false;
  int left_bound = -1;
  while (!is_end_of_line) {
    // Get block
    int right_bound = line_find_first_not_of(valid_char, left_bound + 1);
    int block_length = right_bound - left_bound - 1;
    string parse = line_sub(left_bound + 1, block_length);
    // Check if end of line
    size_t found = line_find_first_of(valid_char, right_bound + 1);
    if (found == string::npos || right_bound == -1) is_end_of_line = true;
    push_parsed_line(parse);
    left_bound = right_bound;
  }
}


