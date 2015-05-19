// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// add_in.h
// *****************************************************************************
// Consists of AddIn class that handles the add-ins. Add-ins are functions those
// are controlled in this header file, implemented in the add_in.cc file, and 
// called somewhere in the code. Add-in function calls should be covered with 
// appropriate comments to  warn the user. In case of a modification to the 
// add-ins, the program should be re-built.

#ifndef MMM_V14_6_AddIn_H_
#define MMM_V14_6_AddIn_H_

#include <string>
#include <vector>

#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class AddIn : public Mediator {
 public:  
  explicit AddIn() : Mediator() {
    // Control switches to add-in functions.
    is_on_static_paper_crack_displacement_ = false;
  }
  ~AddIn() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // AddIn(AddIn&) = delete;
  // AddIn& operator=(const AddIn&) = delete;

  // Add-in functions:
  // Adds displacement in the static paper crack example.
  void StaticPaperCrackDisplacement();
  bool is_on_static_paper_crack_displacement() {
    return is_on_static_paper_crack_displacement_;
  }

 private:
  bool is_on_static_paper_crack_displacement_;
};

#endif  // MMM_V14_6_AddIn_H_

