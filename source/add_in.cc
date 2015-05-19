// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// add_in.cc
// *****************************************************************************

#include "add_in.h"

#include <stdio.h>

#include "atom_group.h";

// Public
// *****************************************************************************

// StaticPaperCrackDisplacement
// *****************************************************************************
// Adds a displacement of +-0.002 at the ends in every 100 iterations. Right end
// consists of atoms 0:38 and left end consists of atoms 1018:1056.
void AddIn::StaticPaperCrackDisplacement() {
  if (!is_on_static_paper_crack_displacement()) return;
  static int step_num = 0;
  if (mmm()->current_iteration() == 0) return;
  if (mmm()->is_iteration_loop_continue() == true) return;  
  if (step_num > 1200) return;
  step_num ++;  
  mmm()->set_is_iteration_loop_continue(true);
  for (int i = 0; i <= 38; i++) {
    atom_group()->position_.add_element(i, 0, -0.002);
  }
  for (int i = 1018; i <= 1056; i++) {
    atom_group()->position_.add_element(i, 0, 0.002);
  }    
}
