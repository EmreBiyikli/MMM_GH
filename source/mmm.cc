// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// mmm.cc
// *****************************************************************************

#include "mmm.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "add_in.h"
#include "atom_group.h"
#include "build_initial_configuration.h"
#include "conjugate_gradient.h"
#include "input.h"
#include "mesh.h"
#include "model.h"
#include "neighbor.h"
#include "output.h"
#include "potential.h"
#include "select.h"
#include "temperature.h"
#include "time.h"
#include "velocity_verlet.h"

using std::max;
using std::string;
using std::vector;

// Independent functions
// *****************************************************************************

namespace {

// Partition
// *****************************************************************************
// (independent)
// Partitions the atoms, ID'd from 1 to input atom_num, into the input number of
// processing elements. Partition is returned in output processing element where
// each element stores the processing element of the corresponding atom. For
// example, processing_element[i] = pe; means that atom i belongs to processing
// element pe.
void Partition(int atom_num, int processing_element_num,
               vector<int>* processing_element) {
  processing_element->resize(atom_num);
  double atom_num_per_processing_element = static_cast<double>(atom_num) /
      static_cast<double>(processing_element_num);
  for (int i = 0; i < processing_element_num; i++) {
    int lower_bound = static_cast<int>(static_cast<double>(i) *
                                           atom_num_per_processing_element);
    for (int j = 0; j < atom_num; j++) {
      if (j >= lower_bound) (*processing_element)[j] = i;
    }
  }
}

}  // namespace

// Constructor
// *****************************************************************************
// It is important to follow the order:
// 1. Construct objects
// 2. Call InitPointer of those inherit from Mediator
MMM::MMM() {
  set_is_dynamic(false);
  set_is_MPI(false);
  add_in_ = new AddIn();
  atom_group_ = new AtomGroup();
  conjugate_gradient_ = new ConjugateGradient();
  input_ = new Input();
  mesh_ = new Mesh;
  model_ = new Model();
  neighbor_ = new Neighbor;
  output_ = new Output();
  potential_ = new Potential();
  select_ = new Select();
  temperature_ = new Temperature();
  time_ = new Time();
  velocity_verlet_ = new VelocityVerlet();
  add_in()->InitPointer(this);
  conjugate_gradient()->InitPointer(this);
  input()->InitPointer(this);
  model()->InitPointer(this);
  output()->InitPointer(this);
  potential()->InitPointer(this);
  temperature()->InitPointer(this);
  time()->InitPointer(this);
  velocity_verlet()->InitPointer(this);
}

// Destructor
// *****************************************************************************
MMM::~MMM() {
  delete add_in_;
  delete atom_group_;
  delete conjugate_gradient_;
  delete input_;
  delete mesh_;
  delete model_;
  delete neighbor_;
  delete output_;
  delete potential_;
  delete select_;
  delete temperature_;
  delete time_;
  delete velocity_verlet_;
}

// Public
// *****************************************************************************

// Final
// *****************************************************************************
void MMM::Final() {
  output()->Final();
}

// Init
// *****************************************************************************
void MMM::Init(string simulation_name) {
  set_name(simulation_name);
  MPI_Comm_size(MPI_COMM_WORLD, processing_element_num_address());
  MPI_Comm_rank(MPI_COMM_WORLD, processing_element_address());  
  if (processing_element_num() > 1) set_is_MPI(true);  
  output()->Init();  
  input()->ReadInputFile();
  output()->OutputInputFileToLog();
}

// Simulate
// *****************************************************************************
// Executes input file command_line-by-command_line.
void MMM::Simulate() {
  time()->start_simulation_time();
  for (int i = 0; i < input()->command_table().size(); i++) {
    input()->set_command_line(input()->command_table_at(i));
    string key = input()->command_line_at(0);
    if (key.compare("dimension") == 0) {
      CommandDimension();
    } else if (key.compare("initial_configuration") == 0) {
      CommandInitialConfiguration();
    } else if (key.compare("load") == 0) {
      CommandLoad();
    } else if (key.compare("mesh") == 0) {
      CommandMesh();
    } else if (key.compare("model") == 0) {
      CommandModel();
    } else if (key.compare("neighbor") == 0) {
      CommandNeighbor();
    } else if (key.compare("output") == 0) {
      CommandOutput();
    } else if (key.compare("potential") == 0) {
      CommandPotential();
    } else if (key.compare("run") == 0) {
      CommandRun();
    } else if (key.compare("select") == 0) {
      CommandSelect();
    } else if (key.compare("temperature") == 0) {
      CommandTemperature();
    } else if (key.compare("type") == 0) {
      CommandType();
    } else {
      output()->ThrowError("incorrect command key");
    }
  }
  time()->finish_simulation_time();
}

// Private
// *****************************************************************************

// CheckLostAtom
// *****************************************************************************
void MMM::CheckLostAtom() {
  // Allowed distance
  double extension_factor = 10;
  double allowed_distance = 0;
  vector<double> length(3, 0);
  for (int i = 0; i < dimension(); i++) {
    length[i] = domain_boundary_at(2 * i + 1) - domain_boundary_at(2 * i);
    allowed_distance = max(allowed_distance, extension_factor * length[i]);
  }
  // Check
  bool is_stop = false;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    for (int j = 0; j < dimension(); j++) {
      double global_displacement = atom_group()->position_.get_element(i, j) -
          atom_group()->initial_position_.get_element(i, j);
      is_stop = global_displacement > allowed_distance ? true : false;
      if (is_stop) break;
    }
  }
  // Stop
  if (is_stop) output()->ThrowError("lost atom");
}

// CommandDimension
// *****************************************************************************
// Input format:
// int
// dimension
void MMM::CommandDimension() {
  set_dimension(atof(input()->command_line_at(1).c_str()));
  output()->OutputToLog("dimension is set to " +
                        output()->ToString(dimension()));
}

// CommandInitialConfiguration
// *****************************************************************************
// Input format:
// double int string double []
// mass configuration_style initial_spacing domain_size
// configuration_style_: file/uniform/hexagonal/fcc
// i.e., mass file configuration_filename
// i.e., mass uniform domain_size(0)
// i.e., mass hexagonal domain_size(0) domain_size(1)
// i.e., mass fcc domain_size(0) domain_size(1) domain_size(2)
void MMM::CommandInitialConfiguration() {
  set_mass(atof(input()->command_line_at(1).c_str()));
  double initial_spacing;
  if (atof(input()->command_line_at(3).c_str()) == 0) {
    initial_spacing = potential()->equilibrium_spacing();
  } else {
    initial_spacing = atof(input()->command_line_at(3).c_str());
  }
  vector<double> initial_configuration;  
  if (input()->command_line_at(2).compare("file") != 0) {
    vector<double> input_domain_size;
    input_domain_size.push_back(atof(input()->command_line_at(4).c_str()));
    if (input()->command_line().size() > 5) {
      input_domain_size.push_back(atof(input()->command_line_at(5).c_str()));
      if (input()->command_line().size() > 6) {
        input_domain_size.push_back(atof(input()->command_line_at(6).c_str()));
      }
    }
    if (input()->command_line_at(2).compare("uniform") == 0) {
      Build1DUniform(input_domain_size, initial_spacing,
                     &initial_configuration);
    } else if (input()->command_line_at(2).compare("hexagonal") == 0) {
      Build2DHexagonal(input_domain_size, initial_spacing,
                       &initial_configuration);
    } else if (input()->command_line_at(2).compare("fcc") == 0) {
      Build3DFCC(input_domain_size, initial_spacing, &initial_configuration);
    }
  } else {
    BuildFromFile(input()->command_line_at(3), &initial_configuration);
  }
  atom_group()->Init(initial_configuration.size() / 3);
  for (int i = 0; i < atom_group()->degree_of_freedom_num(); i++) {
    atom_group()->position_.set(initial_configuration);
  }
  select()->Init(atom_group()->position_.get());
  set_domain_boundary(atom_group()->GetBoundary());
  vector<int> processing_element;
  Partition(atom_group()->atom_num(), processing_element_num(),
            &processing_element);
  atom_group()->processing_element_.set(processing_element);
  output()->OutputToLog("initial configuration is built with " +
                            output()->ToString(atom_group()->atom_num()) +
                            " atoms on a domain of size " +
                            output()->ToString(domain_boundary_at(0)) + " " +
                            output()->ToString(domain_boundary_at(1)) + " " +
                            output()->ToString(domain_boundary_at(2)) + " " +
                            output()->ToString(domain_boundary_at(3)) + " " +
                            output()->ToString(domain_boundary_at(4)) + " " +
                            output()->ToString(domain_boundary_at(5)));
}

// CommandLoad
// *****************************************************************************
// Input format:
// string int string double double double
// load_type_1 selection load_type_2 load_value[0] load_value[1] load_value[2]
// load_type_1: displacement/force
// load_type_2: set/add
// note: NULL should be used for axes without load
// e.g., displacement 1 set 0 NULL NULL
// e.g., force 2 add NULL -1.0 NULL
void MMM::CommandLoad() {
  // Command input line
  string load_type_parsel_1 = input()->command_line_at(1);
  int input_selection = atoi(input()->command_line_at(2).c_str());
  string load_type_parsel_2 = input()->command_line_at(3);
  vector<bool> is_axes_loaded(3, false);
  vector<double> load_value(3, 0);
  if (input()->command_line_at(4).compare("NULL") != 0) {
    is_axes_loaded[0] = true;
    load_value[0] = atof(input()->command_line_at(4).c_str());
  }
  if (input()->command_line_at(5).compare("NULL") != 0) {
    is_axes_loaded[1] = true;
    load_value[1] = atof(input()->command_line_at(5).c_str());
  }
  if (input()->command_line_at(6).compare("NULL") != 0) {
    is_axes_loaded[2] = true;
    load_value[2] = atof(input()->command_line_at(6).c_str());
  }
  // Express load type in number
  int load_type;
  if (load_type_parsel_1.compare("displacement") == 0) {
    if (load_type_parsel_2.compare("set") == 0) {
      load_type = 1;
    } else if (load_type_parsel_2.compare("add") == 0) {
      load_type = 2;
    }
  } else if (load_type_parsel_1.compare("force") == 0) {
    if (load_type_parsel_2.compare("set") == 0) {
      load_type = 3;
    } else if (load_type_parsel_2.compare("add") == 0) {
      load_type = 4;
    }
  }
  // Impose loading to atom_group
  for (int i = 0; i < select()->selection(input_selection).size(); i++) {
    int atom = select()->selection_at(input_selection, i);
    if (!model()->is_rep_.get_1d_element(atom)) {
      output()->ThrowError("type of atom " + output()->ToString(atom) +
                                " is not valid for loading");
    }
    for (int j = 0; j < 3; j++) {
      if (is_axes_loaded[j]) {
        atom_group()->is_load_.set_1d_element(atom, true);
        atom_group()->load_type_.set_element(atom, j, load_type);
        atom_group()->load_value_.set_element(atom, j, load_value[j]);
      }
    }
  }
  // Log
  output()->OutputToLog(load_type_parsel_1 + " loading of " + 
                            input()->command_line_at(4) +
                            " " + input()->command_line_at(5) + " " +
                            input()->command_line_at(6) + " is " +
                            load_type_parsel_2 + " to selection " +
                            output()->ToString(input_selection));
}

// CommandMesh
// *****************************************************************************
// Input format:
// string int/string
// style style_inputs
// style: qhull/file
// i.e., qhull selection
// i.e., file filename
// e.g., mesh qhull 1
// e.g., mesh file coarse_mesh()->txt
void MMM::CommandMesh() {
  if (input()->command_line_at(1).compare("qhull") == 0) {
    int selection_ID = atoi(input()->command_line_at(2).c_str());
    vector<double> selection_position;
    select()->GetPosition(selection_ID, &selection_position);
    mesh()->InitQhull(is_MPI(), selection_position,
                      select()->selection(selection_ID), dimension());
  } else if (input()->command_line_at(1).compare("file") == 0) {
    mesh()->InitFile(is_MPI(),
                     dimension(), "input/" + input()->command_line_at(2));
  }
  mesh()->BuildMesh();
  output()->OutputToLog("mesh is built with " +
                            output()->ToString(mesh()->element_num()) +
                            " elements and " +
                            output()->ToString(mesh()->node_num()) + " nodes");
}

// CommandModel
// *****************************************************************************
// Input format:
// string string
// style scheme
// style: full_atomistic/MMM
// scheme: no_ssmp/all_ssmp/ssmp_around_irep/ssmp_around_psmp
void MMM::CommandModel() {
  if (input()->command_line_at(1).compare("full_atomistic") == 0) {
    model()->InitFullAtomistic();
    model()->BuildModel();
    output()->OutputToLog("model is built with full_atomistic style");
  } else if (input()->command_line_at(1).compare("mmm") == 0) {
    model()->InitMMM(input()->command_line_at(2));
    model()->BuildModel();
    string temp;
    temp = "model is built with mmm style " + model()->scheme() +
           " scheme, atom type numbers are ";
    for (int i = 0; i < 5; i++) {
      temp.append(
          output()->ToString(model()->atom_type_num_.get_1d_element(i)) + " ");
    }
    output()->OutputToLog(temp);
  }
}

// CommandNeighbor
// *****************************************************************************
// Input format:
// double string []
// neighbor_cutoff_radius update_style update_parameter
// update_style: auto/every
// i.e., neighbor_cutoff_radius auto
// i.e., neighbor_cutoff_radius every update_frequency
void MMM::CommandNeighbor() {
  if (input()->command_line_at(2).compare("auto") == 0) {
    neighbor()->InitAuto(atof(input()->command_line_at(1).c_str()),
                         potential()->cutoff_radius());
  } else if (input()->command_line_at(2).compare("every") == 0) {
    neighbor()->InitEvery(atof(input()->command_line_at(1).c_str()),
                          potential()->cutoff_radius(),
                          atoi(input()->command_line_at(3).c_str()));
  }
  if (input()->command_line_at(2).compare("auto") == 0) {
    output()->OutputToLog("neighbor is initialized with auto style " +
        output()->ToString(neighbor()->cutoff_radius()) + " cut-off radius");
  } else if (input()->command_line_at(2).compare("every") == 0) {
    output()->OutputToLog("neighbor is initialized with every " + 
        output()->ToString(neighbor()->update_frequency()) + " style " +
        output()->ToString(neighbor()->cutoff_radius()) + " cut-off radius");
  }
}

// CommandOutput
// *****************************************************************************
// Input format:
// int int
// file_update_frequency screen_update_frequency
void MMM::CommandOutput() {
  output()->set_file_update_frequency(
      atoi(input()->command_line_at(1).c_str()));
  output()->set_screen_update_frequency(
      atoi(input()->command_line_at(2).c_str()));
  output()->OutputToLog("file update frequency is " +
      output()->ToString(output()->file_update_frequency()) +
      " and screen update frequency is " +
      output()->ToString(output()->screen_update_frequency()));
}

// CommandPotential
// *****************************************************************************
// Input format:
// string double []
// potential potential_cutoff_radius potential_parameter
void MMM::CommandPotential() {
  vector<double> parameter;
  if (input()->command_line_at(1).compare("lennard_jones") == 0) {
    parameter.push_back(atof(input()->command_line_at(3).c_str()));
    parameter.push_back(atof(input()->command_line_at(4).c_str()));
  } else if (input()->command_line_at(1).compare("morse") == 0) {
    parameter.push_back(atof(input()->command_line_at(3).c_str()));
    parameter.push_back(atof(input()->command_line_at(4).c_str()));
    parameter.push_back(atof(input()->command_line_at(5).c_str()));
  } else if (input()->command_line_at(1).compare("spring") == 0) {
    parameter.push_back(atof(input()->command_line_at(3).c_str()));
    parameter.push_back(atof(input()->command_line_at(4).c_str()));
  }
  potential()->Init(parameter, atof(input()->command_line_at(2).c_str()),
                    input()->command_line_at(1));
  output()->OutputToLog(input()->command_line_at(1) +
      " potential is initialized with a cut-off radius of " +
      output()->ToString(potential()->cutoff_radius()) + " and ");
  output()->OutputToLog(false, "parameter(s)");
  for (int i = 0; i < parameter.size(); i++) {
    output()->OutputToLog(false, 
                          " " + output()->ToString(potential()->parameter(i)));
  }
  output()->OutputToLog("");
}

// CommandRun
// *****************************************************************************
// Input format:
// string []
// solver solver_parameter
// i.e., conjugate_gradient tolerance max_iteration max_displacement(optional)
// i.e., velocity_verlet timestep total_iteration
void MMM::CommandRun() {
  set_solver(input()->command_line_at(1));
  if (solver().compare("conjugate_gradient") == 0) {
    if (input()->command_line().size() == 4) {
      conjugate_gradient()->Init(0.1, atof(input()->command_line_at(2).c_str()),
                                 atoi(input()->command_line_at(3).c_str()));
    } else if (input()->command_line().size() == 5) {
      conjugate_gradient()->Init(atof(input()->command_line_at(4).c_str()),
                                 atof(input()->command_line_at(2).c_str()),
                                 atoi(input()->command_line_at(3).c_str()));
    }
    output()->OutputToLog("conjugate gradient is initialized with " +
        output()->ToString(conjugate_gradient()->max_displacement()) +
        " maximum displacement " + 
        output()->ToString(conjugate_gradient()->tolerance()) + " tolerance " + 
        output()->ToString(conjugate_gradient()->max_iteration()) + 
        " maximum iteration");
    conjugate_gradient()->Run();
  } else if (solver().compare("velocity_verlet") == 0) {
    velocity_verlet()->Init(atof(input()->command_line_at(2).c_str()),
                            atoi(input()->command_line_at(3).c_str()));
    set_is_dynamic(true);
    output()->OutputToLog("velocity verlet is initialized with " +
        output()->ToString(velocity_verlet()->timestep()) + " timestep and " +
        output()->ToString(velocity_verlet()->total_iteration()) +
                           " total iteration");
    velocity_verlet()->Run();
  }
}

// CommandSelect
// *****************************************************************************
// Input format:
// int string []
// ID selection_name selection_parameters
// i.e., ID block in dof boundary[]
// i.e., ID file filename
// i.e., ID grid in grid[] (should input a very large number for single layer
// dimensions)
// i.e., ID IDs
// i.e., ID radial in selection radius
// i.e., ID subtract selection selection
// i.e., ID surface surface depth
// i.e., ID unite selection selection
void MMM::CommandSelect() {
  int ID = atoi(input()->command_line_at(1).c_str());
  if (input()->command_line_at(2).compare("block") == 0) {
    vector<double> boundary(6, 0);
    for (int i = 3; i < input()->command_line().size(); i++) {
      boundary[i - 3] = atof(input()->command_line_at(i).c_str());
    }
    select()->SelectBlock(ID, boundary);
  } else if (input()->command_line_at(2).compare("file") == 0) {
    select()->SelectFile(ID, input()->command_line_at(3));
  } else if (input()->command_line_at(2).compare("grid") == 0) {
    int in = atoi(input()->command_line_at(3).c_str());
    vector<double> grid(3);
    grid[0] = atof(input()->command_line_at(4).c_str());
    grid[1] = atof(input()->command_line_at(5).c_str());
    grid[2] = atof(input()->command_line_at(6).c_str());
    vector<double> in_position;
    select()->SelectGrid(ID, in, grid);
  } else if (input()->command_line_at(2).compare("ID") == 0) {
    vector<string> rest_of_line(input()->command_line().begin() + 3,
                                input()->command_line().end());
    select()->SelectID(ID, rest_of_line);
  } else if (input()->command_line_at(2).compare("radial") == 0) {
    int in = atoi(input()->command_line_at(3).c_str());
    int input_selection = atoi(input()->command_line_at(4).c_str());
    double radius = atof(input()->command_line_at(5).c_str());
    select()->SelectRadial(radius, ID, in, input_selection);
  } else if (input()->command_line_at(2).compare("subtract") == 0) {
    int input_selection_1 = atoi(input()->command_line_at(3).c_str());
    int input_selection_2 = atoi(input()->command_line_at(4).c_str());
    select()->SelectSubtract(ID, input_selection_1, input_selection_2);
  } else if (input()->command_line_at(2).compare("surface") == 0) {
    string surface = input()->command_line_at(3);
    double depth = atof(input()->command_line_at(4).c_str());
    select()->SelectSurface(ID, depth, surface);
  } else if (input()->command_line_at(2).compare("unite") == 0) {
    int input_selection_1 = atoi(input()->command_line_at(3).c_str());
    int input_selection_2 = atoi(input()->command_line_at(4).c_str());
    select()->SelectUnite(ID, input_selection_1, input_selection_2);
  }
  output()->OutputToLog("selection " + output()->ToString(ID) + " of style " +
                            input()->command_line_at(2) + " is created with " +
                            output()->ToString(select()->selection_size(ID)) +
                            " atoms ");
}

// CommandTemperature
// *****************************************************************************
// Input format:
// string double double int/double[]
// thermostat initial_temperature target_temperature thermostat_parameter
// i.e., berendsen initial_temperature target_temperature
//       dissipation_coefficient
// i.e., langevin initial_temperature target_temperature
//       dissipation_coefficient
// i.e., velocity_rescaling initial_temperature target_temperature
//       thermostat_update_frequency
void MMM::CommandTemperature() {
  if (input()->command_line_at(1).compare("berendsen") == 0) {
    temperature()->InitBerendsen(atof(input()->command_line_at(4).c_str()),
                                 atof(input()->command_line_at(2).c_str()),
                                 atof(input()->command_line_at(3).c_str()));
    output()->OutputToLog("thermostat is initialized with berendsen style " +
        output()->ToString(temperature()->initial_temperature()) + 
        " initial temperature " +
        output()->ToString(temperature()->target_temperature()) +
        " target temperature " +
        output()->ToString(temperature()->dissipation_coefficient()) +
        " dissipation coefficient");
  } else if (input()->command_line_at(1).compare("langevin") == 0) {
    temperature()->InitLangevin(atof(input()->command_line_at(4).c_str()),
                                atof(input()->command_line_at(2).c_str()),
                                atof(input()->command_line_at(3).c_str()));
    output()->OutputToLog("thermostat is initialized with langevin style " +
        output()->ToString(temperature()->initial_temperature()) + 
        " initial temperature " +
        output()->ToString(temperature()->target_temperature()) +
        " target temperature " +
        output()->ToString(temperature()->dissipation_coefficient()) +
        " dissipation coefficient ");
  } else if (input()->command_line_at(1).compare("velocity_rescaling") == 0) {
    temperature()->InitVelocityRescaling(
        atof(input()->command_line_at(2).c_str()),
        atof(input()->command_line_at(3).c_str()),
        atoi(input()->command_line_at(4).c_str()));
    output()->OutputToLog("thermostat is initialized with velocity_rescaling \
        style " +
        output()->ToString(temperature()->initial_temperature()) + 
        " initial temperature " +
        output()->ToString(temperature()->target_temperature()) +
        " target temperature " +
        output()->ToString(temperature()->update_frequency()) +
        " update frequency ");
  }
}

// CommandType
// *****************************************************************************
// Input format:
// int int
// selection type
void MMM::CommandType() {
  int input_selection = atoi(input()->command_line_at(1).c_str());
  model()->SetAtomListType(select()->selection(input_selection),
                           atoi(input()->command_line_at(2).c_str()));
  string temp;
  temp = "selection " + output()->ToString(input_selection) +
         " is set to type " + input()->command_line_at(2) +
         ", atom types are ";
  for (int i = 0; i < 5; i++) {
    temp.append(
        output()->ToString(model()->atom_type_num_.get_1d_element(i)) + " ");
  }
  output()->OutputToLog(temp);
}

