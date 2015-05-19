// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// mesh.cc
// *****************************************************************************

#include "mesh.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using std::abs;
using std::string;
using std::vector;

// Constructor
// *****************************************************************************
Mesh::Mesh() {
  set_build_count(0);
  element_num_ = 0; // MPI needs initialization
}

// Public
// *****************************************************************************

// BuildMesh
// *****************************************************************************
void Mesh::BuildMesh() {
  increase_build_count();
  if (style().compare("qhull") == 0) {
    Qhull();
  } else if (style().compare("file") == 0) {
    ReadFromFile();
  }
}

// InitQhull
// *****************************************************************************
void Mesh::InitQhull(bool is_MPI,
                     const vector<double>& node_position,
                     const vector<int>& node_ID,
                     int dimension) {
  set_is_MPI(is_MPI);
  set_node_position(node_position);
  set_node_ID(node_ID);
  set_dimension(dimension);
  set_element_node_num(dimension + 1);
  set_node_num(node_ID.size());
  set_style("qhull");
}

// InitFile
// *****************************************************************************
void Mesh::InitFile(bool is_MPI, int dimension, string input_file) {
  set_is_MPI(is_MPI);
  set_dimension(dimension);
  set_input_file(input_file);
  set_element_node_num(dimension + 1);
  set_style("file");
}

// Private
// *****************************************************************************

// IsThinElement
// *****************************************************************************
// The algorithm checks if positions of the element nodes are too close to each
// other in any dimension. This is implemented such that:
// For a dimension
//   If the distance in the current dimension for every pair of the current
//   element nodes is smaller than 1.0, then the element is thin.
// Note that input element consists of local IDs.
bool Mesh::IsThinElement(const vector<int>& element) {
  for (int i = 0; i < dimension(); i++) {
    int counter_1 = 0;
    int counter_2 = 0;
    for (int j = 0; j < element_node_num() - 1; j++) {
      int node_1 = element[j];
      for (int k = j + 1; k < element_node_num(); k++) {
        int node_2 = element[k];
        counter_1++;
        if (abs(node_position_at(node_1, i) - node_position_at(node_2, i)) <
            1.0) {
          counter_2++;
        }
      }
    }
    if (counter_1 == counter_2) {
      return true;
    }
  }
  return false;
}

// ListNode
// *****************************************************************************
void Mesh::ListNode() {
  vector<int> node_list = element_node();
  sort(node_list.begin(), node_list.end());
  node_list.erase(unique(node_list.begin(), node_list.end()), node_list.end());
  set_node_ID(node_list);
  set_node_num(node_list.size());
}

// OutputToQhull
// *****************************************************************************
// Outputs in the format:
// dimension
// node_number
// node_1_x node_1_y node_1_z
// node_2_x node_2_y node_2_z
// ...
// Dimensions y & z or z are omitted in 1-D and 2-D cases, respectively.
void Mesh::OutputToQhull() {
  FILE* to_Qhull_ptr = fopen("./object/to_Qhull.txt", "w");
  if (to_Qhull_ptr == NULL) {
    printf("cannot open ./object/to_Qhull.txt\n");
    exit(1);
  }
  fprintf(to_Qhull_ptr, "%d\n%d\n", dimension(), node_num());  
  for (int i = 0; i < node_num(); i++) {
    for (int j = 0; j < dimension(); j++) {
      fprintf(to_Qhull_ptr, "%.8f ", node_position_at(i, j));      
    }
    fprintf(to_Qhull_ptr, "\n");
  }
  fclose(to_Qhull_ptr);  
}

// Qhull
// *****************************************************************************
void Mesh::Qhull() {
  int processing_element;
  MPI_Comm_rank(MPI_COMM_WORLD, &processing_element);
  if (processing_element == 0) {
    OutputToQhull();
    RunQhull();
    ReadFromQhull();
    ListNode();  // Because removing elements may remove some nodes.
  }
  if (is_MPI()) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(element_num_address(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    resize_element_node(element_num() * element_node_num());
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(element_node_address(), element_num() * element_node_num(), 
              MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// ReadFromFile
// *****************************************************************************
// Reads in the format:
// element_number
// element_1_node_1 element_1_node_2 ...
// ...
// The line continues until element_node_number and the file until
// element_number.
// Reading is performed by processing element 0 and then broadcasted to other
// processing elements.
void Mesh::ReadFromFile() {
  int processing_element;
  FILE* input_file_ptr = fopen(input_file().c_str(), "r");
  if (input_file_ptr == NULL) {
    printf("cannot open the input mesh file: %s\n", input_file().c_str());
    exit(1);
  }
  int temp_element_num;
  fscanf(input_file_ptr, "%d", &temp_element_num);
  set_element_num(temp_element_num);
  for (int i = 0; i < element_num(); i++) {
    for (int j = 0; j < element_node_num(); j++) {
      int temp_element_node;
      fscanf(input_file_ptr, "%d", &temp_element_node);
      push_element_node(temp_element_node);
    }
  }
  fclose(input_file_ptr);
  ListNode();
}

// ReadFromQhull
// *****************************************************************************
// Reads in the format:
// element_number
// element_1_node_1 element_1_node_2 ...
// ...
// The line continues until element_node_number and the file until
// element_number.
// Reading is performed by processing element 0 and broadcasted to other
// processing elements later at an outer level function.
void Mesh::ReadFromQhull() {
  FILE* from_Qhull_ptr = fopen("./object/from_Qhull.txt", "r");
  if (from_Qhull_ptr == NULL) {
    printf("cannot open ./object/from_Qhull.txt\n");
    exit(1);
  }
  int temp_element_num;
  fscanf(from_Qhull_ptr, "%d", &temp_element_num);
  set_element_num(0);
  for (int i = 0; i < temp_element_num; i++) {
    vector<int> temp_element;
    for (int j = 0; j < element_node_num(); j++) {
      int temp_element_node;
      fscanf(from_Qhull_ptr, "%d", &temp_element_node);
      temp_element.push_back(temp_element_node);
    }
    if (IsThinElement(temp_element)) continue;
    increase_element_num();
    for (int j = 0; j < element_node_num(); j++) {
      push_element_node(node_ID(temp_element[j]));
    }
  }
  fclose(from_Qhull_ptr);
}

// RunQhull
// *****************************************************************************
void Mesh::RunQhull() {
  string command = "qdelaunay QJ Pp i < ./object/to_Qhull.txt > ";
  command += "./object/from_Qhull.txt";
  int return_val = system(command.c_str());
  if (return_val != 0) {
    printf("error in mesh system command\n");
    exit(1);
  }
}

