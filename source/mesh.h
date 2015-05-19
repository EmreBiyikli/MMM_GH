// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// mesh.h
// ****************************************************************************
// (independent)
// Consists of Mesh class that carries out meshing. In addition to the
// construction of the class, it needs to be initialized by either the
// InitFile or InitQhull function, depending on the choice of meshing style.
// Meshing is performed when the BuildMesh function is called.

#ifndef MMM_V14_6_MESH_H_
#define MMM_V14_6_MESH_H_

#include <algorithm>
#include <string>
#include <vector>

using std::string;
using std::vector;

// See comment at top of file for a complete description.
class Mesh {
 public:
  Mesh();
  ~Mesh() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Mesh(Mesh&) = delete;
  // Mesh& operator=(const Mesh&) = delete;

  // Builds mesh.
  void BuildMesh();
  // Initializes the class to the file meshing style by initializing the input
  // member variables.
  void InitFile(bool is_MPI, int dimension, string input_file);
  // Initializes the class to the qhull meshing style by initializing the input
  // member variables.
  void InitQhull(bool is_MPI,
                 const vector<double>& node_position,
                 const vector<int>& node_ID,
                 int dimension);

  // Accessor and mutator functions:
  // build_count_
  int build_count() const { return build_count_; }
  void set_build_count(int build_count) { build_count_ = build_count; }
  void increase_build_count() { build_count_++; }
  // element_node_
  int* element_node_address() { return &element_node_[0]; }
  const vector<int>& element_node() const { return element_node_; }
  int element_node_at(int element, int node) const {
    return element_node_[element_node_num_ * element + node];
  }
  void set_element_node(const vector<int>& element_node) {
    element_node_ = element_node;
  }
  void push_element_node(int element_node) {
    element_node_.push_back(element_node);
  }
  void resize_element_node(int size) { element_node_.resize(size); }
  // element_node_num_
  int element_node_num() const { return element_node_num_; }
  void set_element_node_num(int element_node_num) {
    element_node_num_ = element_node_num;
  }
  // element_num_
  int* element_num_address() { return &element_num_; }
  int element_num() const { return element_num_; }
  void set_element_num(int element_num) { element_num_ = element_num; }
  void increase_element_num() { element_num_++; }
  // node_ID_
  int node_ID(int node) const { return node_ID_[node]; }
  void set_node_ID(const vector<int>& node_ID) { node_ID_ = node_ID; }
  bool is_node_ID(int atom) const {    
    if (find(node_ID_.begin(), node_ID_.end(), atom) == node_ID_.end()) {
      return false;
    } else {
      return true;
    }
  }
  // node_position_
  const vector<double>& node_position() const { return node_position_; }
  double node_position_at(int node, int dimension) const {
    return node_position_[3 * node + dimension];
  }
  void set_node_position(const vector<double>& node_position) {
    node_position_ = node_position;
  }
  // node_num_
  int node_num() const { return node_num_; }
  void set_node_num(int node_num) { node_num_ = node_num; }
  void increase_node_num() { node_num_++; }

 private:
  // Finds out if the input element is thin.
  bool IsThinElement(const vector<int>& element);
  // Makes a list of nodes of the current mesh.
  void ListNode();
  // Writes input file for Qhull.
  void OutputToQhull();
  // Performs Qhull meshing.
  void Qhull();
  // Reads mesh from input file.
  void ReadFromFile();
  // Reads mesh from Qhull output.
  void ReadFromQhull();
  // Calls Qhull.
  void RunQhull();

  // Accessor and mutator functions:
  // dimension_
  int dimension() const { return dimension_; }
  void set_dimension(int dimension) { dimension_ = dimension; }
  // input_file_
  string input_file() const { return input_file_; }
  void set_input_file(string input_file) { input_file_ = input_file; }
  // is_MPI_
  bool is_MPI() const { return is_MPI_; }
  void set_is_MPI(bool is_MPI) { is_MPI_ = is_MPI; }
  // style_
  string style() const { return style_; }
  void set_style(string style) { style_ = style; }

  bool is_MPI_;                     // whether MPI is on
  int build_count_;                 // count of mesh builds
  int dimension_;                   // dimension
  int element_node_num_;            // node number in an element
  int element_num_;                 // element number
  int node_num_;                    // node number
  string input_file_;               // input filename
  string style_;                    // meshing style
  vector<double> node_position_;    // positions of nodes (indices are local)
  // mesh by means of element nodes (holds global node IDs)
  vector<int> element_node_;
  vector<int> node_ID_;             // global node IDs
};

#endif  // MMM_V14_6_MESH_H_

