// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// matrix.h
// ****************************************************************************
// Consists of Matrix class that defines a type specialized to handle 2-D data.
// The type is defined as a template.

#ifndef MMM_V14_6_MATRIX_H_
#define MMM_V14_6_MATRIX_H_

#include <vector>

using std::vector;

// See comment at top of file for a complete description.
template <class T>
class Matrix {
 public:
  Matrix() {}
  ~Matrix() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Matrix(Matrix&) = delete;
  // Matrix& operator=(const Matrix&) = delete;

  // Accessor and mutator functions:
  // data get
  T* address() { return &data_[0]; }
  const vector<T>& get() const { return data_; }
  T get_element(int row, int column) const {
    return data_[column_size_ * row + column];
  }
  T get_1d_element(int row) const { return data_[row]; }
  // data set
  void set(const vector<T>& value) { data_ = value; }
  void set_element(int row, int column, T value) {
    data_[column_size_ * row + column] = value;
  }
  void set_1d_element(int row, T value) { data_[row] = value; }
  void set_all_to(T value) { data_.assign(data_.size(), value); }
  void reset() { data_.assign(data_.size(), 0); }
  void add(const vector<T>& value) {
    for (int i = 0; i < data_.size(); i++) { data_[i] += value[i]; }
  }
  void add_element(int row, int column, T value) {
    data_[column_size_ * row + column] += value;
  }
  void add_1d_element(int row, T value) { data_[row] += value; }
  void multiply(int row, int column, double value) {
    data_[column_size_ * row + column] *= value;
  }
  // size set
  void set_size(int row_size, int column_size) {
    data_.resize(row_size * column_size);
    column_size_ = column_size;
  }
  void set_1d_size(int row_size) {
    data_.resize(row_size);
    column_size_  = 1;
  }

 private:
  int column_size_;   // column size
  vector<T> data_;    // data
};

#endif  // MMM_V14_6_MATRIX_H_

