// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// time.h
// *****************************************************************************
// Consists of Time class that handles the time. Functions are:
// - ComputeParallelTimeStatistics: Computes the statistics for the time spent by
// each processing element at the parallel parts of the code.
// - ConvertTimeInSecondToHumanReadable: Converts input amount of seconds into
// hour, minute, and seconds.

#ifndef MMM_V14_6_TIME_H_
#define MMM_V14_6_TIME_H_

#include <mpi.h>

#include <vector>

#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Time : public Mediator {
 public:
  explicit Time() : Mediator() {
    parallel_time_ratio_stat_.resize(3);
    reset_parallel_time();
  }
  ~Time() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Time(Time&) = delete;
  // Time& operator=(const Time&) = delete;

  // See comment at top of file for a complete description.
  static void ConvertTimeInSecondToHumanReadable(double time_in_second,
                                                 double* hour,
                                                 double* minute,
                                                 double* second);
  // See comment at top of file for a complete description.
  void ComputeParallelTimeStat();

  // Accessor and mutator functions:
  // average_parallel_time_
  double average_parallel_time() const { return average_parallel_time_; }
  void set_average_parallel_time(double average_parallel_time) {
    average_parallel_time_ = average_parallel_time;
  }
  // parallel_time_ratio_stat_
  double parallel_time_ratio_stat(int index) const {
    return parallel_time_ratio_stat_[index];
  }
  void set_parallel_time_ratio_stat(int index, double value) {
    parallel_time_ratio_stat_[index] = value;
  }
  void reset_parallel_time_ratio_stat() {
    parallel_time_ratio_stat_[0] = 1.0e6;
    parallel_time_ratio_stat_[1] = 0.0;
    parallel_time_ratio_stat_[2] = 0.0;
  }
  void add_parallel_time_ratio_stat(int index, double value) {
    parallel_time_ratio_stat_[index] += value;
  }
  // parallel_time_
  double* parallel_time_address() { return &parallel_time_; }
  double parallel_time() const { return parallel_time_; }
  void reset_parallel_time() { parallel_time_ = 0; }
  void add_parallel_time(double parallel_time) {
    parallel_time_ += parallel_time;
  }
  // simulation_time_
  double simulation_time() const { return simulation_time_; }
  void start_simulation_time() { simulation_time_ = MPI_Wtime(); }
  void finish_simulation_time() {
    simulation_time_ = MPI_Wtime() - simulation_time_;
  }

 private:
  // average of parrallel times spent by processing elements
  // in other words, it is average of parallel_time_
  double average_parallel_time_;
  // parallel time spent by current processing element
  double parallel_time_;
  double simulation_time_;                    // simulation time
  // statistics for parallel time ratio in order of minimum, average, and
  // maximum
  vector<double> parallel_time_ratio_stat_;
};

#endif  // MMM_V14_6_TIME_H_

