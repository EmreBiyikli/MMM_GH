// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// time.cc
// *****************************************************************************

#include "time.h"

#include <mpi.h>

#include <cmath>

#include "mmm.h"

// Independent functions (by means of static type)
// *****************************************************************************

// ConvertTimeInSecondToHumanReadable
// *****************************************************************************
// (independent)
void Time::ConvertTimeInSecondToHumanReadable(double time_in_second,
                                              double* hour,
                                              double* minute,
                                              double* second) {
  *hour = floor(time_in_second / 60 / 60);
  *minute = floor(time_in_second / 60 - *hour * 60);
  *second = time_in_second - *hour * 60 * 60 - *minute * 60;
}

// Private
// *****************************************************************************

// ComputeParallelTimeStat
// *****************************************************************************
void Time::ComputeParallelTimeStat() {
  double temp_parallel_time;
  MPI_Status status;
  reset_parallel_time_ratio_stat();
  for (int i = 0; i < mmm()->processing_element_num(); i++) {
    // Collect the parallel time from each processing element
    if (i > 0) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (i == mmm()->processing_element()) {
        MPI_Send(parallel_time_address(), 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (mmm()->processing_element() == 0) {
        MPI_Recv(&temp_parallel_time, 1, MPI_DOUBLE, i, 10, MPI_COMM_WORLD,
                 &status);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    } else {
      temp_parallel_time = parallel_time();
    }
    // Compute parallel time ratio and update statistics
    if (mmm()->processing_element() == 0) {
      double parallel_time_ratio = temp_parallel_time / simulation_time();
      add_parallel_time_ratio_stat(1, parallel_time_ratio);
      if (parallel_time_ratio < parallel_time_ratio_stat(0)) {
        set_parallel_time_ratio_stat(0, parallel_time_ratio);
      }
      if (parallel_time_ratio > parallel_time_ratio_stat(2)) {
        set_parallel_time_ratio_stat(2, parallel_time_ratio);
      }
    }
  }
  set_parallel_time_ratio_stat(1,
      parallel_time_ratio_stat(1) / 
      static_cast<double>(mmm()->processing_element_num()));
}

