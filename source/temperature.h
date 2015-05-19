// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// temperature.h
// *****************************************************************************
// Consists of Temperature class that handles temperature. In addition to the
// construction of the class, it needs to be initialized by one of the
// following functions, depending on the choice of Thermostat style:
// InitBerendsen, InitLangevin, or InitVelocityRescaling. Thermostat is applied
// by ApplyThermostat function. Functions are:
// - Set: Initial temperature of the system is set by means of velocity.
// - Compute: Temperature and kinetic energy of the system is computed.

#ifndef MMM_V14_6_TEMPERATURE_H_
#define MMM_V14_6_TEMPERATURE_H_

#include <string>
#include <vector>

#include "mediator.h"

// The class inherits from Mediator, see MANUAL.txt for the reason.
// See comment at top of file for a complete description.
class Temperature : public Mediator {
 public:
  explicit Temperature() : Mediator() {
    set_is_on(false);
  }
  ~Temperature() {}

  // This disallows copy and assign of the class.
  // TODO: activate with C++11 compiler
  // Temperature(Temperature&) = delete;
  // Temperature& operator=(const Temperature&) = delete;

  // See comment at top of file for a complete description.
  void Compute();
  // Applies thermostat.
  void ApplyThermostat(string state);
  // Initializes the class to the Berendsen thermostat style by initializing the
  // input member variables.
  void InitBerendsen(double dissipation_coefficient, double initial_temperature,
                     double target_temperature);
  // Initializes the class to the Langevin thermostat style by initializing the
  // input member variables.
  void InitLangevin(double dissipation_coefficient, double initial_temperature,
                    double target_temperature);
  // Initializes the class to the Velocity Rescaling thermostat style by
  // initializing the input member variables.
  void InitVelocityRescaling(double initial_temperature,
                             double target_temperature, int update_frequency);
  // See comment at top of file for a complete description.
  void Set();

  // Accessor and mutator functions:
  // current_temperature_
  double current_temperature() const { return current_temperature_; }
  void set_current_temperature(double current_temperature) {
    current_temperature_ = current_temperature;
  }
  // dissipation_coefficient_
  double dissipation_coefficient() const { return dissipation_coefficient_; }
  void set_dissipation_coefficient(double dissipation_coefficient) {
    dissipation_coefficient_ = dissipation_coefficient;
  }
  // initial_temperature_
  double initial_temperature() const { return initial_temperature_; }
  void set_initial_temperature(double initial_temperature) {
    initial_temperature_ = initial_temperature;
  }
  // is_on_
  bool is_on() const { return is_on_; }
  void set_is_on(bool is_on) { is_on_ = is_on; }
  // style_
  string style() const { return style_; }
  void set_style(string style) { style_ = style; }
  // target_temperature_
  double target_temperature() const { return target_temperature_; }
  void set_target_temperature(double target_temperature) {
    target_temperature_ = target_temperature;
  }
  // update_frequency_
  int update_frequency() const { return update_frequency_; }
  void set_update_frequency(int update_frequency) {
    update_frequency_ = update_frequency;
  }

 private:
  // Applies Berendsen thermostat.
  void ApplyBerendsen();
  // Applies Langevin thermostat.
  void ApplyLangevin();
  // Applies Velocity Rescaling thermostat.
  void ApplyVelocityRescaling();

  bool is_on_;                      // whether the thermostat is on
  double current_temperature_;      // current temperature
  double dissipation_coefficient_;  // dissipation coefficient
  double initial_temperature_;      // initial temperature
  double target_temperature_;       // target temperature
  int update_frequency_;            // frequency of application of thermostat
  string style_;                    // thermostat style
};

#endif  // MMM_V14_6_TEMPERATURE_H_

