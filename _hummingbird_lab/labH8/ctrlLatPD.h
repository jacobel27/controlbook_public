/**
 * \file sensors.h
 * \author Randy Beard <beard@byu.edu>
 *
 * class to implement controller
 */

#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>
#include "sensor_utilities.h"
#include "motor_utilities.h"

struct Gains {
  float kp_phi = 0.0;
  float kd_phi = 0.0;
  float kp_theta = 0.0;
  float kd_theta = 0.0;
  float kp_psi = 0.0;
  float kd_psi = 0.0;
  float ki_psi = 0.0;
  float km = 0.31;
} gains;

#include "tuning_utilities.h"

struct StateVar {
  public:
  float current = 0.0;
  float dot = 0.0;
  private:
  float d1 = 0.0;
  float d2 = 0.0;
  float d3 = 0.0;
  float dot_d1 = 0.0;
  float dot_d2 = 0.0;
  float dot_d3 = 0.0;
  public:

  void reset() {
    current = 0;
    d1 = 0;
    d2 = 0;
    d3 = 0;
    dot = 0;
    dot_d1 = 0;
    dot_d2 = 0;
    dot_d3 = 0;
  }

  // calculate the new value based on history
  void update(float newVal, float Ts) {
    d1 = newVal;
    current = 3 * d1 - 3 * d2 + d3;
    dot_d1 = (d1 - d2) / Ts;
    dot_d1 = (current - d1) / Ts;
    dot = 3 * dot_d1 - 3 * dot_d2 + dot_d3;
  }

  // cache current values to history
  void step() {
    d3 = d2;
    d2 = d1;
    d1 = current;
    dot_d3 = dot_d2;
    dot_d2 = dot_d1;
  }

  // implicit cast to float
  operator float() const { return current; }
};

// physical parameters of the system
static struct Parameters {
  float m1 = 0.108862;
  float ell1 = 0.247;
  float m2 = 0.4717;
  float ell2 = -0.039;
  float m3 = .1905;
  float g = 9.81; // force of gravity
  float ellT = 0.29;
  float ell3x = -.007;
  float ell3y = .018;
  float J1x = 0.000189;
  float J1y = 0.001953;
  float J1z = 0.001894;
  float J2x = 0.00231;
  float J2y = 0.003274;
  float J2z = 0.003416;
  float J3x = 0.0002222;
  float J3y = 0.0001956;
  float J3z = 0.000027;
  float d = 0.12;
  float fe = (m1 * ell1 + m2 * ell2) * g / ellT;
  float force_max = 0.1;
} P;

// reference structure the reference signals for psi and theta
struct Reference {
  float theta = 0.0;
  float psi = 0.0;
  float phi = 0.0;
};

// Lateral controller for hummingbird
class CtrlLatPD {
public:

  StateVar phi; // Roll
  StateVar theta; // Pitch
  StateVar psi; // Yaw

  CtrlLatPD() {
    // tuning parameters
    constexpr float tr_theta = 1.5;             // rise time for pitch
    constexpr float zeta_theta = 0.707;         // damping ratio for pitch
    constexpr float wn_theta = 2.2 / tr_theta;  // natural frequency for pitch

    // lateral dynamics:
    // tuning parameters
    constexpr float tr_phi = 0.25;
    constexpr float zeta_phi = 0.707;
    constexpr float M = 10.0;
    constexpr float tr_psi = M * tr_phi;
    constexpr float zeta_psi = 0.707;

    constexpr float wn_phi = 2.2 / tr_phi;
    constexpr float wn_psi = 2.2 / tr_psi;

    // gain calculation for longitudinal controller
    const float b_theta = P.ellT / (P.m1 * P.ell1 * P.ell1 + P.m2 * P.ell2 * P.ell2 + P.J1y + P.J2y);
    gains.kp_theta = pow(wn_theta, 2) / b_theta;
    gains.kd_theta = 2.0 * zeta_theta * wn_theta / b_theta;

    // Lateral PD design
    // inner loop
    gains.kp_phi = pow(wn_phi, 2) * P.J1x;
    gains.kd_phi = 2 * zeta_phi * wn_phi * P.J1x;
    constexpr float phi_DCGain = 1.0; // DC gain of the inner loop
    const float JT = P.m1 * pow(P.ell1, 2) + P.m2 * pow(P.ell2, 2) + P.J2z + P.m3 * (pow(P.ell3x, 2) + pow(P.ell3y, 2));
    const float bPsi = (P.ellT * P.fe) / (JT + P.J1z);
    // outer loop
    gains.kp_psi = pow(wn_psi, 2) / (bPsi * phi_DCGain);
    gains.kd_psi = 2 * zeta_psi * wn_psi / (bPsi * phi_DCGain);
  }

  void init() {
    phi.reset();
    theta.reset();
    psi.reset();
  }

  void update(const Reference &ref,
              const SensorUtilities &sensors,
              MotorUtilities &rotors,
              float Ts) {

    // tune gains
    tuneGains();

    phi.update(sensors.roll, Ts);
    theta.update(sensors.pitch, Ts);
    psi.update(sensors.yaw, Ts);

    // longitudinal control
    const float error_theta = ref.theta - theta;
    const float f_tilde = gains.kp_theta*(error_theta) - gains.kd_theta*theta.dot;

    // feedback linearized force
    const float force_fl = (P.m1 * P.ell1 + P.m2 * P.ell2) * (P.g / P.ellT) * cos(theta);

    // get inner loop ref from outer loop
    const float error_psi = ref.psi - psi;
    const float phi_ref = gains.kp_psi*(error_psi) - gains.kd_psi*psi.dot;

    // compute torque from inner loop
    const float error_phi = phi_ref - phi;
    const float torque = gains.kp_phi*(error_phi) - gains.kd_phi*phi.dot;

    const float force = force_fl + f_tilde;

    // convert force and torque to pwm and send to motors
    float left_pwm = (force + torque / P.d) / (2.0 * gains.km);
    float right_pwm = (force - torque / P.d) / (2.0 * gains.km);

    left_pwm = saturate(left_pwm, 0.0, 0.7);
    right_pwm = saturate(right_pwm, 0.0, 0.7);

    rotors.update(left_pwm, right_pwm);

    // print commanded values
    Serial.print("psi_ref:");
    Serial.print(ref.psi * 180 / PI);
    Serial.print(",");
    Serial.print("Psi:");
    Serial.println(psi.current * 180 / PI);

    phi.step();
    theta.step();
    psi.step();
  }

  float saturate(float value, float min_value, float max_value) {
    // Implements the saturation function
    return min(max(min_value, value), max_value);
  }
};

#endif
