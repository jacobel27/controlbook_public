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

struct Gains
{
    float kp_phi;
    float kd_phi;
    float kp_theta;
    float ki_theta = 0.03;
    float kd_theta;
    float kp_psi;
    float kd_psi;
    float ki_psi = 0.03;
    float km = 0.30;
} gains;

// physical parameters of the system
static struct Parameters
{
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

#include "tuning_utilities.h"

inline void calculate_gains()
{
    // tuning parameters
    constexpr float tr_theta = 1.4;            // rise time for pitch
    constexpr float tr_phi = 0.25;             // rise time for roll
    
    // Longitudinal dynamics
    constexpr float zeta_theta = 0.707;        // damping ratio for pitch
    constexpr float wn_theta = 2.2 / tr_theta; // natural frequency for pitch

    // lateral dynamics:
    // tuning parameters
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

struct StateVar
{
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
    void reset()
    {
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
    void update(float newVal, float Ts)
    {
        d1 = newVal;
        current = 3 * d1 - 3 * d2 + d3;
        dot_d1 = (d1 - d2) / Ts;
        dot_d1 = (current - d1) / Ts;
        dot = 3 * dot_d1 - 3 * dot_d2 + dot_d3;
    }

    // cache current values to history
    void step()
    {
        d3 = d2;
        d2 = d1;
        d1 = current;
        dot_d3 = dot_d2;
        dot_d2 = dot_d1;
    }

    // implicit cast to float
    operator float() const { return current; }
};

struct Integrator
{
    float integrator;
    float ki;

    void reset(float ki_)
    {
        integrator = 0.0;
        ki = ki_;
        error = 0.0;
        error_d1 = 0.0;
    }

    void integrate(float inError, float Ts)
    {
        error = inError;
        integrator = integrator + (Ts / 2.0) * (error + error_d1);
    }

    void anti_windup(float satError, float Ts)
    {
        if (ki != 0.0)
        {
            integrator = integrator + Ts / ki * (satError);
        }
    }
    
    void step()
    {
        error_d1 = error;
    }

    // implicit cast to float
    operator float() const { return integrator; }

private:
    float error_d1;
    float error;
};

// reference structure the reference signals for psi and theta
struct Reference
{
    float theta = 0.0;
    float psi = 0.0;
    float phi = 0.0;
};

// Lateral controller for hummingbird
class CtrlPID
{
public:
    StateVar phi;   // Roll
    StateVar theta; // Pitch
    StateVar psi;   // Yaw

    Integrator integrator_theta;
    Integrator integrator_psi;

    CtrlPID()
    {
        calculate_gains();
    }

    void init()
    {
        phi.reset();
        theta.reset();
        psi.reset();

        integrator_theta.reset(gains.ki_theta);
        integrator_psi.reset(gains.ki_psi);
    }

    void update(const Reference &ref,
                const SensorUtilities &sensors,
                MotorUtilities &rotors,
                float Ts)
    {

        // tune gains
        tuneGains();

        phi.update(sensors.roll, Ts);
        theta.update(sensors.pitch, Ts);
        psi.update(sensors.yaw, Ts);

        // longitudinal control
        const float error_theta = ref.theta - theta;
        integrator_theta.integrate(error_theta, Ts);
        // feedback linearized force
        const float f_tilde = gains.kp_theta*(error_theta) + gains.ki_theta*integrator_theta  - gains.kd_theta*theta.dot;
        const float force_fl = (P.m1 * P.ell1 + P.m2 * P.ell2) * (P.g / P.ellT) * cos(theta);
        const float force_unsat = force_fl + f_tilde;

        // get inner loop ref from outer loop
        const float error_psi = ref.psi - psi;
        integrator_psi.integrate(error_psi, Ts);
        const float phi_ref = gains.kp_psi*(error_psi) + gains.ki_psi*integrator_psi - gains.kd_psi*psi.dot;

        const float error_phi = phi_ref - phi;
        // compute torque from inner loop
        const float torque_unsat = gains.kp_phi * (error_phi)-gains.kd_phi * phi.dot;
        

        // convert force and torque to pwm and send to motors
        float left_pwm = (force_unsat + torque_unsat / P.d) / (2.0 * gains.km);
        float right_pwm = (force_unsat - torque_unsat / P.d) / (2.0 * gains.km);

        left_pwm = saturate(left_pwm, 0.0, 0.7);
        right_pwm = saturate(right_pwm, 0.0, 0.7);

        float Force = gains.km * (left_pwm + right_pwm);
        float torque = P.d * gains.km * (left_pwm - right_pwm);

        integrator_theta.anti_windup(Force - force_unsat, Ts);
        integrator_psi.anti_windup(torque - torque_unsat, Ts);

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

        integrator_theta.step();
        integrator_psi.step();
    }

    float saturate(float value, float min_value, float max_value)
    {
        // Implements the saturation function
        return min(max(min_value, value), max_value);
    }
};

#endif
