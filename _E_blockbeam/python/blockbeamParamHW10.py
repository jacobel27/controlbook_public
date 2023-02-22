# Inverted Pendulum Parameter File
import numpy as np
import blockBeamParam as P

####################################################
#       PD Control: Time Design Strategy
####################################################

####
# could also use more complicated relationship with zeta below, when not = 0.707, this would be the following:
# wn = 0.5*np.pi/(tr*np.sqrt(1-zeta**2))
# but this requires that zeta be less than 1 (underdamped system) for this equation to hold
####

# tuning parameters
tr_z = 1.2  # rise time for outer loop - first part of problem
zeta_z = 0.707  # damping ratio for outer loop
M = 10  # time scale separation between inner and outer loop
zeta_th  = 0.707  # damping ratio for inner loop
ki_z = -0.1  # integral gain on outer loop

theta_max = 30.0*np.pi/180.0  # Max theta, rads

# PD design for inner loop
ze = P.length/2.0  # equilibrium position - center of beam
b0 = P.length/(P.m2*P.length**2/3.0+P.m1*ze**2)
tr_theta = tr_z/M  # rise time for inner loop
wn_th = 2.2/tr_theta  # natural frequency for inner loop
kp_th = wn_th**2/b0  # kp - inner
kd_th = 2.0*zeta_th*wn_th/b0  # kd - inner

# DC gain for inner loop
k_DC_th = 1.0

# PD design for outer loop
wn_z = 2.2/tr_z  # natural frequency - outer loop
kp_z = -wn_z**2/P.g/k_DC_th  # kp - outer
kd_z = -2.0*zeta_z*wn_z/P.g/k_DC_th  # kd - outer

print('DC_gain', k_DC_th)
print('kp_th: ', kp_th)
print('kd_th: ', kd_th)
print('kp_z: ', kp_z)
print('ki_z: ', ki_z)
print('kd_z: ', kd_z)



