# Hummingbird Parameter File
import numpy as np
# Initial Conditions
phi0 = 0.0 * np.pi / 180  # roll angle in rads
theta0 = 0 * np.pi / 180  # pitch angle in rads
psi0 = 0.0 * np.pi / 180  # yaw angle in rads
phidot0 = 0.0              # roll rate in rads/sec
thetadot0 = 0.0         # pitch rate in rads/sec
psidot0 = 0.0              # yaw rate in rads/sec
# Physical parameters of the hummingbird known to the controller
g = 9.81
ell1 = 0.247
ell2 = -0.039
ell3x = -0.007
ell3y = -0.007
ell3z = 0.018
ellT = 0.355
d = 0.12
m1 = 0.108862
J1x = 0.000189
J1y = 0.001953
J1z = 0.001894
m2 = 0.4717
J2x = 0.000231
J2y = 0.003274
J2z = 0.003416
m3 = 0.1905
J3x = 0.0002222
J3y = 0.0001956
J3z = 0.000027
km = g * (m1 * ell1 + m2 * ell2) / ellT  # need to find this experimentally for hardware


# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 100.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate
# saturation limits
force_max = 2.0                # Max force N
torque_max = 5.0                # Max torque, Nm

