# VTOL Parameter File
import numpy as np
import VTOLParam as P

# tuning parameters
tr_h = 3.0  # rise time for altitude - original (3.0)
zeta_h = 0.95  # damping ratio for altitude - original (0.707)
tr_z = 3.0  # rise time for outer lateral loop (position) - original
M = 10.0  # time separation between inner and outer lateral loops
zeta_z = 0.9  # damping ratio for outer lateral loop
zeta_th = 0.9  # damping ratio for inner lateral loop
ki_h = 1.0  # integrator on altitude

# for ki_z, it actually works best if ki_z = 0.0!! This tells us something
# about the system type. But it is possible to get reasonable performance
# if we increase zeta for theta and z as shown below. The negative is
# strange, but if we examine kp_z and kd_z, they are also negative. This
# is because the plant itself has a negative gain. Making kz_i positive
# actually drives the system to be unstable.
ki_z = 0.0 #-0.02  # integrator on position

# PD gains for longitudinal (altitude) control
#wn_h = 2.2 / tr_h  # original equation for natural frequency (before changing zeta from 0.707)
wn_h = 0.5*np.pi/(tr_h*np.sqrt(1-zeta_h**2))
Delta_cl_d = [1, 2 * zeta_h * wn_h, wn_h ** 2.0]  # desired closed loop char eq
kp_h = Delta_cl_d[2] * (P.mc + 2.0 * P.mr)  # kp - altitude
kd_h = Delta_cl_d[1] * (P.mc + 2.0 * P.mr)  # kd = altitude
Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force

# PD gains for lateral inner loop
b0 = 1.0 / (P.Jc + 2.0 * P.mr * P.d ** 2)
tr_th = tr_z / M
#wn_th = 2.2 / tr_th  # original equation for natural frequency (before changing zeta from 0.707)
wn_th = 0.5*np.pi/(tr_th*np.sqrt(1-zeta_th**2))
kp_th = wn_th ** 2.0 / b0
kd_th = 2.0 * zeta_th * wn_th / b0

# PD gain for lateral outer loop
b1 = -P.Fe / (P.mc + 2.0 * P.mr)
a1 = P.mu / (P.mc + 2.0 * P.mr)
#wn_z = 2.2 / tr_z  # original equation for natural frequency (before changing zeta from 0.707)
wn_z = 0.5*np.pi/(tr_z*np.sqrt(1-zeta_z**2))
kp_z = wn_z ** 2.0 / b1
kd_z = (2.0 * zeta_z * wn_z - a1) / b1

print('kp_z: ', kp_z)
print('ki_z: ', ki_z)
print('kd_z: ', kd_z)
print('kp_h: ', kp_h)
print('ki_z: ', ki_h)
print('kd_h: ', kd_h)
print('kp_th: ', kp_th)
print('kd_th: ', kd_th)
