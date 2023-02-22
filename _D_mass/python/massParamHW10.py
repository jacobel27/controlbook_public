# mass-spring-damper Parameter File
import numpy as np
import massParam as P

Ts = P.Ts  # sample rate of the controller
beta = P.beta  # dirty derivative gain
tau_max = P.F_max  # limit on control signal

#  tuning parameters
tr = 5

# tuned value of the rise time
tr = 2.0
zeta = 0.707
ki = 1.5  # integrator gain


# compute PD gains
# open loop char polynomial and poles
a1 = P.b/P.m
a0 = P.k/P.m

####
# could also use more complicated relationship with zeta below, when not = 0.707, this would be the following:
# wn = 0.5*np.pi/(tr*np.sqrt(1-zeta**2))
# but this requires that zeta be less than 1 (underdamped system) for this equation to hold
####
wn = 2.2/tr
alpha1 = 2.0*zeta*wn
alpha0 = wn**2

kp = P.m*(alpha0-a0)
kd = P.m*(alpha1-a1)

print('kp: ', kp)
print('ki: ', ki)
print('kd: ', kd)


