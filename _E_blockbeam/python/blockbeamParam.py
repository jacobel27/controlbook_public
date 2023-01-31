# Ball on Beam Parameter File
import numpy as np
# import control as cnt

# Physical parameters of the  ballbeam known to the controller
m1 = 0.35  # Mass of the block, kg
m2 = 2.0   # mass of beam, kg
length = 0.5   # length of beam, m
g = 9.8  # gravity at sea level, m/s^2

# parameters for animation
width = 0.05  # width of block
height = width*0.25  # height of block

# Initial Conditions
z0 = 0.0 # initial block position,m
theta0 = 0.0*np.pi/180  # initial beam angle,rads
zdot0 = 0.0  # initial speed of block along beam, m/s
thetadot0 = 0.0 # initial angular speed of the beam,rads/s

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 50.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
Fmax = 100000  # Max Force, N

# dirty derivative parameters
# sigma =   # cutoff freq for dirty derivative
# beta =  # dirty derivative gain

# equilibrium force when block is in center of beam
# ze =
# Fe =
