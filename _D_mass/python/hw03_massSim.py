import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

# instantiate reference input classes
mass = massDynamics()
reference = signalGenerator(amplitude=0.1, frequency=0.5, y_offset=0.2)
thetaRef = signalGenerator(amplitude=1.0, frequency=0.5, y_offset=0.2)
fRef = signalGenerator(amplitude=10.0, frequency=1.0)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    
    
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # set variables
        r = reference.square(t)
        theta = thetaRef.sin(t)
        f = fRef.sin(t)
        y = mass.update(f)
        t = t + P.Ts  # advance time by Ts

    # update animation
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, f)

    #plt.show()
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()