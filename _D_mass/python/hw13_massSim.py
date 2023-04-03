import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlObserver import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

# instantiate reference input classes
mass = massDynamics(alpha=0.0)
controller = ctrlObserver()
reference = signalGenerator(amplitude=1.0, frequency=0.04)
disturbance = signalGenerator(amplitude=0.25)

dataPlot = dataPlotter()
animation = massAnimation()
dataPlotObserver = dataPlotterObserver()

t = P.t_start  # time starts at t_start
y = mass.h()
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        # set variables
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0
        u, xhat = controller.update(r, y+n)
        y = mass.update(u + d)
        t = t + P.Ts  # advance time by Ts

    # update animation
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    dataPlotObserver.update(t, mass.state, xhat)

    plt.pause(P.t_plot)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()