import matplotlib.pyplot as plt
import numpy as np
import blockBeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics
from ctrlPD import ctrlPD

blockBeam = blockBeamDynamics()
controller = ctrlPD()

reference = signalGenerator(amplitude=0.15, frequency=0.05, y_offset=0.25)
disturbance = signalGenerator(amplitude=0.0, frequency=0.0)


dataplot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start
y = blockBeam.h()
while t < P.t_end:
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
    
        # set variables
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0
        x = blockBeam.state
        u = controller.update(r, x)
        y = blockBeam.update(u + d)
        t = t + P.Ts  # advance time by Ts
    
    # update animation
    animation.update(blockBeam.state)
    dataplot.update(t, r, blockBeam.state, u)
    
    #plt.show()
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()