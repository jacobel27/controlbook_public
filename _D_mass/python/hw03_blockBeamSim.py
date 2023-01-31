import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics

beam = blockBeamDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
zRef = signalGenerator(amplitude=0.05, frequency=0.5, y_offset=0.2)
thetaRef = signalGenerator(amplitude=np.pi/8, frequency=0.1, y_offset=0)
fRef = signalGenerator(amplitude=5.0, frequency=0.5)


animation = blockbeamAnimation()
dataplot = dataPlotter()

t = P.t_start
while t < P.t_end:
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
    
        # set variables
        z = zRef.sin(t) # Horizontal position of cart, m
        theta = thetaRef.square(t) # Angle of block-beam, rads
    
        r = reference.square(t)
        f = fRef.sawtooth(t)
        y = beam.update(f)
        t = t + P.Ts  # advance time by Ts
    
    # update animation
    animation.update(beam.state)
    dataplot.update(t, r, beam.state, f)
    
    #plt.show()
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()