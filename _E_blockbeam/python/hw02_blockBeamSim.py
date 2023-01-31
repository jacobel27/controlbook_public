import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter

thetaRef = signalGenerator(amplitude=np.pi/8, frequency=0.1, y_offset=0)
zRef = signalGenerator(amplitude=0.05, frequency=0.5, y_offset=0.2)

reference = signalGenerator(amplitude=0.5, frequency=0.1)
fRef = signalGenerator(amplitude=5, frequency=0.5)

animation = blockbeamAnimation()
dataplot = dataPlotter()

t = P.t_start
while t < P.t_end:
    # set variables
    z = zRef.sin(t) # Horizontal position of cart, m
    theta = thetaRef.square(t) # Angle of block-beam, rads
    
    r = reference.square(t)
    f = fRef.sawtooth(t)
    
    # update animation
    state = np.array([[z],[theta], [0.0], [0.0]])
    animation.update(state)
    dataplot.update(t, r, state, f)
    
    #plt.show()
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()