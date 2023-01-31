import numpy as np
import matplotlib.pyplot as plt
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter


z_reference = signalGenerator(amplitude=0.5, frequency=0.1)
h_reference = signalGenerator(amplitude=0.5, frequency=0.1, y_offset=1.0)

zRef = signalGenerator(amplitude=4, frequency=0.1, y_offset=5)
hRef = signalGenerator(amplitude=2, frequency=0.1, y_offset=2)
thetaRef = signalGenerator(amplitude=np.pi/8, frequency=0.5, y_offset=0)

fRef = signalGenerator(amplitude=5, frequency=0.5)
tauRef = signalGenerator(amplitude=5, frequency=0.5)

animation = VTOLAnimation()
dataPlot = dataPlotter()

t = P.t_start
while t < P.t_end:
    # set variables
    z = zRef.sin(t)
    h = hRef.square(t)
    theta = thetaRef.sin(t)
    
    z_r = z_reference.sin(t)
    h_r = h_reference.square(t)
    f = fRef.sawtooth(t)
    tau = tauRef.sawtooth(t)
    motor_thrusts = P.mixing * np.array([[f], [tau]])
    
    # update animation
    state = np.array([[z],[h],[theta], [0.0], [0.0], [0.0]])
    animation.update(state, z_r)
    dataPlot.update(t, state, z_r, h_r, motor_thrusts[0], motor_thrusts[1])
    
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.1)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()