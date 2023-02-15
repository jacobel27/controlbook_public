import numpy as np
import matplotlib.pyplot as plt
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlPD import ctrlPD

vtol = VTOLDynamics()
controller = ctrlPD()
z_reference = signalGenerator(amplitude=2.5, frequency=0.08, y_offset=3.0)
h_reference = signalGenerator(amplitude=3.0, frequency=0.03, y_offset=5.0)


animation = VTOLAnimation()
dataPlot = dataPlotter()

t = P.t_start
y = vtol.h()
while t < P.t_end:
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        z_ref = z_reference.square(t)
        h_ref = h_reference.square(t)
        r = np.array([[z_ref], [h_ref]])
        d = np.array([[0.0], [0.0]])
        n = 0.0
        x = vtol.state
        u = controller.update(r, x)
        y = vtol.update(u + d)
        t = t + P.Ts  # advance time by Ts

    
    # update animation
    animation.update(vtol.state, z_ref)
    dataPlot.update(t, vtol.state, z_ref, h_ref, u[0][0], u[1][0])
    
    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()