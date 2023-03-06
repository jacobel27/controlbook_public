import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import hummingbirdDynamics
from ctrlPd import ctrlPD

#fl_ref = SignalGenerator(amplitude=1.0, frequency=1.0)
#fr_ref = SignalGenerator(amplitude=1.0, frequency=2.0)

# instantiate the simulation plots and animation
controller = ctrlPD()
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingBird = hummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    u = controller.update()
    
    force = u[0,0]
    torque = u[1,0]
    
    ul = 1/(2*P.km) * (force+torque/P.d)
    ur = 1/(2*P.km) * (force-torque/P.d)
    
    y = hummingBird.update(np.array([[ul],[ur]]))
    
    
    ref = np.array([[0], [0], [P.Fe]])
    animation.update(t, hummingBird.state)
    dataPlot.update(t, hummingBird.state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
