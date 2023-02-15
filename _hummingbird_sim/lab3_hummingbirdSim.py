import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import hummingbirdDynamics

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)

fl_ref = SignalGenerator(amplitude=1.0, frequency=1.0)
fr_ref = SignalGenerator(amplitude=1.0, frequency=2.0)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingBird = hummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    fl = fl_ref.sin(t)
    fr = fr_ref.square(t)
    
    y = hummingBird.update(np.array([[fl],[fr]]))
    
    force = P.km * (fl + fr)
    torque = hummingBird.d * P.km * (fl - fr)
    
    ref = np.array([[0], [0], [0]])
    animation.update(t, hummingBird.state)
    dataPlot.update(t, hummingBird.state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
