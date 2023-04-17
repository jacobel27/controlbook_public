import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import hummingbirdDynamics
from ctrlStateFeedbackIntegrator import ctrlStateFeedbackIntegrator

theta_ref = SignalGenerator(amplitude=np.deg2rad(15), frequency=0.03)
psi_ref = SignalGenerator(amplitude=np.deg2rad(30), frequency=0.05)

# instantiate the simulation plots and animation
controller = ctrlStateFeedbackIntegrator()
dataPlot = DataPlotter()
plt.tight_layout()
animation = HummingbirdAnimation()
hummingBird = hummingbirdDynamics(alpha=0.2)

t = P.t_start  # time starts at t_start
y = hummingBird.h()
while t < P.t_end:  # main simulation loop
    t_next_plot = t + P.t_plot
    while t < t_next_plot: 
        # phi theta psi
        ref = np.array([[0], [theta_ref.square(t)], [psi_ref.square(t)]])
        u = controller.update(ref, hummingBird.state)
        y = hummingBird.update(u)
    
        pwm_left, pwm_right = u[:,0]
        force = P.km * (pwm_left + pwm_right)
        torque = P.d * P.km * (pwm_left - pwm_right)
    
        t = t + P.Ts  # advance time by Ts
    
    animation.update(t, hummingBird.state)
    dataPlot.update(t, hummingBird.state, ref, force, torque)

    plt.pause(P.t_plot)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
