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

theta_ref = SignalGenerator(amplitude=np.deg2rad(15), frequency=0.05)

# instantiate the simulation plots and animation
controller = ctrlPD()
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingBird = hummingbirdDynamics()

t = P.t_start  # time starts at t_start
y = hummingBird.h()
while t < P.t_end:  # main simulation loop
    
    
    t_next_plot = t + P.t_plot
    while t < t_next_plot: 
        ref = np.array([[0], [theta_ref.square(t)], [0]])
    
        # set variables
        u = controller.update(ref, hummingBird.state)
    
        pwm_left = u[0,0]
        pwm_right = u[1,0]
    
        force = P.km * (pwm_left + pwm_right)
        torque = P.d * P.km * (pwm_left - pwm_right)
    
        y = hummingBird.update(np.array([[pwm_left],[pwm_right]]))
        t = t + P.Ts  # advance time by Ts
    
    animation.update(t, hummingBird.state)
    dataPlot.update(t, hummingBird.state, ref, force, torque)

    plt.pause(0.01)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
