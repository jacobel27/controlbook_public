import numpy as np
import blockBeamParam as P

class blockbeamDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],  # z initial block position
            [P.theta0],  # Theta initial beam angle
            [P.zdot0],  # zdot initial block velocity
            [P.thetadot0],  # Thetadot initial beam angular velocity
        ])

        #################################################
        # The parameters for any physical system are never known exactly. Feedback
        # systems need to be designed to be robust to this uncertainty. In the simulation
        # we model uncertainty by changing the physical parameters by a uniform random variable
        # that represents alpha*100 % of the parameter, i.e., alpha = 0.2, means that the parameter
        # may change by up to 20%. A different parameter value is chosen every time the simulation
        # is run. This solution does not require the "alpha" parameter to be defined unless we want
        # to model uncertainty in our model. This is something that comes later in the book when
        # doing feedback control.
        #################################################

        self.Ts = P.Ts
        self.m1 = P.m1 * (1+2*alpha*np.random.rand()-alpha)
        self.m2 = P.m2 * (1+2*alpha*np.random.rand()-alpha)
        self.length = P.length * (1+2*alpha*np.random.rand()-alpha)
        self.g = P.g

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        self.rk4_step(u)  # propagate the state by one time sample

        # separating out "y" by itself is currently not required, but will be in future homework
        y = self.h()  # using a measurement model, return the corresponding output
        return y

    def f(self, state, u):
        # return xdot = f(x,u)
        z = state[0,0]
        theta = state[1,0]
        zdot = state[2,0]
        thetadot = state[3,0]
        F = u
        # The equations of motion.
        zddot = (1.0/self.m1)*(self.m1*z*thetadot**2
                               - self.m1*self.g*np.sin(theta))

        thetaddot = (1.0/((self.m2*self.length**2)/3.0
                          + self.m1*z**2))*(-2.0*self.m1*z*zdot*thetadot
                                            - self.m1*self.g*z*np.cos(theta)
                    - self.m2*self.g*self.length/2.0*np.cos(theta)
                    + self.length*F*np.cos(theta))
        # build xdot and return
        xdot = np.array([[zdot], [thetadot], [zddot], [thetaddot]])
        return xdot

    def h(self):
        # return y = h(x)
        z = self.state[0,0]
        theta = self.state[1,0]
        y = np.array([[z], [theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
