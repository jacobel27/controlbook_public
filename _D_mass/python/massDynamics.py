import numpy as np 
import massParam as P

class massDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],      # initial position
            [P.zdot0]    # initial velocity 
        ]) 

        # vary physical parameters known to the controller
        
        # Mass of the block, kg
        self.m = P.m * (1.+alpha*(2.*np.random.rand()-1.))     
        # spring constant Kg/s^2
        self.k = P.k * (1.+alpha*(2.*np.random.rand()-1.))
        # Damping coefficient, Ns
        self.b = P.b * (1.+alpha*(2.*np.random.rand()-1.))  


        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts  
        
        # set saturation limits
        self.force_limit = P.F_max

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        F = u
        F = self.saturate(F, self.force_limit)
        
        self.rk4_step(F)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output

        return y
 
    def f(self, state, F):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        
        z = state[0][0]
        zdot = state[1][0]
        zddot = (F - self.b*zdot - self.k*z)/self.m
        
        xdot = np.array([[zdot],
                         [zddot]])
        
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        z = self.state[0,0]
        z = np.array([[z]])

        return z

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)

        return u
