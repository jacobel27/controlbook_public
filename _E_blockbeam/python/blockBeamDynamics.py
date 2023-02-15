import numpy as np 
import blockBeamParam as P

class blockBeamDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],             # initial block position,m
            [P.theta0],         # initial beam angle,rads
            [P.zdot0],          # initial speed of block along beam, m/s
            [P.thetadot0]       # initial angular speed of the beam,rads/s
        ]) 

        # vary physical parameters known to the controller
        self.m1 = P.m1 * (1.+alpha*(2.*np.random.rand()-1.))  # Mass of the block, kg
        self.m2 = P.m2 * (1.+alpha*(2.*np.random.rand()-1.))  # mass of beam, kg
        self.length = P.length * (1.+alpha*(2.*np.random.rand()-1.))  # length of beam, m 

        # gravity is well known
        self.g = P.g
        
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts  
        
        # set saturation limits
        self.force_limit = P.F_max

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = self.saturate(u, self.force_limit)
        
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output

        return y
 
    def f(self, state, u):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        z = state[0,0]
        theta = state[1,0]
        zdot = state[2,0]
        thetadot = state[3,0]
        F = u
        
        #The equations of motion.
        
        # M = np.array([[self.m1, 0.0],
        #              [0.0, ((self.m2*(self.length**2.0))/(3.0)) + (self.m1 * (z**2.0))]])
        
        # C = np.array([[self.m1*z*(thetadot**2.0) - (self.m1*self.g*np.sin(theta))],
        #              [(F*self.length*np.cos(theta))-(2.0*self.m1*z*zdot*thetadot)-(self.m1*self.g*z*np.cos(theta))-(((self.m2*self.g*self.length)/2.0)*np.cos(theta))]])
        
        # tmp = np.linalg.inv(M) @ C
        # zddot = tmp[0,0]
        # thetaddot = tmp[1,0]
        
        
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
        # return the output equations
        # could also use input u if needed
        z = self.state[0,0]
        theta = self.state[1,0]
        y = np.array([[z],[theta]])

        return y

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
