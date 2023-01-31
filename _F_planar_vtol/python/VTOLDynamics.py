import numpy as np 
import VTOLParam as P

class VTOLDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],        # initial lateral position
            [P.h0],        # initial altitude
            [P.theta0],    # initial roll angle
            [P.zdot0],     # initial lateral velocity
            [P.hdot0],     # initial climb rate
            [P.thetadot0]  # initial roll rate
        ])

        # simulation time step
        self.Ts = P.Ts

        self.mc = P.mc * (1.+alpha*(2.*np.random.rand()-1.)) # kg
        self.mr = P.mr * (1.+alpha*(2.*np.random.rand()-1.)) # kg
        self.Jc = P.Jc * (1.+alpha*(2.*np.random.rand()-1.)) # kg m^2
        self.d = P.d * (1.+alpha*(2.*np.random.rand()-1.))  # m
        self.mu = P.mu * (1.+alpha*(2.*np.random.rand()-1.)) # kg/s
        self.F_wind = P.F_wind * (1.+alpha*(2.*np.random.rand()-1.))

        # gravity constant is well known, don't change.
        self.g = P.g
        self.force_limit = P.fmax

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = self.saturate(u, self.force_limit)

        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output

        return y

    def f(self, state, u):
        # Return xdot = f(x,u)
        z = state[0,0]
        h = state[1,0]
        theta = state[2,0]
        zdot = state[3,0]
        hdot = state[4,0]
        thetadot = state[5,0]
        fr = u[0,0]
        fl = u[1,0]

        # The equations of motion.
        # M = np.array([[(self.mc + (2.0*self.mr)),0,0],
        #               [0,(self.mc + (2.0*self.mr)),0],
        #               [0,0,(self.Jc + (2.0*self.mr*(self.d**2)))]])
        # C = np.array([[-(fr + fl)*np.sin(theta) - self.mu*zdot],
        #               [-((self.mc + (2.0*self.mr)*self.g) + ((fr + fl)*np.cos(theta)))],
        #               [self.d*(fr-fl)]])
        # tmp = np.linalg.inv(M) @ C
        # zddot = tmp[0,0]
        # hddot = tmp[1,0]
        # thetaddot = tmp[2,0]
        
        zddot = (-(fr + fl)*np.sin(theta) - self.mu*zdot + self.F_wind) / (self.mc + (2.0*self.mr))
        hddot = (-((self.mc + (2.0*self.mr)*self.g) + ((fr + fl)*np.cos(theta)))) / (self.mc + (2.0*self.mr))
        thetaddot = self.d*(fr-fl) / (self.Jc + (2.0*self.mr*(self.d**2)))

        # build xdot and return
        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])

        return xdot

    def h(self):
        # return y = h(x)
        z = self.state[0,0]
        h = self.state[1,0]
        theta = self.state[2,0]
        y = np.array([[z],[h],[theta]])

        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

    def saturate(self, u, limit):
        if abs(u[0,0]) > limit:
            u[0,0] = limit*np.sign(u[0,0])
        if abs(u[1,0]) > limit:
            u[1,0] = limit*np.sign(u[1,0])

        return u
