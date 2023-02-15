import numpy as np
import massParam as P

class ctrlPD:
    
    def __init__(self):
        #tuning paramters
        tr = 2.0
        zeta = 0.7
        #computre PD gains
        # open loop char polynomial and poles
        a1 = P.b/P.m
        a0 = P.k/P.m
        wn = 2.2/tr
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2
        self.kp = P.m * (alpha0 - a0)
        self.kd = P.m * (alpha1 - a1)
        print('kp: ', self.kp)
        print('kd: ', self.kd)
        
    def update(self, z_r, state):
        z = state[0][0]
        zdot = state[1][0]
        
        tau_tilde = self.kp * (z_r - z) - self.kd * zdot
        tau = saturate(tau_tilde, P.F_max)
        return tau

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u