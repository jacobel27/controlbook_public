import numpy as np
import hummingbirdParam as P


class ctrlPD:
    def __init__(self):
        self.sigma = 0.005  # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - P.Ts) \
            / (2 * self.sigma + P.Ts)  
        
        # tuning parameters
        tr_th = 1.5  # rise time
        zeta_th = 0.707  # damping ratio
        
        # PD design
        wn_th = 2.2 / tr_th
        self.kp_th = wn_th**2 / P.bTheta
        self.kd_th = 2.*zeta_th*wn_th / P.bTheta
        
        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)
        
        self.theta_dot = 0.
        self.theta_dot1 = 0.
        self.theta_d1 = 0.
        
        

    def update(self, ref, state):
        phi_r = ref[0][0]
        theta_r = ref[1][0]
        psi_r = ref[2][0]
        
        phi = state[0][0]
        theta = state[1][0]
        psi = state[2][0]
        
        
        # estimate thetadot
        self.theta_dot = self.beta * self.theta_dot1 \
            + (1 - self.beta) * ((theta - self.theta_d1))
            
        self.theta_dot = state[4][0]


        # compute the linearized torque using PD
        F_tilde = self.kp_th * (theta_r - theta) - self.kd_th * self.theta_dot
        Ffl = (P.m1*P.ell1 + P.m2*P.ell2)*(P.g/P.ellT) * np.cos(theta)
        
        force = Ffl + F_tilde
        torque = 0
        
        ul = 1/(2*P.km) * (force+torque/P.d)
        ur = 1/(2*P.km) * (force-torque/P.d)
        
        ul = saturate(ul, 0, 0.7)
        ur = saturate(ur, 0, 0.7)
        
        self.theta_d1 = theta
        self.theta_dot1 = self.theta_dot
        
        # return ul and ur (PWM)
        return np.array([[ul],[ur]])


def saturate(u, minAllowed, maxAllowed):
    if u < minAllowed:
        u = minAllowed
    elif u > maxAllowed:
        u = maxAllowed
    return u








