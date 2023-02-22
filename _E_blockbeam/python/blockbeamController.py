import numpy as np
import blockBeamParam as P
import blockbeamParamHW10 as P10
from PIDControl import PIDControl


class blockbeamController:
    ''' 
        This class uses multiple PID controller objects to perform successive loop closure.
    '''

    def __init__(self):
        # Instantiates the SS_ctrl object
        self.zCtrl = PIDControl(P10.kp_z, P10.ki_z, P10.kd_z, P10.theta_max, P.beta, P.Ts)
        self.thetaCtrl = PIDControl(P10.kp_th, 0.0, P10.kd_th, P.Fmax, P.beta, P.Ts)
        self.limit = P.Fmax

    def update(self, z_r, y):
        z = y[0,0]
        theta = y[1,0]

        # the reference angle for theta comes from the outer loop PD control
        theta_r = self.zCtrl.PID(z_r, z, error_limit=0.1, flag=False)

        # the force applied to the cart comes from the inner loop PD control
        F_tilde = self.thetaCtrl.PID(theta_r, theta, flag=False)

        # feedback linearizing force
        F_fl = P.m1*P.g*(z/P.length) + P.m2*P.g/2.0

        # total force
        F = F_tilde + F_fl

        F_sat = self.saturate(F)

        return F_sat

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u






