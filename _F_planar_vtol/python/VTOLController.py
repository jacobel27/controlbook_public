import numpy as np
import VTOLParam as P
import VTOLParamHW10 as P10
from PIDControl import PIDControl


class VTOLController:
    '''
        This class uses multiple PID controller objects to perform
        successive loop closure and control longitudinal and lateral
        dynamics.
    '''
    def __init__(self):
        self.zCtrl = PIDControl(P10.kp_z, P10.ki_z, P10.kd_z, \
                                P.fmax, P.beta, P.Ts)
        self.hCtrl = PIDControl(P10.kp_h, P10.ki_h, P10.kd_h, \
                                P.fmax, P.beta, P.Ts)
        self.thetaCtrl = PIDControl(P10.kp_th, 0.0, P10.kd_th, \
                                    P.fmax, P.beta, P.Ts)

    def update(self, r, y):
        z_r = r[0,0]
        h_r = r[1,0]
        z = y[0,0]
        h = y[1,0]
        theta = y[2,0]
        F_tilde = self.hCtrl.PID(h_r, h, flag=False)
        F = F_tilde + P.Fe
        theta_ref = self.zCtrl.PID(z_r, z, flag=False)
        tau = self.thetaCtrl.PID(theta_ref, theta, flag=False)
        return np.array([[F], [tau]])
