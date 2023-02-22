import numpy as np
import massParam as P0
import massParamHW10 as P
from PIDControl import PIDControl

class massController:
    def __init__(self):
        # Instantiates the PD object
        self.zCtrl = PIDControl(P.kp, P.ki, P.kd, P0.F_max, P0.beta, P0.Ts)
        self.limit = P0.F_max

    def update(self, z_r, y):
        z = y[0,0]
        force_tilde = self.zCtrl.PID(z_r, z, flag=False)
        force = self.saturate(force_tilde)
        return force

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u







