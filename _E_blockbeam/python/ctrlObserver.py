import numpy as np
import control as cnt
import blockBeamParam as P
from scipy import signal

class ctrlObserver:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        ze = P.length/2.0
        # Tuning parameters
        tr_z = 2.0
        zeta_z = 0.707
        zeta_th = 0.707
        integrator_pole = -5.0
        
        # Inner Loop
        M = 10.0
        tr_theta = tr_z/M
        wn_th = 2.2/tr_theta
        
        # Outer Loop
        wn_z = 2.2/tr_z
        
        # Observer
        tr_z_obs = tr_z/5.0 # rise time for position
        tr_theta_obs = tr_theta / 5.0  # rise time for angle
        wn_z_obs = 2.2 / tr_z_obs
        wn_th_obs = 2.2 / tr_theta_obs
        
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x

        self.A = np.zeros((4,4),dtype=np.float64)
        self.A[0,2] = 1.
        self.A[1,3] = 1.
        self.A[2,1] = -P.g
        self.A[3,0] = -(P.m1*P.g)/((P.m2*P.length**2/3.)+P.m1*ze**2)
        
        # A = np.array([[0.0, 0.0,               1.0,      0.0],
        #               [0.0, 0.0,               0.0,      1.0],
        #               [0.0, -P.g, 0.0, 0.0],
        #               [-(P.m1*P.g)/((P.m2*P.length**2/3.)+P.m1*ze**2), 0.0, 0.0, 0.0]])
        self.B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/((P.m2*P.length**2/3.)+P.m1*ze**2)]])
        self.C = np.array([[1.0, 0.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0, 0.0]])
        
        Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
        # form augmented system
        A1 = np.vstack((
                np.hstack((self.A, np.zeros((4,1)))),
                np.hstack((-Cr, np.zeros((1,1))))))
        B1 = np.vstack((self.B, np.zeros((1,1))))
        
        des_char_poly = np.convolve(
            np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                        [1, 2 * zeta_z * wn_z, wn_z**2]),
            np.poly([integrator_pole])
        )
        
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.acker(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
            
        ### Compute Observer Gains ###
        des_obs_char_poly = np.convolve(
                [1, 2 * zeta_z * wn_z_obs, wn_z_obs**2],
                [1, 2 * zeta_th * wn_th_obs, wn_th_obs**2])
        des_obs_poles = np.roots(des_obs_char_poly)
        # Compute the observer gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print("The system is not observable")
        else:
            self.L = signal.place_poles(self.A.T, self.C.T, 
                                        des_obs_poles).gain_matrix.T
        # print gains to terminal
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('L^T: ', self.L.T)    
        
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample
        
        # estimated state variables
        self.x_hat = np.array([
            [0.0],
            [0.0],
            [0.0],
            [0.0]])
        self.tau_d1 = 0.0  # delayed by one sample

    def update(self, z_ref, y):
        x_hat = self.update_observer(y)
        z_hat = x_hat[0][0]
        
        
        zr_tilde = z_ref - P.z0
        error_z = zr_tilde - z_hat
        # integrate error
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        
        x_tilde = x_hat - np.array([[P.z0],[0],[0],[0]])
        # Compute the state feedback controller
        F_tilde = -self.K @ x_tilde - self.ki * self.integrator_z
        F_fl = P.m1 * P.g * (P.z0 / P.length) + P.m2 * P.g / 2.0
        F_unsat = F_tilde + F_fl
        tau = saturate(F_unsat[0], P.Fmax)
        self.tau_d1 = tau
        return tau, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)
        tau_fl = P.m1 * P.g * (P.z0 / P.length) + P.m2 * P.g / 2.0
        xhat_dot = self.A @ x_hat \
                   + self.B * (self.tau_d1 - tau_fl) \
                   + self.L @ (y_m-self.C @ x_hat)
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


