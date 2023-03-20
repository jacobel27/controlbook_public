import numpy as np
import control as cnt
import blockBeamParam as P

class ctrlStateFeedbackIntegrator:
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
        
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x

        A = np.zeros((4,4),dtype=np.float64)
        A[0,2] = 1.
        A[1,3] = 1.
        A[2,1] = -P.g
        A[3,0] = -(P.m1*P.g)/((P.m2*P.length**2/3.)+P.m1*ze**2)
        
        # A = np.array([[0.0, 0.0,               1.0,      0.0],
        #               [0.0, 0.0,               0.0,      1.0],
        #               [0.0, -P.g, 0.0, 0.0],
        #               [-(P.m1*P.g)/((P.m2*P.length**2/3.)+P.m1*ze**2), 0.0, 0.0, 0.0]])
        B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [P.length/((P.m2*P.length**2/3.)+P.m1*ze**2)]])
        C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        
        Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
        # form augmented system
        A1 = np.vstack((
                np.hstack((A, np.zeros((4,1)))),
                np.hstack((-Cr, np.zeros((1,1))))))
        B1 = np.vstack((B, np.zeros((1,1))))
        
        
        
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
        # print gains to terminal
        print('K: ', self.K)
        print('ki: ', self.ki)
        
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator_z = 0.0  # integrator
        self.error_z_d1 = 0.0  # error signal delayed by 1 sample

    def update(self, z_ref, x):
        z = x[0,0]
        
        zr_tilde = z_ref - P.z0
        error_z = zr_tilde - z
        # integrate error
        self.integrator_z = self.integrator_z \
            + (P.Ts / 2.0) * (error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        
        x_tilde = x - np.array([[P.z0],[0],[0],[0]])
        # Compute the state feedback controller
        F_tilde = -self.K @ x_tilde - self.ki * self.integrator_z
        F_fl = P.m1 * P.g * (P.z0 / P.length) + P.m2 * P.g / 2.0
        F_unsat = F_tilde + F_fl
        tau = saturate(F_unsat[0], P.Fmax)
        return tau


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


