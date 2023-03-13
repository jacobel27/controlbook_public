import numpy as np
import control as cnt
import blockBeamParam as P

class ctrlStateFeedback:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        ze = P.length/2.0
        # Tuning parameters
        tr_z = 2.0
        zeta_z = 0.707
        zeta_th = 0.707
        
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
        
        
        des_char_poly = np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                                    [1, 2 * zeta_z * wn_z, wn_z**2])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = cnt.acker(A, B, des_poles)
            Cr = C[0,:]
            self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)
        # print gains to terminal
        print('K: ', self.K)
        print('kr: ', self.kr)

    def update(self, z_ref, x):
        z = x[0,0]
        
        x_tilde = x - np.array([[P.z0],[0],[0],[0]])
        zr_tilde = z_ref - P.z0
        # Compute the state feedback controller
        F_tilde = -self.K @ x_tilde + self.kr * zr_tilde
        F_fl = P.m1 * P.g * (P.z0 / P.length) + P.m2 * P.g / 2.0
        F_unsat = F_tilde + F_fl
        tau = saturate(F_unsat[0][0], P.Fmax)
        return tau


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


