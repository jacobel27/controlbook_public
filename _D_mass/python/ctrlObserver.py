import numpy as np
import control as cnt
import massParam as P

class ctrlObserver:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        tr_obs = tr/10  # natural frequency for observer
        zeta_obs = 0.707  # damping ratio for observer
        zeta = 0.707
        integrator_pole = -5.
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.A = np.array([[0.0, 1.0],
                      [-P.k/P.m, -P.b/P.m]])
        self.B = np.array([[0.0],
                      [1./P.m]])        
        self.C = np.array([[1.0, 0.0]])
        
        # form augmented system
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A,1),1)))), 
                        np.hstack((-self.C, np.array([[0.0]]))) ))
        B1 = np.vstack( (self.B, 0.0) )
        
        # gain calculation
        wn = 2.2 / tr  # natural frequency
        des_char_poly = np.convolve([1, 2*zeta*wn, wn**2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        
        ### Observer Design ###
        
        wn_obs = 2.2 / tr_obs
        des_obsv_char_poly = [1, 2*zeta_obs*wn_obs, wn_obs**2]
        des_obsv_poles = np.roots(des_obsv_char_poly)
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 2:
            print("The system is not observerable")
        else:
            self.L = cnt.acker(self.A.T, self.C.T, des_obsv_poles).T
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', self.L.T)
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.x_hat = np.array([
            [0.0],  # z_hat_0
            [0.0],  # zdot_hat_0
        ])
        self.F_d1 = 0.0 # delayed 1 sample

    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        z_hat = x_hat[0][0]
               
        error = z_r - z_hat
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        
                
        F_unsat = -self.K @ x_hat - self.ki * self.integrator
        F = saturate(F_unsat[0], P.F_max)
        self.F_d1 = F
        return F, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y_m):
        xhat_dot = self.A @ x_hat \
                   + self.B * (self.F_d1)\
                   + self.L * (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

