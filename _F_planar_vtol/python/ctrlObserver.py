import numpy as np
import control as cnt
import VTOLParam as P
from scipy import signal

class ctrlObserver:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        tr_h = 2.0
        zeta_h = 0.707
        tr_z = 2.0
        M = 10.0
        zeta_z = 0.707
        zeta_th = 0.707
        
        integrator_pole_lon = -.15
        integrator_pole_lat = -.15
        
        # saturation limits
        self.theta_max = 10.0 * np.pi / 180.0
        self.Fe = (P.mc + 2.0 * P.mr) * P.g
        
        # PD gains for longitudinal (altitude) control
        wn_h = 2.2/tr_h
        
        # PD gains for lateral inner loop
        tr_th = tr_z/M
        wn_th = 2.2/tr_th
        
        wn_z = 2.2/tr_z
        
        self.A_lon = np.array([[0.0,1.0],[0.0,0.0]])
        self.B_lon = np.array([[0.0],[1.0/(P.mc+2.0*P.mr)]])
        self.C_lon = np.array([[1.0,0.0]])
        
        
        Cr_lon = np.array([[1.0, 0.0]])
        # form augmented system
        A1_lon = np.vstack((np.hstack((self.A_lon, np.zeros((np.size(self.A_lon,1),1)))), 
                        np.hstack((-Cr_lon, np.array([[0.0]]))) ))
        B1_lon = np.vstack( (self.B_lon, 0.0) )
        
        
        
        self.A_lat = np.array([[0.0, 0.0,               1.0,      0.0],
                      [0.0, 0.0,               0.0,      1.0],
                      [0.0, -P.Fe/(P.mc+2.0*P.mr),-P.mu/(P.mc+2.0*P.mr),0.0],
                      [0.0, 0.0,               0.0,      0.0]])
        self.B_lat = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [1.0/(P.Jc + 2.0*P.mr*P.d**2)]])
        self.C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        
        Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
        # form augmented system
        A1_lat = np.vstack((
                np.hstack((self.A_lat, np.zeros((4,1)))),
                np.hstack((-Cr_lat, np.zeros((1,1))))))
        B1_lat = np.vstack((self.B_lat, np.zeros((1,1))))
        
        
        
        # gain calculation lateral
        des_char_poly_lat = np.convolve(
            np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                        [1, 2 * zeta_z * wn_z, wn_z**2]),
                        np.poly([integrator_pole_lat])
                        )
        des_poles_lat = np.roots(des_char_poly_lat)
        
        # gain calculation longitudinal
        des_char_poly_lon = np.convolve([1, 2 * zeta_h * wn_h, wn_h**2],
                                        np.poly([integrator_pole_lon]))
        des_poles_lon = np.roots(des_char_poly_lon)
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != 3:
            print("The system is not controllable")
        else:
            K1_lon = cnt.acker(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print("The system is not controllable")
        else:
            K1_lat = cnt.acker(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4] 
        
        ### Compute Observer Gains
        # pick observer poles
        wn_h_obs = 10.0 * wn_h
        wn_z_obs = 10.0 * wn_z
        wn_th_obs = 10.0 * wn_th
        
        # Long
        des_obsv_char_poly_lon = [1, 2*zeta_h*wn_h_obs, wn_h_obs**2]
        des_obsv_poles_lon = np.roots(des_obsv_char_poly_lon)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lon.T, self.C_lon.T)) != 2:
            print("The system is not observerable")
        else:
            self.L_lon = cnt.acker(self.A_lon.T, self.C_lon.T, des_obsv_poles_lon).T
        
        # Lat
        des_obs_char_poly_lat = np.convolve(
                [1, 2 * zeta_z * wn_z_obs, wn_z_obs**2],
                [1, 2 * zeta_th * wn_th_obs, wn_th_obs**2])
        des_obs_poles_lat = np.roots(des_obs_char_poly_lat)
        # Compute the observer gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A_lat.T, self.C_lat.T)) != 4:
            print("The system is not observable")
        else:
            self.L_lat = signal.place_poles(self.A_lat.T, self.C_lat.T, 
                                        des_obs_poles_lat).gain_matrix.T
            
        # print gains to terminal
        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('L^T: ', self.L_lon.T)
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('L^T: ', self.L_lat.T)
        
        self.integrator_lon = 0.0  # integrator
        self.error_d1_lon = 0.0  # error signal delayed by 1 sample
        
        self.integrator_lat = 0.0  # integrator
        self.error_d1_lat = 0.0  # error signal delayed by 1 sample
        
        self.x_hat_lon = np.array([
            [0.0],
            [0.0]])
        self.x_hat_lat = np.array([
            [0.0],
            [0.0],
            [0.0],
            [0.0]])
        self.F_d1 = 0.0  # delayed by one sample
        self.tau_d1 = 0.0

    def update(self, reference, y):       
        z,h,theta = y[:,0]
        y_lon = np.array([[h]])
        y_lat = np.array([[z],[theta]])
        
        x_lon = self.update_observer_lon(y_lon)
        x_lat = self.update_observer_lat(y_lat)
        
        h_hat = x_lon[0][0]
        
        z_hat = x_lat[0][0]
        theta_hat = x_lat[2][0]
        
        z_r = reference[0][0]
        h_r = reference[1][0]
        
        # altitude
        error_lon = h_r - h_hat
        self.integrator_lon = self.integrator_lon \
                          + (P.Ts / 2.0) * (error_lon + self.error_d1_lon)
        self.error_d1_lon = error_lon
        
        # position
        error_lat = z_r - z_hat
        self.integrator_lat = self.integrator_lat \
                          + (P.Ts / 2.0) * (error_lat + self.error_d1_lat)
        self.error_d1_lat = error_lat
        
        
        # Compute the state feedback controller
        F_tilde = -self.K_lon @ x_lon - self.ki_lon * self.integrator_lon
        f_fl = self.Fe/np.cos(theta_hat)
        F = saturate(F_tilde[0] + f_fl, 2*P.max_thrust)
        
        tau_tilde = -self.K_lat @ x_lat - self.ki_lat * self.integrator_lat
        tau = tau_tilde[0]
        motor_thrusts = P.mixing @ np.array([[F], [tau]])
        self.F_d1 = F
        self.tau_d1 = tau
        return motor_thrusts, x_lat, x_lon

    def update_observer_lat(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.x_hat_lat, y_m)
        F2 = self.observer_f_lat(self.x_hat_lat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.x_hat_lat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.x_hat_lat + P.Ts * F3, y_m)
        self.x_hat_lat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat_lat
    
    def update_observer_lon(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lon(self.x_hat_lon, y_m)
        F2 = self.observer_f_lon(self.x_hat_lon + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lon(self.x_hat_lon + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lon(self.x_hat_lon + P.Ts * F3, y_m)
        self.x_hat_lon += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat_lon

    def observer_f_lat(self, x_hat, y_m):
        # xhatdot = A*xhat + B*u + L(y-C*xhat)    
        
        xhat_dot_lat = self.A_lat @ x_hat \
                   + self.B_lat * (self.tau_d1) \
                   + self.L_lat @ (y_m-self.C_lat @ x_hat)
        
        return xhat_dot_lat
    
    def observer_f_lon(self, x_hat, y_m):
        theta_hat = self.x_hat_lat[1][0]
        f_fl = self.Fe / np.cos(theta_hat)      
        
        xhat_dot_lon = self.A_lon @ x_hat \
                   + self.B_lon * (self.F_d1 - f_fl) \
                   + self.L_lon @ (y_m-self.C_lon @ x_hat)
        
        return xhat_dot_lon
    
    


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


