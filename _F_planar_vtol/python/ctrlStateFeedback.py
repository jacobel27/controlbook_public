import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedback:
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
        
        # saturation limits
        self.theta_max = 10.0 * np.pi / 180.0
        self.Fe = (P.mc + 2.0 * P.mr) * P.g
        
        # PD gains for longitudinal (altitude) control
        wn_h = 2.2/tr_h
        
        # PD gains for lateral inner loop
        tr_th = tr_z/M
        wn_th = 2.2/tr_th
        
        wn_z = 2.2/tr_z
        
        A_lon = np.array([[0.0,1.0],[0.0,0.0]])
        B_lon = np.array([[0.0],[1.0/(P.mc+2.0*P.mr)]])
        C_lon = np.array([[1.0,0.0]])
        
        
        A_lat = np.array([[0.0, 0.0,               1.0,      0.0],
                      [0.0, 0.0,               0.0,      1.0],
                      [0.0, -P.Fe/(P.mc+2.0*P.mr),-P.mu/(P.mc+2.0*P.mr),0.0],
                      [0.0, 0.0,               0.0,      0.0]])
        B_lat = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [1.0/(P.Jc + 2.0*P.mr*P.d**2)]])
        C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        # gain calculation
        des_char_poly_lat = np.convolve([1, 2 * zeta_th * wn_th, wn_th**2],
                                    [1, 2 * zeta_z * wn_z, wn_z**2])
        des_poles_lat = np.roots(des_char_poly_lat)
        
        des_char_poly_lon = [1, 2 * zeta_h * wn_h, wn_h**2]
        des_poles_lon = np.roots(des_char_poly_lon)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A_lat, B_lat)) != 4:
            print("The system is not controllable")
        else:
            self.K_lon = cnt.acker(A_lon, B_lon, des_poles_lon)
            Cr_lon = C_lon[0,:]
            self.kr_lon = -1.0 / (Cr_lon @ np.linalg.inv(A_lon - B_lon @ self.K_lon) @ B_lon)
            
            self.K_lat = cnt.acker(A_lat, B_lat, des_poles_lat)
            Cr_lat = C_lat[0,:]
            self.kr_lat = -1.0 / (Cr_lat @ np.linalg.inv(A_lat - B_lat @ self.K_lat) @ B_lat)
        # print gains to terminal
        print('K_lon: ', self.K_lon)
        print('kr_lon: ', self.kr_lon)
        print('K_lat: ', self.K_lat)
        print('kr_lat: ', self.kr_lat)

    def update(self, reference, x):
        z,h,theta,zdot,hdot,thetadot = x[:,0]
        
        x_lon = np.array([[h],[hdot]])
        x_lat = np.array([[z],[theta],[zdot],[thetadot]])
        
        z_r = reference[0][0]
        h_r = reference[1][0]
        # Compute the state feedback controller
        F_tilde = -self.K_lon @ x_lon + self.kr_lon * h_r
        F = saturate(F_tilde[0][0] + self.Fe, 2*P.max_thrust)
        tau_tilde = -self.K_lat @ x_lat + self.kr_lat * z_r
        tau = saturate(tau_tilde[0][0], 2*P.max_thrust*P.d)
        motor_thrusts = P.mixing @ np.array([[F], [tau]])
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


