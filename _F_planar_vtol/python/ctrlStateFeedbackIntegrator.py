import numpy as np
import control as cnt
import VTOLParam as P

class ctrlStateFeedbackIntegrator:
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
        
        integrator_pole_lon = -5.
        integrator_pole_lat = -5.
        
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
        
        
        Cr_lon = np.array([[1.0, 0.0]])
        # form augmented system
        A1_lon = np.vstack((np.hstack((A_lon, np.zeros((np.size(A_lon,1),1)))), 
                        np.hstack((-Cr_lon, np.array([[0.0]]))) ))
        B1_lon = np.vstack( (B_lon, 0.0) )
        
        
        
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
        
        Cr_lat = np.array([[1.0, 0.0, 0.0, 0.0]])
        # form augmented system
        A1_lat = np.vstack((
                np.hstack((A_lat, np.zeros((4,1)))),
                np.hstack((-Cr_lat, np.zeros((1,1))))))
        B1_lat = np.vstack((B_lat, np.zeros((1,1))))
        
        
        
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
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print("The system is not controllable")
        else:
            K1_lon = cnt.acker(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]
            
            K1_lat = cnt.acker(A1_lat, B1_lat, des_poles_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4]
            
        # print gains to terminal
        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        
        self.integrator_lon = 0.0  # integrator
        self.error_d1_lon = 0.0  # error signal delayed by 1 sample
        
        self.integrator_lat = 0.0  # integrator
        self.error_d1_lat = 0.0  # error signal delayed by 1 sample

    def update(self, reference, x):
        z,h,theta,zdot,hdot,thetadot = x[:,0]
        
        x_lon = np.array([[h],[hdot]])
        x_lat = np.array([[z],[theta],[zdot],[thetadot]])
        
        z_r = reference[0][0]
        h_r = reference[1][0]
        
        # altitude
        error_lon = h_r - h
        self.integrator_lon = self.integrator_lon \
                          + (P.Ts / 2.0) * (error_lon + self.error_d1_lon)
        self.error_d1_lon = error_lon
        
        # position
        error_lat = z_r - z
        self.integrator_lat = self.integrator_lat \
                          + (P.Ts / 2.0) * (error_lat + self.error_d1_lat)
        self.error_d1_lat = error_lat
        
        
        # Compute the state feedback controller
        F_tilde = -self.K_lon @ x_lon - self.ki_lon * self.integrator_lon
        F = saturate(F_tilde[0] + self.Fe/np.cos(theta), 2*P.max_thrust)
        
        tau_tilde = -self.K_lat @ x_lat - self.ki_lat * self.integrator_lat
        tau = tau_tilde[0]
        motor_thrusts = P.mixing @ np.array([[F], [tau]])
        return motor_thrusts


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


