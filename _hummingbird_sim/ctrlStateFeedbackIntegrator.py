import numpy as np
import control as cnt
import hummingbirdParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        tr_theta = 1.0
        zeta_theta = 0.707
        tr_psi = 1.0
        M = 10.0
        zeta_psi = 0.707
        zeta_phi = 0.707
        
        # PD gains for longitudinal (altitude) control
        wn_theta = 2.2/tr_theta
        
        # PD gains for lateral inner loop
        tr_phi = tr_psi/M
        wn_phi = 2.2/tr_phi
        
        wn_psi = 2.2/tr_psi
        
        # integrator poles
        integrator_pole_lat = -wn_psi/2.0
        integrator_pole_lon = -wn_theta/2.0
        
        A_lon = np.array([[0.0,1.0],[0.0,0.0]])
        B_lon = np.array([[0.0],[P.bTheta]])
        C_lon = np.array([[1.0, 0.0]])
        Cr_lon = np.array([C_lon[0,:]])
        
        # form augmented system
        A1_lon = np.vstack((np.hstack((A_lon, np.zeros((np.size(A_lon,1),1)))), 
                        np.hstack((-Cr_lon, np.array([[0.0]]))) ))
        B1_lon = np.vstack( (B_lon, 0.0) )
    
        
        A_lat = np.array([[0.0,                            0.0, 1.0, 0.0],
                          [0.0,                            0.0, 0.0, 1.0],
                          [0.0,                            0.0, 0.0, 0.0],
                          [((P.ellT*P.Fe)/(P.JT + P.J1z)), 0.0, 0.0, 0.0]])
        
        B_lat = np.array([[0.0],
                          [0.0],
                          [1.0/P.J1x],
                          [0.0]])
        
        C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        
        Cr_lat = np.array([C_lat[1,:]])
        # form augmented system
        A1_lat = np.vstack((
                np.hstack((A_lat, np.zeros((4,1)))),
                np.hstack((-Cr_lat, np.zeros((1,1))))))
        B1_lat = np.vstack((B_lat, np.zeros((1,1))))
        
        
        
        # gain calculation longitudinal
        des_char_poly_lon = np.convolve([1, 2 * zeta_theta * wn_theta, wn_theta**2],
                                        np.poly([integrator_pole_lon]))
        des_poles_lon = np.roots(des_char_poly_lon) 
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != 3:
            print(np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)))
            print("The longitudinal system is not controllable")
        else:
            K1_lon = cnt.acker(A1_lon, B1_lon, des_poles_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]
        
        
        # gain calculation lateral
        des_char_poly_lat = np.convolve(
            np.convolve([1, 2 * zeta_phi * wn_phi, wn_phi**2],
                        [1, 2 * zeta_psi * wn_psi, wn_psi**2]),
                        np.poly([integrator_pole_lat])
                        )
        des_poles_lat = np.roots(des_char_poly_lat) 
        
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != 5:
            print(np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)))
            print("The lateral system is not controllable")
        else:            
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
        
        phi,theta,psi,phidot,thetadot,psidot = x[:,0]
        
        x_lon = np.array([[theta],[thetadot]])
        x_lat = np.array([[phi],[psi],[phidot],[psidot]])
        
        phi_r = reference[0][0]
        theta_r = reference[1][0]
        psi_r = reference[2][0]  
        
        # altitude
        error_lon = theta_r - theta
        self.integrator_lon = self.integrator_lon \
                          + (P.Ts / 2.0) * (error_lon + self.error_d1_lon)
        self.error_d1_lon = error_lon
        
        # position
        error_lat = psi_r - psi
        self.integrator_lat = self.integrator_lat \
                          + (P.Ts / 2.0) * (error_lat + self.error_d1_lat)
        self.error_d1_lat = error_lat
        
        
        # longitudinal
        # compute the linearized torque using PID
        Ffl = (P.m1*P.ell1 + P.m2*P.ell2)*(P.g/P.ellT) * np.cos(theta)
        # OR USE P.Fe directly??? 
        
        
        # Compute the state feedback controller
        F_tilde = -self.K_lon @ x_lon - self.ki_lon*self.integrator_lon
        F = F_tilde[0] + Ffl
        
        tau_tilde = -self.K_lat @ x_lat - self.ki_lat*self.integrator_lat
        tau = tau_tilde[0]
        
        ul = 1/(2*P.km) * (F+tau/P.d)
        ur = 1/(2*P.km) * (F-tau/P.d)

        # return ul and ur (PWM)
        return np.array([[ul],[ur]])


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


