import numpy as np
import hummingbirdParam as P


class ctrlPID:
    def __init__(self):
        self.sigma = 0.005  # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - P.Ts) \
            / (2 * self.sigma + P.Ts)  
        
        # longitudinal dynamics:
        # tuning parameters
        tr_th = 1.0  # rise time
        zeta_th = 0.707  # damping ratio
        self.ki_th = 0.03 # integrator for longitudinal dyanamics
        
        # lateral dynamics:
        # tuning parameters
        tr_phi = 0.2
        zeta_phi = 0.707
        M = 10.0
        tr_psi = M*tr_phi
        zeta_psi = 0.707
        
        wn_phi = 2.2 / tr_phi
        wn_th = 2.2 / tr_th
        wn_psi = 2.2 / tr_psi
        
        # Longitudinal PD design
        self.kp_th = wn_th**2 / P.bTheta
        self.kd_th = 2*zeta_th*wn_th / P.bTheta
        
        # Lateral PD design
        # inner loop
        self.kp_phi = wn_phi**2 * P.J1x
        self.kd_phi = 2*zeta_phi*wn_phi * P.J1x
        self.phi_DCGain = 1.0
        # outer loop
        self.kp_psi = wn_psi**2 / (P.bPsi * self.phi_DCGain)
        self.kd_psi = 2*zeta_psi*wn_psi / (P.bPsi * self.phi_DCGain)
        self.ki_psi = 0.01
        
        print('kp_th: ', self.kp_th)
        print('ki_th: ', self.ki_th)
        print('kd_th: ', self.kd_th)
        
        print('kp_phi: ', self.kp_phi)
        print('kd_phi: ', self.kd_phi)
        
        print('kp_psi: ', self.kp_psi)
        print('ki_psiL ', self.ki_psi)
        print('kd_psi: ', self.kd_psi)
        
        self.theta_dot = 0.
        self.theta_d1 = 0.
        self.error_th_d1 = 0.0  # Error delayed by one sample 
        self.integrator_th = 0.0  # integrator 
        
        self.phi_dot = 0.
        self.phi_d1 = 0.
        
        self.psi_dot = 0.
        self.psi_d1 = 0.  
        self.error_psi_d1 = 0.0  # Error delayed by one sample 
        self.integrator_psi = 0.0  # integrator  

    def update(self, ref, state):
        phi_r, theta_r, psi_r = ref[:,0]  
        phi, theta, psi = state[:,0]    
        
        # estimate thetadot
        self.theta_dot = (
            self.beta * self.theta_dot
            + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
            )
        
        # estimate phidot
        self.phi_dot = (
            self.beta * self.phi_dot
            + (1 - self.beta) * ((phi - self.phi_d1) / P.Ts)
            )   
        
        # estimate psidot
        self.psi_dot = (
            self.beta * self.psi_dot
            + (1 - self.beta) * ((psi - self.psi_d1) / P.Ts)
            )
        
        error_th = theta_r - theta
        # integrate error
        self.integrator_th = self.integrator_th \
            + (P.Ts / 2) * (error_th + self.error_th_d1)
            
        # outer loop
        error_psi = psi_r - psi
        # integrate error
        self.integrator_psi = self.integrator_psi \
            + (P.Ts / 2) * (error_psi + self.error_psi_d1)
        
        # compute the linearized torque using PID
        F_tilde = self.kp_th * (error_th) + self.ki_th * self.integrator_th - self.kd_th * self.theta_dot
        Ffl = (P.m1*P.ell1 + P.m2*P.ell2)*(P.g/P.ellT) * np.cos(theta)
        F_unsat = Ffl + F_tilde    
        
        # lateral PD
        # outer loop gives phi_r
        phi_r = self.kp_psi * (error_psi) + self.ki_psi * self.integrator_psi - self.kd_psi * self.psi_dot
        ref[0,0] = phi_r # send back to reference array
        # inner loop
        torque_unsat = self.kp_phi * (phi_r - phi) - self.kd_phi * self.phi_dot
        
        ul_unsat = 1/(2*P.km) * (F_unsat+torque_unsat/P.d)
        ur_unsat = 1/(2*P.km) * (F_unsat-torque_unsat/P.d)
        
        ul = saturate(ul_unsat, 0, 0.7)
        ur = saturate(ur_unsat, 0, 0.7)
        
        Force = P.km * (ul + ur)
        torque = P.d * P.km * (ul - ur)
        
        # theta integrator anti - windup
        if self.ki_th != 0.0:
            self.integrator_th = self.integrator_th \
                + P.Ts / self.ki_th * (Force - F_unsat)
        
        # outer loop integrator anti - windup
        if self.ki_psi != 0.0:
            self.integrator_psi = self.integrator_psi \
                + P.Ts / self.ki_psi * (torque_unsat - torque)
        
        self.theta_d1 = theta
        self.error_th_d1 = error_th
        self.phi_d1 = phi
        self.psi_d1 = psi
        self.error_psi_d1 = error_psi
        
        # return ul and ur (PWM)
        return np.array([[ul],[ur]])


def saturate(u, minAllowed, maxAllowed):
    if u < minAllowed:
        u = minAllowed
    elif u > maxAllowed:
        u = maxAllowed
    return u