# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython


from libc.math cimport sqrt # sin, cos, acos, exp, sqrt, fabs, M_PI

cpdef C_torque_limits(double R, double V, double n2, double n1, double m, double b, double h, double rho,
                      double cd, double a, double mu):
    kn = n2 / n1 * R
    to_max_wh = kn * (m * 9.81 * b / h - V ** 2 * rho * cd * a / 2)
    to_max_fr = mu * m * 9.81 * kn
    to_max = min(to_max_wh, to_max_fr)
    return [to_max, to_max_wh, to_max_fr]


cpdef inline chain_eff_single(int n2, int n1, double p, double m_ch, double cd, double w_d, double torque, double mu_p, double r_b):
    # returns chain efficiency = e, chain speed = v.
    # N1 = driven sprocket, N2 = driving.
    # m_ch = chain mass per unit length, kg/m
    # p = pitch, mm
    # Cd = Sprocket center distance, m
    # w_d = Angular speed of N1, rad s^-1
    # torque = driving sprocket torque, Nm
    # mu_p = 0.11;            % Pin-bush friction co-efficient, from paper
    # r_b = 0.00508;          % Internal diameter of chain bush
    cdef double pi,den,p_loss,e,r_d,f_c,f_cf,f_chaintension #

    pi = 3.14159265

    w_d = abs(w_d)
    torque = abs(torque)

    r_d = n2 * p / 2000.0 / pi
    # r_o = r_d * n1 / n2
    # w_o = w_d*N1/N2;
    f_c = torque / r_d  # Tension force
    f_cf = m_ch * r_d * r_d * w_d * w_d  # Tension due to centripetal accel.

    f_chaintension = f_c + f_cf
    # if any(t > 20000 for t in f_chaintension): # seems slower
    if f_chaintension > 20000:
        print('Warning, chain tension very high at ', f_chaintension, 'N, ', torque)

    # alpha_d = 2 * pi / n2
    # alpha_o = 2 * pi / n1
    # alpha_m = alpha_d + alpha_o

    den = mu_p * r_b * ((2 * pi / n2) + (2 * pi / n1)) / (sqrt(1 + mu_p * mu_p))
    # wl = f_chaintension * den  # loaded side
    # wu = f_cf * den  # slack side
    # work = den * (f_chaintension + f_cf)

    p_loss = n2 * (w_d / 2.0 / pi) * den * (f_chaintension + f_cf)  # N*w*SUM(W)

    e = (w_d * torque - p_loss) / (w_d * torque)

    return [e, w_d * r_d, f_chaintension]


cpdef motorbike_mech4_backup(double t,double v,double r,double rho,double cd,double jr,double area,double m,double p_tyre,double motor_torque,int n2,int n1,double gradient):
    cdef double rmg, torque, torque_air, torque_roll, torque_gradient
    # Bush diameter guessed
    torque = chain_eff_single(n2, n1, 12.7, 1.21 / 1.27, 1, v / r * n1 / n2, motor_torque, 0.11, 0.00508 * 1.2)[0] * n1 / n2 * motor_torque
    rmg = r * m * 9.81
    # Losses *add bearing and transmission losses*
    torque_air = 0.5 * r * v * v * rho * cd * area
    if v < 45.83:  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = rmg * (0.0085 + 0.018 / p_tyre + 2.06064e-5 / p_tyre * v * v)
    else:
        torque_roll = rmg * (0.018 + 3.77136e-05 * v * v) / p_tyre

    torque_gradient = rmg * gradient

    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt


cpdef double motorbike_mech_base(double t,double v,double r,double rho,double cd,double jr,double area,double m,double p_tyre,double motor_torque,int n2,int n1,double e_chain,double gradient):
    cdef double rmg, torque, torque_air, torque_roll, torque_gradient
    # Bush diameter guessed
    torque = e_chain * n1 / n2 * motor_torque
    rmg = r * m * 9.81
    # Losses *add bearing and transmission losses*
    torque_air = 0.5 * r * v * v * rho * cd * area
    if v < 45.83:  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = rmg * (0.0085 + 0.018 / p_tyre + 2.06064e-5 / p_tyre * v * v)
    else:
        torque_roll = rmg * (0.018 + 3.77136e-05 * v * v) / p_tyre

    torque_gradient = rmg * gradient

    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt


cpdef motorbike_mech4(double t,double v,double r,double rho,double cd,double jr,double area,double m,double p_tyre,double motor_torque,int n2,int n1,double gradient):
    return motorbike_mech_base(t, v, r, rho, cd, jr, area, m, p_tyre, motor_torque, n2, n1, chain_eff_single(n2, n1, 12.7, 1.21 / 1.27, 1, v / r * n1 / n2, motor_torque, 0.11, 0.00508 * 1.2)[0], gradient)

from C_interp import fast_interp
import numpy as np
cimport numpy as np

@cython.boundscheck(False) # turn off bounds checking for this func.
def motorbike_mech2(double t,double v,double r,double rho,double cd,double jr,double area,double m,double p_tyre,np.ndarray[np.double_t, ndim=1, negative_indices=False] t_mot,np.ndarray[np.double_t, ndim=1, negative_indices=False] t_mott,int n2,int n1,double e_chain,double gradient):
    return motorbike_mech_base(t, v, r, rho, cd, jr, area, m, p_tyre, fast_interp(t, t_mott, t_mot, 0, 0), n2, n1,
                               e_chain, gradient)