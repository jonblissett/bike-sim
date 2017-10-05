# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython


from libc.math cimport cos, M_PI # sin, cos, acos, exp, sqrt, fabs, M_PI


cpdef C_motor_losses(double stator_current,double speed_rpm, double resistance=6.333e-3, double k_rpm_2=2.6853e-5,
                   double k_rpm_1=0.01528):
    cdef double resistive = stator_current ** 2 * resistance
    cdef double moving = k_rpm_2 * speed_rpm ** 2 + k_rpm_1 * speed_rpm
    cdef double total = resistive + moving
    return [total, resistive, moving]


cpdef C_inverter_loss(double v_bus, double v_oll, double i_o_rms, double power_factor, double l, double f_sw, double u_ce0, double u_d0, double r_c, double r_d,
                      double e_ton, double e_toff, double e_d):
    # e.g. [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]
    # = inverter_losses(450,230,350,0.95,75e-6,1e4,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3)
    cdef double i_ripple = (v_bus - 2 ** 0.5 * v_oll) * v_oll / (2 * l * v_bus * f_sw)
    cdef double i_opk = 2 ** 0.5 * i_o_rms
    cdef double m = 0.57735 * v_oll * 2 ** 0.5 / v_bus # 0.57735026919 = sqrt(1/3)
    cdef double p_ct = u_ce0 * i_opk * (1 / (2 * M_PI) + m * cos(power_factor) / 8) + r_c * i_opk ** 2 * (
        1 / 8 + m * cos(power_factor) / (3 * M_PI))
    cdef double p_cd = u_d0 * i_opk * (1 / (2 * M_PI) - m * cos(power_factor) / 8) + r_d * i_opk ** 2 * (
        1 / 8 - m * cos(power_factor) / (3 * M_PI))
    cdef double i_dc = i_opk / M_PI  # DC equivalent to actual AC output current
    cdef double i_c_on = i_dc - i_ripple / 2
    cdef double i_c_off = i_dc + i_ripple / 2
    cdef double p_st = (e_ton * i_c_on + e_toff * i_c_off) * f_sw / 300.0 * v_bus / 550.0  # Test at Vce=300V, Ic=550A
    cdef double p_sd = e_d * f_sw / 300.0 * v_bus / 550.0 * i_dc
    cdef double p_total = 6.0 * (p_ct + p_cd + p_st + p_sd)
    return [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]

import numpy as np
cimport numpy as np
cpdef C_inverter_losses(np.ndarray[np.double_t, ndim=1, negative_indices=False] v_bus, np.ndarray[np.double_t, ndim=1, negative_indices=False] v_oll, np.ndarray[np.double_t, ndim=1, negative_indices=False] i_o_rms, np.ndarray[np.double_t, ndim=1, negative_indices=False] power_factor, double l, double f_sw, double u_ce0, double u_d0, double r_c, double r_d,
                      double e_ton, double e_toff, double e_d):
    # e.g. [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]
    # = inverter_losses(450,230,350,0.95,75e-6,1e4,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3)
    i_ripple = (v_bus - 2 ** 0.5 * v_oll) * v_oll / (2 * l * v_bus * f_sw)
    i_opk = 2 ** 0.5 * i_o_rms
    m = 0.57735 * v_oll * 2 ** 0.5 / v_bus # 0.57735026919 = sqrt(1/3)
    p_ct = u_ce0 * i_opk * (1 / (2 * M_PI) + m * cos(power_factor) / 8) + r_c * i_opk ** 2 * (
        1 / 8 + m * cos(power_factor) / (3 * M_PI))
    p_cd = u_d0 * i_opk * (1 / (2 * M_PI) - m * cos(power_factor) / 8) + r_d * i_opk ** 2 * (
        1 / 8 - m * cos(power_factor) / (3 * M_PI))
    i_dc = i_opk / M_PI  # DC equivalent to actual AC output current
    i_c_on = i_dc - i_ripple / 2
    i_c_off = i_dc + i_ripple / 2
    p_st = (e_ton * i_c_on + e_toff * i_c_off) * f_sw / 300.0 * v_bus / 550.0  # Test at Vce=300V, Ic=550A
    p_sd = e_d * f_sw / 300.0 * v_bus / 550.0 * i_dc
    p_total = 6.0 * (p_ct + p_cd + p_st + p_sd)
    return [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]