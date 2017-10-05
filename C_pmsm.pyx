# -*- coding: utf-8 -*-
# cython: profile=False
# bike-sim
# Copyright (C) 2017  Jonathan Blissett
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact, jonathan@blissett.me.uk
cimport cython

from libc.math cimport sqrt, atan, M_PI, cos # sin, cos, acos, exp, sqrt, fabs

cpdef C_id_pmsm(double ke,int poles, double r_s, double l_d, double l_q, double i_q, double w_m, double v_s_max):
    cdef double w_e = w_m * (poles / 2)
    cdef double a = w_e ** 2 * l_d ** 2 + r_s ** 2
    cdef double b = 2 * w_e * (w_e * l_d * ke + r_s * i_q * (l_d - l_q))
    cdef double c = (r_s * i_q + w_e * ke) ** 2 + w_e ** 2 * l_q ** 2 * i_q ** 2 - v_s_max ** 2
    return (-b + (b ** 2 - 4 * a * c) ** 0.5) / 2 / a

cpdef C_w_pmsm(double ke,int poles, double r_s, double l_d, double l_q, double i_d, double i_q, double v_s_max):
    cdef double a = l_q ** 2 * i_q ** 2 + (l_d * i_d + ke) ** 2
    cdef double b = 2 * r_s * (i_d * i_q * (l_d - l_q) + i_q * ke)
    cdef double c = r_s ** 2 * (i_d ** 2 + i_q ** 2) - v_s_max ** 2
    cdef double w_e = (-b + (b ** 2 - 4 * a * c) ** 0.5) / 2 / a
    return 2 * w_e / poles


cpdef C_vs_pmsm(double ke,int poles, double r_s, double l_d, double l_q, double i_d, double i_q, double w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    cdef double w_e = w_m * (poles / 2)
    cdef double v_d = r_s * i_d - w_e * l_q * i_q
    cdef double v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    cdef double v_s = (v_d ** 2 + v_q ** 2) ** 0.5
    return v_s

cpdef C_v_dq_pmsm(double ke,int poles, double r_s, double l_d, double l_q, double i_d, double i_q, double w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    cdef double w_e = w_m * (poles / 2)
    cdef double v_d = r_s * i_d - w_e * l_q * i_q
    cdef double v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    cdef double v_s = (v_d ** 2 + v_q ** 2) ** 0.5
    #
    cdef double v_alpha = atan(v_d / (v_q + 0.0001))
    cdef double i_alpha = atan(i_d / (i_q + 0.0001))
    cdef double power_factor = i_alpha - v_alpha
    return [v_s, v_d, v_q, power_factor]

import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t


def C_v_dq_pmsm_array(double ke,int poles, double r_s, double l_d, double l_q, np.ndarray i_d, np.ndarray i_q, np.ndarray w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    # assert i_d.dtype == DTYPE and i_q.dtype == DTYPE and w_m.dtype == DTYPE

    cdef np.ndarray w_e = w_m * (poles / 2)
    cdef np.ndarray v_d = r_s * i_d - w_e * l_q * i_q
    cdef np.ndarray v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    cdef np.ndarray v_s = (v_d ** 2 + v_q ** 2) ** 0.5
    #
    cdef np.ndarray v_alpha = np.arctan(v_d / (v_q + 0.0001))
    cdef np.ndarray i_alpha = np.arctan(i_d / (i_q + 0.0001))
    cdef np.ndarray power_factor = i_alpha - v_alpha
    return [v_s, v_d, v_q, power_factor]

cpdef C_motor_current_newton(co, double torque, double rated_torque, double error):
    cdef double i = torque / co[2]
    # print('Initial guess = ' + str(i))
    if torque > rated_torque:
        iterate = True
        count = 0
        while iterate:
            count += 1
            fx = (co[0] * i ** 3 + co[1] * i ** 2 + co[2] * i - torque)
            i -= fx / (3 * co[0] * i ** 2 + 2 * co[1] * i + co[2])
            iterate = abs(fx) > error
            # print(i, fx)
            if count > 100:
                print('NEWTON ITERATION STUCK ' + str(torque) + ' ' + str(rated_torque))
                iterate = False
    return i


cpdef Copy_motor_losses(double stator_current,double speed_rpm, double resistance=6.333e-3, double k_rpm_2=2.6853e-5,
                   double k_rpm_1=0.01528):
    cdef double resistive = stator_current ** 2 * resistance
    cdef double moving = k_rpm_2 * speed_rpm ** 2 + k_rpm_1 * speed_rpm
    cdef double total = resistive + moving
    return [total, resistive, moving]


cpdef Copy_inverter_loss(double v_bus, double v_oll, double i_o_rms, double power_factor, double l, double f_sw, double u_ce0, double u_d0, double r_c, double r_d,
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


cpdef C_torque_fw(double t_motor, double t_motor_prev, double w_m_now, double vdc, double voc, variables, co):
    # t_motor_prev = t_motor
    Ke = variables[0]
    Poles = variables[1]
    Rs = variables[2]
    Ld = variables[3]
    Lq = variables[4]
    IR = variables[5]
    L_core = variables[6]
    T_con = variables[7]
    drive_n = variables[8]
    p_batt_estimate = t_motor_prev * w_m_now
    vdc = voc - p_batt_estimate / vdc * IR
    iq_estimate = C_motor_current_newton(co, t_motor_prev, T_con, 0.01)
    if C_vs_pmsm(Ke, Poles, Rs, Ld, Lq, 0, iq_estimate, w_m_now) > (vdc * 0.9 * 0.866):
        p_weak_now = t_motor * C_w_pmsm(Ke, Poles, Rs, Ld, Lq, 0, iq_estimate, vdc * 0.9 * 0.866)
        t_motor = p_weak_now / w_m_now
        # print(str(p_batt_estimate / p_weak_now))
        iq_now = C_motor_current_newton(co, t_motor, T_con, 0.01)
        id_now = C_id_pmsm(Ke, Poles, Rs, Ld, Lq, iq_now, w_m_now, vdc * 0.9 * 0.866)
        is_now = (id_now ** 2 + iq_now ** 2) ** 0.5
        motor_loss_now = Copy_motor_losses(is_now, w_m_now * 30 / M_PI, Rs, 2.6853e-5 / 150 * L_core,
                                        0.01528 / 150 * L_core)[0]
        [vs, vd, vq, PF] = C_v_dq_pmsm(Ke, Poles, Rs, Ld, Lq, id_now, iq_now, w_m_now)
        p_drive_loss_now = drive_n * Copy_inverter_loss(vdc, vs, is_now / (2 ** 0.5) / drive_n, PF, 82e-6, 13e3, 0.8, 1,
                                                    0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
        # print(str(p_batt_now / (t_motor * w_m_now + motor_loss_now + p_drive_loss_now)))
        p_batt_now = t_motor * w_m_now + motor_loss_now + p_drive_loss_now
        vdc = voc - p_batt_now / vdc * IR
        id = C_id_pmsm(Ke, Poles, Rs, Ld, Lq, iq_now, w_m_now, vdc * 0.9 * 0.866)
    else:
        id = 0
    return [t_motor, id, vdc]