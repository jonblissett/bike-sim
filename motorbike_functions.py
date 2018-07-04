# -*- coding: utf-8 -*-
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
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass

import scipy.io as sio
from scipy import integrate, optimize
from scipy.integrate import ode
from C_interp import fast_interp
from C_motorbike import chain_eff_single, motorbike_mech_base, motorbike_mech4, motorbike_mech2, C_torque_limits, C_lean_calc, C_dist_calc
from C_pmsm import C_id_pmsm, C_w_pmsm, C_vs_pmsm, C_v_dq_pmsm, C_motor_current_newton, C_torque_fw
from C_losses import C_motor_losses, C_inverter_losses, C_inverter_loss
from math import cos, atan
# from numba import jit


#def motorbike_mech4(t, v, r, rho, cd, jr, area, m, p_tyre, motor_torque, n2, n1, gradient):
#    return motorbike_mech_base(t, v, r, rho, cd, jr, area, m, p_tyre, motor_torque, n2, n1,
#                               chain_eff_single(n2, n1, 12.7, 1.21 / 1.27, 1, v / r * n1 / n2, motor_torque, 0.11,
#                                                0.00508 * 1.2)[0], gradient)
#
# def motorbike_mech2(t, v, r, rho, cd, jr, area, m, p_tyre, t_mot, t_mott, n2, n1, e_chain, gradient):
#    return motorbike_mech_base(t, v, r, rho, cd, jr, area, m, p_tyre, np.interp(t, t_mott, t_mot, 0, 0), n2, n1,
#                               e_chain, gradient)

def spec_igbt(igbt):
    # ADD MODULE LEAD RESISTANCE
    print('Todo: add IGBT module lead resistance')
    if igbt == 'FS800':
        ret = {'Uce0': 0.8, 'Ud0': 1.0, 'Rc': 0.95e-3, 'Rd': 0.54e-3, 'Eon': 12e-3, 'Eoff': 25e-3,
                          'Ed': 9.5e-3,
                          'Vce_test': 300.0, 'Ic_test': 550.0, 'Fsw': 12e3}
    elif igbt == 'FF600':
        ret = {'Uce0': 0.8, 'Ud0': 0.8, 'Rc': 1.85e-3, 'Rd': 1.3e-3, 'Eon': 83e-3, 'Eoff': 72e-3,
                          'Ed': 44e-3,
                          'Vce_test': 600.0, 'Ic_test': 600.0, 'Fsw': 12e3}
    elif igbt == 'SEMiX603_SiC':
        ret = {'Uce0': 0.8, 'Ud0': 1.0, 'Rc': 1.65e-3, 'Rd': 2.6e-3, 'Eon': 17e-3, 'Eoff': 72e-3, 'Ed': 0,
                          'Vce_test': 600.0, 'Ic_test': 600.0, 'Fsw': 12e3}
    elif igbt == 'FS900':  # copied diodes from f800 as no data
        ret = {'Uce0': 0.7, 'Ud0': 1.0, 'Rc': 0.67e-3, 'Rd': 0.54e-3, 'Eon': 34e-3, 'Eoff': 36.5e-3,
                          'Ed': 17.5e-3,
                          'Vce_test': 400.0, 'Ic_test': 550.0, 'Fsw': 12e3}
    elif igbt == 'CAS325':
        ret = {'Uce0': 0.0, 'Ud0': 0.75, 'Rc': 6.5e-3, 'Rd': 4.5e-3, 'Eon': 5.6e-3, 'Eoff': 3.7e-3, 'Ed': 2.8,
                          'Vce_test': 600.0, 'Ic_test': 300.0, 'Fsw': 12e3}
    elif igbt == 'FS450':
        ret = {'Uce0': 0.8, 'Ud0': 0.8, 'Rc': 2.5e-3, 'Rd': 1.7e-3, 'Eon': 40.5e-3, 'Eoff': 56.5e-3,
                          'Ed': 39.5e-3,
                          'Vce_test': 600.0, 'Ic_test': 450.0, 'Fsw': 12e3}
    return ret


def chain_eff(n2, n1, p, m_ch, cd, w_d, torque, mu_p, r_b):
    # returns chain efficiency = e, chain speed = v.
    # N1 = driven sprocket, N2 = driving.
    # m_ch = chain mass per unit length, kg/m
    # p = pitch, mm
    # Cd = Sprocket center distance, m
    # w_d = Angular speed of N1, rad s^-1
    # torque = driving sprocket torque, Nm
    # mu_p = 0.11;            % Pin-bush friction co-efficient, from paper
    # r_b = 0.00508;          % Internal diameter of chain bush
    w_d = abs(w_d)
    torque = abs(torque)

    r_d = n2 * p / 2000 / np.pi
    r_o = r_d * n1 / n2
    # w_o = w_d*N1/N2;
    f_c = torque / r_d  # Tension force
    f_cf = m_ch * np.square(r_d) * np.square(w_d)  # Tension due to centripetal accel.

    f_chaintension = f_c + f_cf
    if any(t > 20000 for t in f_chaintension): # seems slower
    #if f_chaintension.max > 20000:
        print('Warning, chain tension very high at ', f_chaintension, 'N, ', torque)

    alpha_d = 2 * np.pi / n2
    alpha_o = 2 * np.pi / n1
    alpha_m = alpha_d + alpha_o

    den = 1 / (np.sqrt(1 + np.square(mu_p))) * mu_p * r_b * alpha_m
    wl = f_chaintension * den  # loaded side
    wu = f_cf * den  # slack side
    work = wl + wu

    p_loss = n2 * (w_d / 2 / np.pi) * work  # N*w*SUM(W)

    e = (w_d * torque - p_loss) / (w_d * torque)

    #  e(isnan(e)) = 0.9;  # for divide/0 errors - will set e to 0.9 for w=0

    speed = w_d * r_d
    return [e, speed, f_chaintension]


def chain_eff_only(n2, n1, p, m_ch, cd, w_d, torque, mu_p, r_b):
    # returns chain efficiency = e, chain speed = v.
    # N1 = driven sprocket, N2 = driving.
    # m_ch = chain mass per unit length, kg/m
    # p = pitch, mm
    # Cd = Sprocket center distance, m
    # w_d = Angular speed of N1, rad s^-1
    # torque = driving sprocket torque, Nm
    # mu_p = 0.11;            % Pin-bush friction co-efficient, from paper
    # r_b = 0.00508;          % Internal diameter of chain bush
    w_d = abs(w_d)
    torque = abs(torque)

    r_d = n2 * p / 2000 / np.pi
    r_o = r_d * n1 / n2
    # w_o = w_d*N1/N2;
    f_c = torque / r_d  # Tension force
    f_cf = m_ch * np.square(r_d) * np.square(w_d)  # Tension due to centripetal accel.

    f_chaintension = f_c + f_cf
    #  if any(t > 20000 for t in f_chaintension): # seems slower
    if f_chaintension > 20000:
        print('Warning, chain tension very high at ', f_chaintension, 'N, ', torque)

    alpha_d = 2 * np.pi / n2
    alpha_o = 2 * np.pi / n1
    alpha_m = alpha_d + alpha_o

    wl = f_chaintension / (np.sqrt(1 + np.square(mu_p))) * mu_p * r_b * alpha_m  # loaded side
    wu = f_cf / (np.sqrt(1 + np.square(mu_p))) * mu_p * r_b * alpha_m  # slack side
    work = wl + wu

    p_loss = n2 * (w_d / 2 / np.pi) * work  # N*w*SUM(W)

    e = (w_d * torque - p_loss) / (w_d * torque)

    #  e(isnan(e)) = 0.9;  # for divide/0 errors - will set e to 0.9 for w=0

    speed = w_d * r_d
    return e


def motor_torque_speed(torque_max, w_max, power_max, w_limit, n, show_figure):
    w = np.linspace(0.00001, w_max, n)
    torque = power_max / w

    #  torque(torque > torque_max) = torque_max
    #  torque = [torque_max for i in torque if i > torque_max]
    torque = [torque_max if x > torque_max else x for x in torque]
    # torque = filter(lambda x: x >= torque_max, torque)

    w = np.append(w, w_limit)
    torque = np.append(torque, 0)

    power = w * torque

    if show_figure == 1:
        fig, ax1 = plt.subplots()
        ax1.plot(w * 30.0 / np.pi, torque, 'b-')
        ax1.set_xlabel('Angular Speed (rpm)')
        # Make the y-axis label and tick labels match the line color.
        ax1.set_ylabel('Torque (Nm)', color='b')
        for tl in ax1.get_yticklabels():
            tl.set_color('b')

        ax2 = ax1.twinx()
        ax2.plot(w * 30.0 / np.pi, power / 1000, 'g-')
        ax2.set_ylabel('Power (kW)', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        plt.title('Motor Torque and Power vs Speed')
        # plt.show()
        fig.show()

    return [w, torque, power]


def braking_regen(w, torque_wheel, t_brake_lim, k):
    torque = torque_wheel
    torque[torque < 0] = -w[torque < 0] * k
    torque[torque < -t_brake_lim] = -t_brake_lim
    return torque


def lean_calc(vel, dh, dt):  # double check that the dh value sign passed to this is correct
    if dh > np.pi / 2:
        dh -= np.pi
    if dh < -np.pi / 2:
        dh += np.pi
    return atan(vel * dh / dt / 9.81)  # w = dH / dt, a_lateral = V[-1] * w


def dist_calc(d_prev, v_prev, v_prevprev, dt):
    return d_prev + (v_prev + v_prevprev) * dt / 2.0


def tyre_r_from_lean(r_nominal, w, angle):
    return r_nominal - w * (1 - cos(angle))


def gear_optimise(TT_Sim, ref_race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map, filename_ref_brake,
                  structure_map, var_name_brake, enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed):
    N0 = TT_Sim['N'][1]
    N = [N0 - 1, N0, N0 + 1]
    tmax = []

    # better algorithm is to first check if point is the best (i.e. one sim either side). if not embark on iterating in the correct direction
    # TODO need to return results of the best simulation, not the last one!
    for n in N:
        TT_Sim['N'][1] = n
        Sim = lap_analyse3(TT_Sim, ref_race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map,
                              filename_ref_brake, structure_map, var_name_brake, enable_warnings, verbosity,
                              calibration_mode, full_data_exp, battery_fixed)

        if full_data_exp:
            if Sim['ERROR']['LVC'] == 0:
                tmax.append(Sim['t'][sum(Sim['Distance_race'] < end_dist)-1])
            else:
                tmax.append(1e6)
        else:
            if Sim[14]['LVC'] == 0:
                tmax.append(np.array(Sim[6]))
            else:
                tmax.append(1e6)
        if verbosity > 1:
            print('N = ' + str(Sim['N']) + ' t = ' + str(tmax[-1]))

    if tmax[0] > tmax[1] < tmax[2]:
        # print('Minimum found')
        TT_Sim['N'][1] = N[1]
        lim = 999
        dN = 0
    elif tmax[0] < tmax[1]:
        TT_Sim['N'][1] = N[0]
        dN = -1
        lim = 1
    elif tmax[2] < tmax[1]:
        TT_Sim['N'][1] = N[2]
        dN = 1
        lim = 1
    else:
        # maximum found ... wat?
        TT_Sim['N'][1] = N[2]
        dN = 1
        lim = 1

    TT_Sim['N'][1] += dN

    while lim < 8:
        Sim = lap_analyse3(TT_Sim, ref_race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map,
                              filename_ref_brake, structure_map, var_name_brake, enable_warnings, verbosity,
                              calibration_mode, full_data_exp, battery_fixed)
        if full_data_exp:
            if Sim['ERROR']['LVC'] == 0:
                tma = Sim['t'][sum(Sim['Distance_race'] < end_dist)-1]
            else:
                tma = 1e6
        else:
            if Sim[14]['LVC'] == 0:
                tma = np.array(Sim[6])
            else:
                tma = 1e6

        if dN > 0:  # Maintaining proper list order
            tmax.append(tma)
            N.append(TT_Sim['N'][1])
            if tmax[-1] > tmax[-2]:
                dN *= -1
        else:
            tmax.insert(0, tma)
            N.insert(0, TT_Sim['N'][1])
            if tmax[0] > tmax[1]:
                dN *= -1

        if verbosity > 2:
            print(tmax)
            print(N)

        TT_Sim['N'][1] += dN

        if TT_Sim['N'][1] in N:
            break
        else:
            lim += 1

    TT_Sim = lap_analyse3(TT_Sim, ref_race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map,
                          filename_ref_brake, structure_map, var_name_brake, enable_warnings, verbosity,
                          calibration_mode, full_data_exp, battery_fixed)
    if full_data_exp:
        tma = TT_Sim['t'][sum(TT_Sim['Distance_race'] < end_dist) - 1]
    else:
        tma = np.array(TT_Sim[6])
    #if verbosity > 0:
    print('After ' + str(lim) + ' iterations, found N = ' + str(TT_Sim['N']) + ' with t = ' + str(tma))
    return TT_Sim


def lap_analyse3(TT_Sim, ref_race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map, filename_ref_brake, structure_map,
                 var_name_brake, enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed):
    ramp_start = 0.3
    ramp_time = 1.5    # Throttle ramp settings
    ramp = [ramp_start, ramp_time]
    dt = [0.01, 0.025]

    [TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = motor_torque_speed(TT_Sim['motor']['T_max'],
                                                                                            TT_Sim['motor']['W_speed_lim'],
                                                                                            TT_Sim['motor']['P_max'],
                                                                                            TT_Sim['motor']['W_lim'],
                                                                                            50, 0)
    # Import course map and reference lap data
    mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
    Course_map = mat_contents[structure_map]
    corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
    if laps > 1:
        startlap = corners[var_name_brake]  # -1 as matlab indexing starts at 1
        onelap = startlap[range(first_corner + 1, min(startlap.size - 1, 1 + last_corner))]
        startlap = startlap[range(first_corner, min(startlap.size - 1, 1 + last_corner))]
        locsmin = startlap
        for l in range(1, laps):
            locsmin = np.append(locsmin, onelap)
            last_corner += onelap.size
    else:
        locsmin = corners[var_name_brake]  # -1 as matlab indexing starts at 1
    TT_Sim = motor_saturation_interp(TT_Sim)

    locsmin = np.delete(locsmin, corner_delete) # [11, 60, 93, 107]

    # Preallocate model bike data structure
    TT_Sim['ERROR'] = {}
    TT_Sim['ERROR']['LVC'] = 0
    TT_Sim['t'] = np.array([])
    TT_Sim['v'] = np.array([])
    TT_Sim['Id'] = np.array([])
    TT_Sim['Iq'] = np.array([])
    TT_Sim['torque'] = np.array([])
    TT_Sim['Rpm'] = np.array([])
    TT_Sim['Distance'] = np.array([])
    TT_Sim['Distance_race'] = np.array([])
    if not battery_fixed:
        TT_Sim['Vdc_sim'] = np.zeros(locsmin.size + 1)  # Will store battery V at START of each corner
        TT_Sim['Vdc_sim'][first_corner] = battery_simple(TT_Sim, 0, verbosity)[0]
    TT_Sim['corner'] = {}
    TT_Sim['corner']['ID'] = []
    TT_Sim['corner']['index'] = []
    TT_Sim['Vdc'] = np.array([])
    TT_Sim['Vdc2'] = np.array([])
    TT_Sim['P'] = {}
    TT_Sim['P']['Mech'] = np.array([])
    TT_Sim['P']['Motor'] = np.array([])
    TT_Sim['P']['Drive'] = np.array([])
    TT_Sim['P']['MotorLosses'] = np.array([])
    TT_Sim['P']['DriveLosses'] = np.array([])
    TT_Sim['P']['motRLoss'] = np.array([])
    TT_Sim['P']['motwLoss'] = np.array([])
    TT_Sim['PF'] = np.array([])
    TT_Sim['Energy'] = {}
    TT_Sim['Energy']['Drive'] = np.array([])
    TT_Sim['Vd'] = np.array([])
    TT_Sim['Vq'] = np.array([])
    TT_Sim['Vs'] = np.array([])
    TT_Sim['lean'] = np.array([])
    TT_Sim['Temperature'] = {}
    TT_Sim['Temperature']['Battery'] = np.array([])

    TT_Sim['J']['linear'] = TT_Sim['constants']['m'] * TT_Sim['constants']['r'] ** 2
    j_mot_to_wheel = np.square(TT_Sim['N'][0] / TT_Sim['N'][1]) * TT_Sim['J']['motor']
    TT_Sim['J']['r'] = 2 * TT_Sim['J']['wheel'] + TT_Sim['J']['linear'] + j_mot_to_wheel  # J referred to wheel

    fw_vars = [TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'], TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'],
                    TT_Sim['motor']['Lq'], TT_Sim['battery']['IR'], TT_Sim['motor']['L_core'], TT_Sim['motor']['T_con'],
                    TT_Sim['drive']['n']]

    TT_Sim['v_flag'] = -1
    current_lap = 0
    dist_len = [0]
    if laps > 1:
        corner_range = range(first_corner, min(locsmin.size - 1, last_corner))
        for c in corner_range:
            if locsmin[c] < locsmin[c + 1]:  # if end of lap, wrap round
                TT_Sim = corner_sim_single_fw(c, locsmin, v, dt, ramp, Course_map, ref_race, TT_Sim, fw_vars,
                                              enable_warnings, calibration_mode, verbosity, battery_fixed)
            else:
                TT_Sim['Vdc_sim'][c + 1] = TT_Sim['Vdc_sim'][c]
                current_lap += 1
                dist_len.append(TT_Sim['Distance'].size)
    else:
        corner_range = range(first_corner, min(locsmin.size - 1, last_corner))
        for c in corner_range:
            TT_Sim = corner_sim_single_fw(c, locsmin, v, dt, ramp, Course_map, ref_race, TT_Sim, fw_vars,
                                          enable_warnings, calibration_mode, verbosity, battery_fixed)

    dist_len.append(TT_Sim['Distance'].size)

    if verbosity > 0:
        print('lap indexes: ' + str(dist_len))
    TT_Sim['Distance_race'] = TT_Sim['Distance'][0:(dist_len[1])]
    if locsmin.size > 1:
        d1 = ref_race.Distance[int(locsmin[1])]  # startline to first corner distance
    else:
        d1 = 0

    for i in range(1, len(dist_len) - 1):
        # print(d1, TT_Sim['Distance_race'][-1])
        TT_Sim['Distance_race'] = np.append(TT_Sim['Distance_race'], TT_Sim['Distance_race'][-1] - d1 + TT_Sim['Distance'][dist_len[i]:(dist_len[i+1])])
            #TT_Sim['Distance_race'] = np.append(TT_Sim['Distance_race'], TT_Sim['Distance'][-length:-1]+TT_Sim['Distance_race'][-1])

    if verbosity > 0:
        print('Final Energy=%.0fWh' % (TT_Sim['Energy']['Drive'] / 3600))
        print('Final cell voltage=%.2fV' % (TT_Sim['Vdc_sim'][c + 1] / TT_Sim['battery']['series']))
        print('ERROR status:', TT_Sim['ERROR'])
    i_cell = TT_Sim['Idc'] / TT_Sim['battery']['parallel']
    m_cell = TT_Sim['battery']['cellAh'] * TT_Sim['battery']['cellVnom'] / TT_Sim['battery']['E_density']
    TT_Sim['Temperature']['Battery'] = battery_heat(TT_Sim['battery']['cellIR'], m_cell, i_cell, TT_Sim['t'],
                                                    verbosity)
    TT_Sim['gradient'] = np.interp(TT_Sim['Distance'], Course_map.dist, Course_map.gradient, 0, 0)
    TT_Sim['constants']['R'] = TT_Sim['v'] * 30 / np.pi / TT_Sim['Rpm'] * TT_Sim['N'][0] / TT_Sim['N'][1]  # v/w = r

    if full_data_exp:
        return TT_Sim
    else:
        end_dist = TT_Sim['Distance_race'] < end_dist
        tmax = TT_Sim['t'][sum(end_dist)-1]
        out = ['Pmax,Tmax,Eused,Erated,Einit,tmax,N1,N2,Vend,m,dTbat,vmax,RPMmax,ERROR,Lcore,turns,Ns,Np,drives',
               TT_Sim['motor']['P_max']/1e3, max(TT_Sim['torque']), TT_Sim['Energy']['Drive'] / 3600,
               TT_Sim['battery']['E_rated'], TT_Sim['battery']['E_init'], tmax, TT_Sim['N'][0],
               TT_Sim['N'][1], TT_Sim['Vdc_sim'][c + 1], TT_Sim['constants']['m'], TT_Sim['Temperature']['Battery'][-1],
               max(TT_Sim['v']), max(TT_Sim['Rpm']), TT_Sim['ERROR'], TT_Sim['motor']['L_core'], TT_Sim['motor']['N'],
               TT_Sim['battery']['series'], TT_Sim['battery']['parallel'], TT_Sim['drive']['n']]
        if verbosity == 0:
            print(out)
        return out


def corner_sim_single(c, locsmin, v, dt, ramp, Course_map, Ref_Race, TT_Sim, enable_warnings, verbosity):
    corner_index = np.arange(locsmin[c], locsmin[c + 1] - 1)
    # corner_index_end = locsmin[c + 1]
    dt_a = dt[0]
    # dt_b = dt[1]
    ramp_start = ramp[0]
    ramp_time = ramp[1]    # Throttle ramp settings

    if TT_Sim['v_flag'] == -1:  # Set start vel to corner vel iff bike made that speed on previous
        v0 = v[corner_index[0]]
    else:
        v0 = TT_Sim['v_flag']
    d0 = Ref_Race.Distance[corner_index[0]]

    t = Ref_Race.t[corner_index[0]] + np.linspace(0, corner_index.size * dt_a - dt_a, corner_index.size)
    # print('dtA = ', str(t[0]), 'to ', str(t[-1]), 'step = ', str(t[1] - t[0]))

    wfw = np.linspace(0, TT_Sim['motor']['W_max'], 1000)
    v_s = v_dq_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'], TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'],
                    TT_Sim['motor']['Lq'], 0, TT_Sim['motor']['T_max'] / TT_Sim['constants']['Km'], wfw)[0]

    wfwlimited = wfw[v_s > 0.9 * TT_Sim['Vdc_sim'][c]]
    if wfwlimited.size == 0:
        p_weak = TT_Sim['motor']['P_max']
    else:
        p_weak = TT_Sim['motor']['T_max'] * wfwlimited[0]
    if p_weak >= TT_Sim['motor']['P_max']:
        p_weak = TT_Sim['motor']['P_max']
    else:
        if verbosity > 1:
            print('Power limited to', str(p_weak), 'by FIELD weakening')

    if TT_Sim['ERROR']['LVC'] == 1:
        if verbosity > 1:
            print('Out of battery, corner simulation ABORTED')
        return TT_Sim
    else:
        TT_Sim['ERROR']['LVC'] = 0
        TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p'] = motor_torque_speed(TT_Sim['motor']['T_max'],
                                                                                              TT_Sim['motor']['W_speed_lim'],
                                                                                              p_weak,
                                                                                              TT_Sim['motor']['W_lim'],
                                                                                              50, 0)
        # TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
        # TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
        # TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']

        V = []  # list to hold solutions
        D = []  # distances
        G = []  # course gradients
        A = []  # lean angle f(D)
        lean = []
        T_motor = []  # motor torque
        J_l = []
        J_r = []
        R = []
        H = []  # Course heading f(D)
        V.append(v0)  # Put y0 into solution
        D.append(d0)
        G.append(np.interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
        A.append(np.interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
        lean.append(0.0)
        R.append(tyre_r_from_lean(TT_Sim['constants']['r'], 0.11, A[-1]))  # BAAAAAAD - VARIABLE NEEDED FOR TYRE
        H.append(np.interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
        T_motor.append(np.interp(V[-1] / R[-1] * TT_Sim['N'][0] / TT_Sim['N'][1], TT_Sim['motor']['w'],
                                 TT_Sim['motor']['t'], 0.001, 0.001))
        J_l.append(TT_Sim['constants']['m'] * R[-1] ** 2)
        J_r.append(2 * TT_Sim['J']['wheel'] + J_l[-1] + j_mot_to_wheel)  # J referred to wheel

        solver = ode(motorbike_mech4)
        solver.set_integrator('vode', with_jacobian=False)
        solver.set_initial_value(V[0], t[0])
        solver.set_f_params(R[-1], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], J_r[-1],
                            TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'],
                            T_motor[-1], TT_Sim['N'][1], TT_Sim['N'][0], G[-1])

        # print(t[0], V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
        # corner_index2 = np.arange(locsmin[c], locsmin[c + 3] - 1)
        wheelie = False
        # heading_interp = interpolate.interp1d(Course_map.dist, Course_map.heading)

        for time in t[1:]:
            if D[-1] < Ref_Race.Distance[corner_index[-1]]:
                V.append(solver.integrate(time))
                G.append(fast_interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
                A.append(fast_interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
                R.append(tyre_r_from_lean(TT_Sim['constants']['r'], 0.11, A[-1]))  # BAAAAAAD - VARIABLE NEEDED FOR TYRE
                #  R[-1]=TT_Sim['constants']['r']
                T_motor.append(fast_interp(V[-1] / R[-1] * TT_Sim['N'][0] / TT_Sim['N'][1], TT_Sim['motor']['w'],
                                         TT_Sim['motor']['t'], T_motor[-1] / 2, T_motor[-1] / 2))
                if (time - t[0] < ramp_time) & (c != 0):  # ramp of 1.5 second
                    T_motor[-1] *= ramp_start + (1 - ramp_start) * (time - t[0]) / ramp_time

                to_max = wheelie_torque_limit(TT_Sim, R[-1], V[-1])
                #  print('to_max', to_max, T_motor[-1])
                if T_motor[-1] > to_max:
                    T_motor[-1] = to_max
                    wheelie = True

                #  T_motor[-1] = Ref_Race['constants']Km * fast_interp(time, Ref_Race.t, Ref_Race.Iq, -1, -1)

                J_l.append(TT_Sim['constants']['m'] * R[-1] ** 2)
                J_r.append(2 * TT_Sim['J']['wheel'] + J_l[-1] + j_mot_to_wheel)  # J ref to wheel
                solver.set_f_params(R[-1], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], J_r[-1],
                                    TT_Sim['constants']['area'], TT_Sim['constants']['m'],
                                    TT_Sim['constants']['p_tyre'], T_motor[-1], TT_Sim['N'][1], TT_Sim['N'][0], G[-1])
                D.append(C_dist_calc(D[-1], V[-1], V[-2], dt_a))
                # D.append(D[-1] + np.squeeze(V[-1] + V[-2]) * dt_a / 2.0)

                H.append(fast_interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
                lean.append(C_lean_calc(V[-1], H[-1] - H[-2], dt_a))


                # print(D[-1], H[-1], V[-1], lean[-1]*180/np.pi, A[-1]*180/np.pi)
                # print(time, V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
                # todo put field weakening calculation in here? include IR drop at least.
                # Could just calc. motor curve for two adjacent points
                if not solver.successful():
                    print('Warning: integration not successful')
            else:
                # print('Simulating too far')
                # better way would be to look at 'time' and fill these arrays in one go, then quit loop
                D.append(D[-1])
                H.append(H[-1])
                V.append(V[-1])
                G.append(G[-1])
                R.append(R[-1])
                lean.append(lean[-1])
                T_motor.append(T_motor[-1])
        if wheelie & (verbosity > 1):
            print('Wheelie alert!')
        # V = np.squeeze(V)
        # D = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(V, t, initial=0))
        V = np.squeeze(V)
        R = np.squeeze(R)
        T_motor = np.squeeze(T_motor)
        # A = np.squeeze(A)
        lean = np.squeeze(lean)
        D = np.squeeze(D)
        # gradient_save = gradient
        # gradient = np.interp(D, TT_map.dist, TT_map.gradient, 0.0, 0.0)

        MOTORSPEED = V / (1.0 / 60 * TT_Sim['N'][1] / TT_Sim['N'][0] * 2 * np.pi * R)  # in rpm
        MOTORTORQUE = T_motor

        T = t
        # t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[0]],
        #                 corner_index.size)  # as previous but flipped

        # print('dtB = ', str(t[1] - t[0]))
        t = Ref_Race.t[corner_index[-1]] + np.linspace(0, 2 * dt_b * (1 - corner_index.size), 2 * corner_index.size)
        dt_b = t[0] - t[1]
        # print('dtB = ', str(t[1] - t[0]))

        TBrake_t = t
        TBrake = -TT_Sim['brake']['PeakTorque'] * TT_Sim['N'][1] / TT_Sim['N'][0] * np.ones(t.shape)
        rampIndex = TBrake_t > (TBrake_t[0] - TT_Sim['brake']['RampTime'])
        TBrake[rampIndex] = np.linspace(0, TBrake[-1], sum(rampIndex))
        TBrake_t = np.flipud(TBrake_t)
        TBrake = np.flipud(TBrake)

        # plt.close()
        # plt.plot(Tbraket, Tbrake)
        # plt.show()

        v0 = v[corner_index[-1]]
        d0 = Ref_Race.Distance[corner_index[-1]]

        #  gradientt = t
        #  gradient = np.interp(Ref_Race.Distance[corner_index], Course_map.dist, Course_map.gradient, 0,
        #                       0)  # Initially gradient(t) = same as lap data
        #  gradientt = np.flipud(gradientt)

        # print('a=', a)
        e_chain = 1
        # V2 = odeint(motorbike_mech2, v0, t,
        #            args=(TT_Sim['constants']['r'], TT_Sim['constants']['r']ho, TT_Sim['constants']['cd'], TT_Sim['J']['r'],
        #                  TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'], TBrake, TBrake_t,
        #                  TT_Sim['N'][1], TT_Sim['N'][0], e_chain, gradient, gradientt))

        V2 = []  # list to hold solutions
        D2 = []
        H2 = []
        G2 = []  # course gradients
        lean2 = []
        V2.append(v0)  # Put y0 into solution
        D2.append(d0)
        G2.append(fast_interp(D2[-1], Course_map.dist, Course_map.gradient, 0, 0))
        H2.append(fast_interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
        lean2.append(0.0)
        # print(str(a))
        solver = ode(motorbike_mech2)
        # solver.set_integrator('dopri5')     #   5th order runge-kutta
        solver.set_integrator('vode', with_jacobian=False)
        solver.set_initial_value(v0, t[0])
        solver.set_f_params(TT_Sim['constants']['r'], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], TT_Sim['J']['r'],
                            TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'], TBrake, TBrake_t,
                            TT_Sim['N'][1], TT_Sim['N'][0], e_chain, G[-1])

        v_max = max(V) * 1.1  # Need to stop simulating a bit after intersection, for gradient delta V
        # print('v_max = ', v_max, 't0 = ', t[0])
        for time in t[1:]:
            if V2[-1] >= v_max:
                # print('Simulating too far')
                # better way would be to look at 'time' and fill these arrays in one go, then quit loop
                V2.append(v_max)  # would be neater to break the loop here, but messes up the array lengths
                D2.append(D2[-1] + np.squeeze(V2[-1] + V2[-2]) * -dt_b / 2.0)
                H2.append(H2[-1])
                lean2.append(lean2[-1])
            else:
                V2.append(solver.integrate(time))
                if not solver.successful():
                    print('Warning: integration not successful')
                D2.append(D2[-1] + np.squeeze(V2[-1] + V2[-2]) * -dt_b / 2.0)
                G2.append(fast_interp(D2[-1], Course_map.dist, Course_map.gradient, 0, 0))
                H2.append(fast_interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
                lean2.append(C_lean_calc(V2[-1], H2[-2] - H2[-1], -dt_b))
                # print('V2 = ', V2[j], 't = ', t[j])
                solver.set_f_params(TT_Sim['constants']['r'], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], TT_Sim['J']['r'],
                                    TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'], TBrake,
                                    TBrake_t, TT_Sim['N'][1], TT_Sim['N'][0], e_chain, G[-1])

        V2 = np.squeeze(V2)
        D2 = np.squeeze(D2)
        # D2 = np.squeeze(Ref_Race.Distance[corner_index[-1]] + integrate.cumtrapz(V2, t, initial=0))
        lean2 = np.squeeze(lean2)
        # T2 = np.flipud(t)
        D2 = np.flipud(D2)
        V2 = np.flipud(V2)
        MOTORSPEED2 = V2 / (1.0 / 60 * TT_Sim['N'][1] / TT_Sim['N'][0] * 2 * np.pi * TT_Sim['constants']['r'])  # in rpm
        lean2 = np.flipud(lean2)
        # gradientt = np.flipud(gradientt)
        # gradient = np.flipud(graidient)

        #  fig3 = plt.figure(3)
        # ax = fig3.add_subplot(1, 1, 1)
        # ax.plot(D, V, TT_Race.Distance, v, TT_map.dist, TT_map.gradient * 100, D2, V2)
        #  plt.xlabel('Distance')
        #  plt.ylabel('v')
        #  plt.xlim(D[0], D[-1])
        # plt.plot(T2, V2, time[corner_index], v[corner_index], gradientt, gradient * 100)
        # plt.xlim(gradientt[0], gradientt[-1])
        # plt.ylim(-10, 100)
        # fig3.show()
        # plt.show()

        # interp values of V2 on D

        if np.all(np.diff(D2) > 0):
            V2i = np.interp(D, D2, V2, -10, -20)
        else:
            D2, indices = np.unique(D2, return_index=True)
            V2 = V2[indices]
            print('ERROR, D2 not increasing - i.e. bike stopped or reversing, duplicate points deleted')
            V2i = np.interp(D, D2, V2, -10, -20)

        dout = np.argwhere(np.isclose(V2i, V, atol=0.1))  # match within 0.5m/s
        if dout.size == 0:  # try again with bigger tolerance
            dout = np.argwhere(np.isclose(V2i, V, atol=0.5))
            if verbosity > 1:
                print('No. of intersections = %d' % dout.size)

        if dout.size == 0:
            if enable_warnings:
                plt.close()
                fig4 = plt.figure(4)
                ax = fig4.add_subplot(1, 1, 1)
                # ax.plot(Ref_Race.t, v, T, V, '-o', gradientt, gradient * 100)#, Ref_Race.t, 100*np.interp(Ref_Race.Distance,TT_map.dist,TT_map.gradient))#, np.flipud(t), V2, '-o'
                ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o', Course_map.dist, Course_map.gradient * 100)
                plt.xlim(D[0], D[-1])
                plt.ylim(0, V.max() * 1.2)
                ax.plot(D[dout], V[dout], 'ro')
                fig4.show()
                plt.show()
            dout = Ref_Race.Distance[corner_index[-1]]  # BIT BAD THIS - causes jump in Vel
            TT_Sim['v_flag'] = V[(D < dout)]
            TT_Sim['v_flag'] = TT_Sim['v_flag'][-1]
            if verbosity > 1:
                print('############################################################')
                print('Bike too slow on corner %d' % c, ', perhaps you used an Agni motor?')
                print('Higher motor torque required to achieve desired corner speed')
                print('End vel = %.1f' % TT_Sim['v_flag'])
                print('############################################################')
        else:
            dout = np.median([D[dout]])
            TT_Sim['v_flag'] = -1

        # print('dout =', str(dout))
        # print('Braking point = %d' % dout, 'm')
        D2i = np.squeeze((D2 > dout) & (D2 < 60725))
        Di = np.squeeze(D <= dout)

        #  Vboth = np.hstack((V[:dout], V2i[dout:]))

        Vboth = np.hstack((V[Di], V2[D2i]))
        Dboth = np.hstack((D[Di], D2[D2i]))
        ToBoth = np.hstack((MOTORTORQUE[Di], TBrake[D2i]))
        RPMBoth = np.hstack((MOTORSPEED[Di], MOTORSPEED2[D2i]))
        leanBoth = np.hstack((lean[Di], lean2[D2i]))

        dt = dt_b  # [1] - T[0]
        # Tboth = T[0] + np.arange(0, Vboth.size * dt - dt, dt)

        Tboth = np.hstack((T[Di], T[Di][-1] + dt * np.arange(1, (D2[D2i].size + 1))))
        # print(Tboth.size, ' =? ', Vboth.size, Ref_Race.Distance[corner_index[-1]])
        # Dboth = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(Vboth, Tboth, initial=0))

        # fig5 = plt.figure(5)
        # ax = fig5.add_subplot(1, 1, 1)
        # ax.plot(TT_Race.Distance, v, Dboth, Vboth)
        # plt.xlim(Dboth[0], Dboth[-1])
        # fig5.show()

        # PLOT THIS ONE
        if enable_warnings:
            fig6 = plt.figure(6)
            ax = fig6.add_subplot(1, 1, 1)
            ax.plot(Ref_Race.t, v, Tboth, Vboth, TT_Sim['t'], TT_Sim['v'])  # , Tboth, ToBoth
            plt.xlim(Tboth[0], Tboth[-1])
            fig6.show()
            plt.show()

        if verbosity > 1:
            tgained = Ref_Race.t[corner_index[-1]] - Tboth[-1]
            print('Time gained = %.1f s' % tgained, ' on corner %d' % c)

        if TT_Sim['t'].size == 0:  # if not [] == true
            TT_Sim['t'] = np.hstack((TT_Sim['t'], Tboth))
        else:
            TT_Sim['t'] = np.hstack((TT_Sim['t'], Tboth - Tboth[0] + TT_Sim['t'][-1]))
        TT_Sim['v'] = np.hstack((TT_Sim['v'], Vboth))
        TT_Sim['Rpm'] = np.hstack((TT_Sim['Rpm'], RPMBoth))
        TT_Sim['torque'] = np.hstack((TT_Sim['torque'], ToBoth))
        TT_Sim['Distance'] = np.hstack((TT_Sim['Distance'], Dboth))
        TT_Sim['lean'] = np.hstack((TT_Sim['lean'], leanBoth))

        # Now postprocess to find energy
        TT_Sim['t'], indices = np.unique(TT_Sim['t'], return_index=True)  # Remove stops in time
        TT_Sim['v'] = TT_Sim['v'][indices]
        TT_Sim['Rpm'] = TT_Sim['Rpm'][indices]
        TT_Sim['torque'] = TT_Sim['torque'][indices]
        TT_Sim['Distance'] = TT_Sim['Distance'][indices]
        TT_Sim['lean'] = TT_Sim['lean'][indices]

        # Limit braking torque to rear wheel regenerative
        torque = braking_regen(TT_Sim['Rpm'] / 30 * np.pi * TT_Sim['N'][1] / TT_Sim['N'][0],
                               TT_Sim['torque'] * TT_Sim['N'][0] / TT_Sim['N'][1],
                               TT_Sim['brake']['LimitTorque'], TT_Sim['brake']['k_wt'])

        torque_motor = torque / TT_Sim['N'][0] * TT_Sim['N'][1]
        # TT_Sim['Iq'] = torque_motor / TT_Sim['constants']['Km']
        TT_Sim['Iq'] = motor_currents(TT_Sim, torque_motor)
        [total_loss, resistive, moving] = motor_losses(TT_Sim['Is'], TT_Sim['Rpm'], TT_Sim['motor']['Rs'],
                                                       TT_Sim['motor']['k_rpm'][2], TT_Sim['motor']['k_rpm'][1])
        TT_Sim['P']['MotorLosses'] = total_loss
        TT_Sim['P']['motRLoss'] = resistive
        TT_Sim['P']['motwLoss'] = moving
        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'],
                                                    TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'], TT_Sim['motor']['Lq'],
                                                    0, TT_Sim['Iq'], TT_Sim['Rpm'] / 30 * np.pi)
        TT_Sim['Vs'] = v_s
        TT_Sim['Vd'] = v_d
        TT_Sim['Vq'] = v_q
        TT_Sim['PF'] = power_factor
        TT_Sim['Vdc'] = np.hstack((TT_Sim['Vdc'], TT_Sim['Vdc_sim'][c] * np.ones(TT_Sim['t'].size - TT_Sim['Vdc'].size)))

        TT_Sim['P']['DriveLosses'] = TT_Sim['drive']['n'] * inverter_losses(TT_Sim['Vdc2'], TT_Sim['Vs'],
                                                                            TT_Sim['Is'] / np.sqrt(2) /
                                                                            TT_Sim['drive']['n'], TT_Sim['PF'],
                                                                            TT_Sim['motor']['Lq']*3/2, TT_Sim)[0]
        TT_Sim['P']['Mech'] = TT_Sim['Rpm'] / 30 * np.pi * torque_motor
        TT_Sim['P']['Motor'] = TT_Sim['P']['Mech'] + TT_Sim['P']['MotorLosses']
        TT_Sim['P']['Drive'] = TT_Sim['P']['Motor'] + TT_Sim['P']['DriveLosses']
        # TT_Sim['Idc'] = TT_Sim['P']['Drive'] / TT_Sim['Vdc']
        TT_Sim['Energy']['Drive'] = np.trapz(TT_Sim['P']['Drive'], TT_Sim['t'])
        TT_Sim['Vdc_sim'][c + 1], battery_status = battery_simple(TT_Sim, TT_Sim['Energy']['Drive'] / 3600, verbosity)
        if battery_status == 'LVC':
            TT_Sim['ERROR']['LVC'] = 1
        if verbosity > 1:
            print('Energy=%dWh' % (TT_Sim['Energy']['Drive'] / 3600))
            print('Voltage=%.1fV' % (TT_Sim['Vdc_sim'][c + 1]))
        return TT_Sim


def corner_sim_single_fw(c, locsmin, v, dt, ramp, Course_map, Ref_Race, TT_Sim, fw_vars, enable_warnings, calibration_mode, verbosity, battery_fixed):
    if TT_Sim['ERROR']['LVC'] == 1:
        if verbosity > 1:
            print('Out of battery, corner simulation ABORTED')
        return TT_Sim
    else:
        j_mot_to_wheel = np.square(TT_Sim['N'][0] / TT_Sim['N'][1]) * TT_Sim['J']['motor']
        TT_Sim['ERROR']['LVC'] = 0
        corner_index = np.arange(locsmin[c], locsmin[c + 1] - 1)
        # corner_index_end = locsmin[c + 1]
        #np.append(TT_Sim['corner']['ID'], np.array([c]))
        #np.append(TT_Sim['corner']['index'], np.array([len(TT_Sim['v'])]))
        TT_Sim['corner']['ID'].append(c)
        TT_Sim['corner']['index'].append(len(TT_Sim['v']))

        dt_a = dt[0]
        dt_b = dt[1]
        ramp_start = ramp[0]
        ramp_time = ramp[1]  # Throttle ramp settings

        # TT_Sim['v_flag'] = -1
        if TT_Sim['v_flag'] == -1:  # Set start vel to corner vel iff bike made that speed on previous
            v0 = v[int(corner_index[0])]
        else:
            v0 = TT_Sim['v_flag']
        d0 = Ref_Race.Distance[int(corner_index[0])]

        t = Ref_Race.t[int(corner_index[0])] + np.linspace(0, 13*corner_index.size * dt_a - dt_a, 13*corner_index.size)
        # print('dtA = ', str(t[0]), 'to ', str(t[-1]), 'step = ', str(t[1] - t[0]))

        if w_fw(TT_Sim['motor']['W_max'], TT_Sim, TT_Sim['Vdc_sim'][c], motor_currents(TT_Sim, TT_Sim['motor']['T_max'])) > 0:
            p_weak = TT_Sim['motor']['T_max'] * C_w_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'],
                                                         TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'],
                                                         TT_Sim['motor']['Lq'], 0,
                                                         motor_currents(TT_Sim, TT_Sim['motor']['T_max']),
                                                         TT_Sim['Vdc_sim'][c])
            if p_weak >= TT_Sim['motor']['P_max']:
                p_weak = TT_Sim['motor']['P_max']
            else:
                if verbosity > 1:
                    print('Power limited to', str(p_weak), 'by FIELD weakening')
            TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p'] = \
                motor_torque_speed(TT_Sim['motor']['T_max'], TT_Sim['motor']['W_speed_lim'], p_weak, TT_Sim['motor']['W_lim'],
                                   50, 0)

        # TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
        # TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
        # TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']

        V = []  # list to hold solutions
        D = []  # distances
        G = []  # course gradients
        A = []  # lean angle f(D)
        lean = []
        T_motor = []  # motor torque
        J_l = []
        J_r = []
        R = []
        H = []  # Course heading f(D)
        Vdc = []
        Id = []
        Id.append(0)
        Vdc.append(TT_Sim['Vdc_sim'][c])
        V.append(v0)  # Put y0 into solution
        # print('ZZZZZ ' + str(Vdc[-1]) + ' ' + str(V[-1]))
        D.append(d0)
        G.append(fast_interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
        A.append(fast_interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
        lean.append(0.0)
        R.append(tyre_r_from_lean(TT_Sim['constants']['r'], 0.11, A[-1]))  # BAAAAAAD - VARIABLE NEEDED FOR TYRE
        H.append(fast_interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))

        T_motor.append(fast_interp(V[-1] / R[-1] * TT_Sim['N'][0] / TT_Sim['N'][1], TT_Sim['motor']['w'],
                                   TT_Sim['motor']['t'], 0.001, 0.001))
        to_max = C_torque_limits(R[-1], V[-1], TT_Sim['N'][1], TT_Sim['N'][0], TT_Sim['constants']['m'],
                                 TT_Sim['constants']['b'], TT_Sim['constants']['h'], TT_Sim['constants']['rho'],
                                 TT_Sim['constants']['cd'], TT_Sim['constants']['area'],
                                 TT_Sim['constants']['mu_tyre'])[0]
        if T_motor[-1] > to_max:
            T_motor[-1] = to_max

        J_l.append(TT_Sim['constants']['m'] * R[-1] ** 2)
        J_r.append(2 * TT_Sim['J']['wheel'] + J_l[-1] + j_mot_to_wheel)  # J referred to wheel

        solver = ode(motorbike_mech4)
        solver.set_integrator('vode', with_jacobian=False)
        solver.set_initial_value(V[0], t[0])
        solver.set_f_params(R[-1], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], J_r[-1],
                            TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'],
                            T_motor[-1], TT_Sim['N'][1], TT_Sim['N'][0], G[-1])

        # print(t[0], V[-1], D[-1], G[-1], ,T_motor[-1], J_r[-1], A[-1], R[-1])
        # corner_index2 = np.arange(locsmin[c], locsmin[c + 3] - 1)
        wheelie = False
        # heading_interp = interpolate.interp1d(Course_map.dist, Course_map.heading)
        timr = 0
        for times in t[1:]:
            if D[-1] < Ref_Race.Distance[int(corner_index[-1])]:
                V.append(solver.integrate(times))
                G.append(fast_interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
                A.append(fast_interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
                if c == 0:  # No lean on start line
                    R.append(TT_Sim['constants']['r'])
                else:
                    R.append(tyre_r_from_lean(TT_Sim['constants']['r'], 0.11, A[-1]))
                w_m_now = V[-1] / R[-1] * TT_Sim['N'][0] / TT_Sim['N'][1]
                T_motor.append(fast_interp(w_m_now, TT_Sim['motor']['w'], TT_Sim['motor']['t'],
                                           T_motor[-1] / 2, T_motor[-1] / 2))

                # Limit motor torque for wheelies
                to_max = C_torque_limits(R[-1], V[-1], TT_Sim['N'][1], TT_Sim['N'][0], TT_Sim['constants']['m'],
                                         TT_Sim['constants']['b'], TT_Sim['constants']['h'], TT_Sim['constants']['rho'],
                                         TT_Sim['constants']['cd'], TT_Sim['constants']['area'],
                                         TT_Sim['constants']['mu_tyre'])[0]
                # print(T_motor[-1], to_max)
                if T_motor[-1] > to_max:
                    T_motor[-1] = to_max
                    wheelie = True

                # Ramp motor torque by rider, unless startline or missed corner
                #timrs = time.time()
                # this if takes about 0.5s exec time over TT lap! code inside neglegible
                if (times - t[0] < ramp_time) & (c != 0) & (TT_Sim['v_flag'] == -1):
                    T_motor[-1] *= ramp_start + (1 - ramp_start) * (times - t[0]) / ramp_time
                #timr += time.time()-timrs
                # Limit motor torque for field weakening
                # T_motor[-1], Id_new, Vdc_new = torque_fw(TT_Sim, T_motor[-1], T_motor[-2], w_m_now, Vdc[-1], TT_Sim['Vdc_sim'][c])
                # print('T pre fw:' + str(T_motor[-1]) + 'Rpm' + str(w_m_now*30/np.pi))
                T_motor[-1], id_new, vdc_new = C_torque_fw(T_motor[-1], T_motor[-2], w_m_now, Vdc[-1],
                                                           TT_Sim['Vdc_sim'][c], fw_vars, TT_Sim['motor']['co'])
                Id.append(id_new)
                Vdc.append(vdc_new)
                # print('T postfw:' + str(T_motor[-1]))

                if calibration_mode:
                    T_motor[-1] = motor_torque(TT_Sim['motor']['co'], fast_interp(times, Ref_Race.t, Ref_Race.Iq, 0, 0))
                    Id[-1] = 0
                    # print(T_motor[-1])
                # print('Torque = ' + str(T_motor[-1]))

                J_l.append(TT_Sim['constants']['m'] * R[-1] ** 2)
                J_r.append(2 * TT_Sim['J']['wheel'] + J_l[-1] + j_mot_to_wheel)  # J ref to wheel
                solver.set_f_params(R[-1], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'], J_r[-1],
                                    TT_Sim['constants']['area'], TT_Sim['constants']['m'],
                                    TT_Sim['constants']['p_tyre'], T_motor[-1], TT_Sim['N'][1], TT_Sim['N'][0], G[-1])
                D.append(C_dist_calc(D[-1], V[-1], V[-2], dt_a))

                H.append(fast_interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
                lean.append(C_lean_calc(V[-1], H[-1] - H[-2], dt_a))

                # print(D[-1], H[-1], V[-1], lean[-1]*180/np.pi, A[-1]*180/np.pi)
                # print(times, V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
                if not solver.successful():
                    print('Warning: integration not successful')
            else:
                # print('Simulating too far')
                # better way would be to look at 'time' and fill these arrays in one go, then quit loop
                D.append(D[-1])
                H.append(H[-1])
                V.append(V[-1])
                G.append(G[-1])
                R.append(R[-1])
                lean.append(lean[-1])
                T_motor.append(T_motor[-1])
                Vdc.append(TT_Sim['Vdc_sim'][c])
                Id.append(0)
        if wheelie & (verbosity > 1):
            print('Wheelie alert!')
        # print('duration = %.4f s' % timr)#(time.time() - timr))


        # V = np.squeeze(V)
        # D = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(V, t, initial=0))
        V = np.squeeze(V)
        Vdc = np.squeeze(Vdc)
        Id = np.squeeze(Id)
        R = np.squeeze(R)
        T_motor = np.squeeze(T_motor)
        # A = np.squeeze(A)
        lean = np.squeeze(lean)
        D = np.squeeze(D)
        # gradient_save = gradient
        # gradient = np.interp(D, TT_map.dist, TT_map.gradient, 0.0, 0.0)

        MOTORSPEED = V / (1.0 / 60 * TT_Sim['N'][1] / TT_Sim['N'][0] * 2 * np.pi * R)  # in rpm
        MOTORTORQUE = T_motor

        T = t
        # t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[0]],
        #                 corner_index.size)  # as previous but flipped

        # print('dtB = ', str(t[1] - t[0]))
        t = Ref_Race.t[int(corner_index[-1])] + np.linspace(0, 3.5 * dt_b * (1 - corner_index.size), 5 * corner_index.size)
        dt_b = t[0] - t[1]
        # print('dtB = ', dt_b)

        TBrake_t = t
        TBrake = -TT_Sim['brake']['PeakTorque'] * TT_Sim['N'][1] / TT_Sim['N'][0] * np.ones(t.shape)
        rampIndex = TBrake_t > (TBrake_t[0] - TT_Sim['brake']['RampTime'])
        TBrake[rampIndex] = np.linspace(0, TBrake[-1], sum(rampIndex))
        TBrake_t = np.flipud(TBrake_t)
        TBrake = np.flipud(TBrake)

        # plt.close()
        # plt.plot(TBrake_t, TBrake)
        # plt.show()

        v0 = v[int(corner_index[-1])]
        d0 = Ref_Race.Distance[int(corner_index[-1])]

        #  gradientt = t
        #  gradient = fast_interp(Ref_Race.Distance[corner_index], Course_map.dist, Course_map.gradient, 0,
        #                       0)  # Initially gradient(t) = same as lap data
        #  gradientt = np.flipud(gradientt)

        # print('a=', a)
        e_chain = 1
        # V2 = odeint(motorbike_mech2, v0, t,
        #            args=(TT_Sim['constants']['r'], TT_Sim['constants']['r']ho, TT_Sim['constants']['cd'], TT_Sim['J']['r'],
        #                  TT_Sim['constants']['area'], TT_Sim['constants']['m'], TT_Sim['constants']['p_tyre'], TBrake, TBrake_t,
        #                  TT_Sim['N'][1], TT_Sim['N'][0], e_chain, gradient, gradientt))
        V2 = []  # list to hold solutions
        D2 = []
        H2 = []
        G2 = []  # course gradients
        lean2 = []
        V2.append(v0)  # Put y0 into solution
        D2.append(d0)
        G2.append(fast_interp(D2[-1], Course_map.dist, Course_map.gradient, 0, 0))
        H2.append(fast_interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
        lean2.append(0.0)
        # print(str(a))
        solver = ode(motorbike_mech2)
        solver.set_integrator('vode', with_jacobian=False)
        solver.set_initial_value(v0, t[0])
        # TBrake = -TT_Sim['brake']['PeakTorque'] * TT_Sim['N'][1] / TT_Sim['N'][0] * np.ones(t.shape)
        solver.set_f_params(TT_Sim['constants']['r'], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'],
                            TT_Sim['J']['r'], TT_Sim['constants']['area'], TT_Sim['constants']['m'],
                            TT_Sim['constants']['p_tyre'], TBrake, TBrake_t, TT_Sim['N'][1], TT_Sim['N'][0], e_chain,
                            G2[-1])

        v_max = max(V) * 1.1  # Need to stop simulating a bit after intersection, for gradient delta V
        # print('v_max = ', v_max, 't0 = ', t[0])

        for times in t[1:]:
            if V2[-1] >= v_max:
                # print('Simulating too far')
                # better way would be to look at 'time' and fill these arrays in one go, then quit loop
                V2.append(v_max)  # would be neater to break the loop here, but messes up the array lengths
                D2.append(C_dist_calc(D2[-1], V2[-1], V2[-2], -dt_b))
                H2.append(H2[-1])
                lean2.append(lean2[-1])
            else:
                V2.append(solver.integrate(times))
                if not solver.successful():
                    print('Warning: integration not successful')
                D2.append(C_dist_calc(D2[-1], V2[-1], V2[-2], -dt_b))
                G2.append(fast_interp(D2[-1], Course_map.dist, Course_map.gradient, 0, 0))
                H2.append(fast_interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
                lean2.append(C_lean_calc(V2[-1], H2[-2] - H2[-1], dt_b))
                # print('V2 = ', V2[j], 't = ', t[j])
                # TWO THINGS TO LIMIT by: stoppie and friction
                # [torque_MaxRear, torque_MaxFront] = wheel_forces(TT_Sim['Distance'], 1.415, TT_Sim['constants']['h'],
                #                                                 TT_Sim['constants']['b'], TT_Sim['constants']['r'],
                #                                                 TT_Sim['constants']['m'], V2[-1] ** 2 *
                #                                                 TT_Sim['constants']['rho'] * TT_Sim['constants']['cd'] *
                #                                                 TT_Sim['constants']['area'] / 2.0,
                #                                                 0 * TT_Sim['N'][0] / TT_Sim['N'][1], 0)[2:4]
                # print(str(torque_MaxRear) + str(torque_MaxFront))
                solver.set_f_params(TT_Sim['constants']['r'], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'],
                                    TT_Sim['J']['r'], TT_Sim['constants']['area'], TT_Sim['constants']['m'],
                                    TT_Sim['constants']['p_tyre'], TBrake, TBrake_t, TT_Sim['N'][1], TT_Sim['N'][0],
                                    e_chain, G2[-1])
                # print('t =', times, 'v = ', np.squeeze(V2[-1]), 'T =', fast_interp(times, TBrake_t, TBrake, 0, 0), 'Dist', D2[-1])
                #print(TT_Sim['constants']['r'], TT_Sim['constants']['rho'], TT_Sim['constants']['cd'],
                #                    TT_Sim['J']['r'], TT_Sim['constants']['area'], TT_Sim['constants']['m'],
                #                    TT_Sim['constants']['p_tyre'], TT_Sim['N'][1], TT_Sim['N'][0],
                #                    e_chain, G2[-1])


        V2 = np.squeeze(V2)
        D2 = np.squeeze(D2)
        # D2 = np.squeeze(Ref_Race.Distance[corner_index[-1]] + integrate.cumtrapz(V2, t, initial=0))
        lean2 = np.squeeze(lean2)
        # T2 = np.flipud(t)
        D2 = np.flipud(D2)
        V2 = np.flipud(V2)
        MOTORSPEED2 = V2 / (1.0 / 30 * TT_Sim['N'][1] / TT_Sim['N'][0] * np.pi * TT_Sim['constants']['r'])  # in rpm
        lean2 = np.flipud(lean2)
        # gradientt = np.flipud(gradientt)
        # gradient = np.flipud(graidient)

        #  fig3 = plt.figure(3)
        # ax = fig3.add_subplot(1, 1, 1)
        # ax.plot(D, V, TT_Race.Distance, v, TT_map.dist, TT_map.gradient * 100, D2, V2)
        #  plt.xlabel('Distance')
        #  plt.ylabel('v')
        #  plt.xlim(D[0], D[-1])
        # plt.plot(T2, V2, times[corner_index], v[corner_index], gradientt, gradient * 100)
        # plt.xlim(gradientt[0], gradientt[-1])
        # plt.ylim(-10, 100)
        # fig3.show()
        # plt.show()

        # interp values of V2 on D

        if np.all(np.diff(D2) > 0):
            V2i = np.interp(D, D2, V2, -10, -20)
        else:
            D2, indices = np.unique(D2, return_index=True)
            V2 = V2[indices]
            print('ERROR, D2 not increasing - i.e. bike stopped or reversing, duplicate points deleted')
            V2i = np.interp(D, D2, V2, -10, -20)

        dout = np.argwhere(np.isclose(V2i, V, atol=0.1))  # match within 0.5m/s
        if dout.size == 0:  # try again with bigger tolerance
            dout = np.argwhere(np.isclose(V2i, V, atol=0.5))
            if verbosity > 1:
                print('No. of intersections = %d' % dout.size)

        if dout.size == 0:
            if enable_warnings:
                plt.close()
                fig4 = plt.figure(4)
                ax = fig4.add_subplot(1, 1, 1)
                # ax.plot(Ref_Race.t, v, T, V, '-o', gradientt, gradient * 100)#, Ref_Race.t, 100*np.interp(Ref_Race.Distance,TT_map.dist,TT_map.gradient))#, np.flipud(t), V2, '-o'
                ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o', Course_map.dist, Course_map.gradient * 100)
                plt.xlim(D[0], D[-1])
                plt.ylim(0, V.max() * 1.2)
                ax.plot(D[dout], V[dout], 'ro')
                fig4.show()
                plt.show()
            dout = Ref_Race.Distance[int(corner_index[-1])]  # BIT BAD THIS - causes jump in Vel
            TT_Sim['v_flag'] = V[(D < dout)]
            TT_Sim['v_flag'] = TT_Sim['v_flag'][-1]
            if verbosity > 1:
                print('############################################################')
                print('Bike too slow on corner %d' % c, ', perhaps you used an Agni motor?')
                print('Higher motor torque required to achieve desired corner speed')
                print('End vel = %.1f' % TT_Sim['v_flag'])
                print('############################################################')
        else:
            dout = np.median([D[dout]])
            TT_Sim['v_flag'] = -1

        # print('dout =', str(dout))
        # print('Braking point = %d' % dout, 'm')
        D2i = np.squeeze((D2 > dout) & (D2 < 60725))
        Di = np.squeeze(D <= dout)

        #  Vboth = np.hstack((V[:dout], V2i[dout:]))

        Vboth = np.hstack((V[Di], V2[D2i]))
        Dboth = np.hstack((D[Di], D2[D2i]))
        ToBoth = np.hstack((MOTORTORQUE[Di], TBrake[D2i]))
        RPMBoth = np.hstack((MOTORSPEED[Di], MOTORSPEED2[D2i]))
        leanBoth = np.hstack((lean[Di], lean2[D2i]))

        dt = dt_b  # [1] - T[0]
        # Tboth = T[0] + np.arange(0, Vboth.size * dt - dt, dt)

        Tboth = np.hstack((T[Di], T[Di][-1] + dt * np.arange(1, (D2[D2i].size + 1))))
        Idboth = np.hstack((Id[Di], np.zeros_like(D2[D2i])))
        # print(Tboth.size, ' =? ', Vboth.size, Ref_Race.Distance[corner_index[-1]])
        # Dboth = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(Vboth, Tboth, initial=0))

        # fig5 = plt.figure(5)
        # ax = fig5.add_subplot(1, 1, 1)
        # ax.plot(TT_Race.Distance, v, Dboth, Vboth)
        # plt.xlim(Dboth[0], Dboth[-1])
        # fig5.show()

        # PLOT THIS ONE
        if enable_warnings:
            fig6 = plt.figure(6)
            ax = fig6.add_subplot(1, 1, 1)
            ax.plot(Ref_Race.t, v, Tboth, Vboth, TT_Sim['t'], TT_Sim['v'])  # , Tboth, ToBoth
            plt.xlim(Tboth[0], Tboth[-1])
            fig6.show()
            plt.show()

        if verbosity > 1:
            tgained = Ref_Race.t[int(corner_index[-1])] - Tboth[-1]
            print('Time gained = %.1f s' % tgained, ' on corner %d' % c)

        if TT_Sim['t'].size == 0:  # if not [] == true
            TT_Sim['t'] = np.hstack((TT_Sim['t'], Tboth))
        else:
            TT_Sim['t'] = np.hstack((TT_Sim['t'], Tboth - Tboth[0] + TT_Sim['t'][-1]))
        TT_Sim['v'] = np.hstack((TT_Sim['v'], Vboth))
        TT_Sim['Rpm'] = np.hstack((TT_Sim['Rpm'], RPMBoth))
        TT_Sim['torque'] = np.hstack((TT_Sim['torque'], ToBoth))
        TT_Sim['Distance'] = np.hstack((TT_Sim['Distance'], Dboth))
        TT_Sim['lean'] = np.hstack((TT_Sim['lean'], leanBoth))
        TT_Sim['Id'] = np.hstack((TT_Sim['Id'], Idboth))

        # Now postprocess to find energy
        TT_Sim['t'], indices = np.unique(TT_Sim['t'], return_index=True)  # Remove stops in time
        TT_Sim['v'] = TT_Sim['v'][indices]
        TT_Sim['Rpm'] = TT_Sim['Rpm'][indices]
        TT_Sim['torque'] = TT_Sim['torque'][indices]
        TT_Sim['Distance'] = TT_Sim['Distance'][indices]
        TT_Sim['lean'] = TT_Sim['lean'][indices]
        TT_Sim['Id'] = TT_Sim['Id'][indices]
        # Limit braking torque to rear wheel regenerative
        torque = braking_regen(TT_Sim['Rpm'] / 30 * np.pi * TT_Sim['N'][1] / TT_Sim['N'][0],
                               TT_Sim['torque'] * TT_Sim['N'][0] / TT_Sim['N'][1],
                               TT_Sim['brake']['LimitTorque'], TT_Sim['brake']['k_wt'])
        # this is the desired motor torque
        torque_motor = torque / TT_Sim['N'][0] * TT_Sim['N'][1]
        # todo limit it due to FW
        # TT_Sim['Iq'] = torque_motor / TT_Sim['constants']['Km']
        TT_Sim['Iq'] = motor_currents(TT_Sim, torque_motor)
        TT_Sim['Is'] = np.sqrt(TT_Sim['Id'] ** 2 + TT_Sim['Iq'] ** 2)  # negligeble execution time
        [total_loss, resistive, moving] = motor_losses(TT_Sim['Is'], TT_Sim['Rpm'], TT_Sim['motor']['Rs'],
                                                       TT_Sim['motor']['k_rpm'][2], TT_Sim['motor']['k_rpm'][1])
        TT_Sim['P']['MotorLosses'] = total_loss
        TT_Sim['P']['motRLoss'] = resistive
        TT_Sim['P']['motwLoss'] = moving

        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'],
                                                  TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'], TT_Sim['motor']['Lq'],
                                                  TT_Sim['Id'], TT_Sim['Iq'], TT_Sim['Rpm'] / 30 * np.pi)
        TT_Sim['Vs'] = v_s
        TT_Sim['Vd'] = v_d
        TT_Sim['Vq'] = v_q
        TT_Sim['PF'] = power_factor
        # TT_Sim['Vdc'] = np.hstack((TT_Sim['Vdc'], TT_Sim['Vdc_sim'][c] * np.ones(TT_Sim['t'].size - TT_Sim['Vdc'].size)))

        TT_Sim['P']['Mech'] = TT_Sim['Rpm'] / 30 * np.pi * torque_motor
        TT_Sim['P']['Motor'] = TT_Sim['P']['Mech'] + TT_Sim['P']['MotorLosses']
        Vdcboth = np.hstack((Vdc[Di], Vdc[Di][-1] + np.zeros_like(D2[D2i])))
        TT_Sim['Vdc2'] = np.hstack((TT_Sim['Vdc2'], Vdcboth))
        TT_Sim['Vdc2'] = TT_Sim['Vdc2'][indices]
        TT_Sim['P']['DriveLosses'] = TT_Sim['drive']['n'] * inverter_losses(TT_Sim['Vdc2'], TT_Sim['Vs'],
                                                                            TT_Sim['Is'] / np.sqrt(2) /
                                                                            TT_Sim['drive']['n'], TT_Sim['PF'],
                                                                            TT_Sim['motor']['Lq']*3/2, TT_Sim)[0]
        TT_Sim['P']['Drive'] = TT_Sim['P']['Motor'] + TT_Sim['P']['DriveLosses']
        TT_Sim['Idc'] = TT_Sim['P']['Drive'] / TT_Sim['Vdc2']

        TT_Sim['Energy']['Drive'] = np.trapz(TT_Sim['P']['Drive'], TT_Sim['t'])
        if not battery_fixed:
            TT_Sim['Vdc_sim'][c + 1], battery_status = battery_simple(TT_Sim, TT_Sim['Energy']['Drive'] / 3600, verbosity)
            # TODO battery_simple should really use power as input so it can calculate IR drop losses inside battery
            if battery_status == 'LVC':
                TT_Sim['ERROR']['LVC'] = 1
        if verbosity > 1:
            print('Energy=%.0fWh' % (TT_Sim['Energy']['Drive'] / 3600))
            print('Voltage=%.1fV' % (TT_Sim['Vdc_sim'][c + 1]))
        return TT_Sim


def inverter_losses(v_bus, v_oll, i_o_rms, power_factor, l, sim):
    #def inverter_losses(v_bus, v_oll, i_o_rms, power_factor, l, f_sw, u_ce0=0.8, u_d0=1, r_c=0.95e-3, r_d=0.54e-3,
    #                    e_ton=12e-3, e_toff=25e-3, e_d=9.5e-3, vce_test=300.0, ic_test=550):
    # e.g. [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]
    # = inverter_losses(450,230,350,0.95,75e-6,1e4,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3)
    f_sw = sim['IGBT']['Fsw']
    u_ce0 = sim['IGBT']['Uce0']
    u_d0 = sim['IGBT']['Ud0']
    r_c = sim['IGBT']['Rc']
    r_d = sim['IGBT']['Rd']
    e_ton = sim['IGBT']['Eon']
    e_toff = sim['IGBT']['Eoff']
    e_d = sim['IGBT']['Ed']
    vce_test = sim['IGBT']['Vce_test']
    ic_test = sim['IGBT']['Ic_test']

    i_ripple = (v_bus - 1.41459 * v_oll) * v_oll / (2.0 * l * v_bus * f_sw)
    i_opk = 1.41459 * i_o_rms
    i_opk_sq = i_o_rms ** 2
    m = 0.81650 * v_oll / v_bus  # 0.81650 = sqrt(2/3)
    m_cos_pf = m * np.cos(power_factor)
    p_ct = u_ce0 * i_opk * (1.0 / (2.0 * np.pi) + m_cos_pf / 8.0) + r_c * i_opk_sq * (0.125 + m_cos_pf / (3.0 * np.pi))
    p_cd = u_d0 * i_opk * (1.0 / (2.0 * np.pi) - m_cos_pf / 8.0)  + r_d * i_opk_sq * (0.125 - m_cos_pf / (3.0 * np.pi))
    i_dc = i_opk / np.pi  # DC equivalent to actual AC output current
    i_c_on = i_dc - i_ripple / 2.0
    i_c_off = i_dc + i_ripple / 2.0
    p_st = (e_ton * i_c_on + e_toff * i_c_off) * f_sw / vce_test * v_bus / ic_test
    p_sd = e_d * f_sw / vce_test * v_bus / ic_test * i_dc
    p_total = 6.0 * (p_ct + p_cd + p_st + p_sd)
    return [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]


def power_mech_losses(v, rCdA, m, p_tyre, gradient):
    g = 9.81
    # Losses *add bearing and transmission losses*
    p_air = v / 2.0 * np.square(v) * rCdA
    p_roll = v * 0
    th = 165.0 / 3.6
    p_roll[v < th] = v[v < th] * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v[v < th] * 3.6))
    p_roll[v > th] = v[v > th] * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v[v > th] * 3.6))
    # torque_gradient = torque_motor - torque_roll - torque_air
    p_gradient = v * m * g * gradient
    return [p_air, p_gradient, p_roll]


def power_mech_losses_no_grad(v, rCdA, m, p_tyre):
    g = 9.81
    # Losses *add bearing and transmission losses*
    p_air = v / 2.0 * np.square(v) * rCdA
    p_roll = v * 0
    th = 165.0 / 3.6
    p_roll[v < th] = v[v < th] * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v[v < th] * 3.6))
    p_roll[v > th] = v[v > th] * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v[v > th] * 3.6))
    return [p_air, p_roll]


def v_dq_pmsm(ke, poles, r_s, l_d, l_q, i_d, i_q, w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    w_e = w_m * (poles / 2)
    v_d = r_s * i_d - w_e * l_q * i_q
    v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    v_s = (v_d ** 2 + v_q ** 2) ** 0.5
    #
    v_alpha = np.arctan(v_d / (v_q + 0.0001))
    i_alpha = np.arctan(i_d / (i_q + 0.0001))
    power_factor = i_alpha - v_alpha
    # power_factor[np.isinf(power_factor)] = 1
    return [v_s, v_d, v_q, power_factor]


def motor_losses(stator_current, speed_rpm, resistance=6.333e-3, k_rpm_2=2.6853e-5, k_rpm_1=0.01528):
    resistive = stator_current ** 2 * resistance  # check this, 3*?
    moving = k_rpm_2 * speed_rpm ** 2 + k_rpm_1 * speed_rpm
    total = resistive + moving
    return [total, resistive, moving]


def battery_simple(TT_Sim, energy_used, verbosity):
    ref_cell_energy = TT_Sim['battery']['cell_ve']['e'] * TT_Sim['battery']['cellAh'] * TT_Sim['battery']['cellVnom']
    ref_cell_voltage = TT_Sim['battery']['cell_ve']['v']

    n_cells = TT_Sim['battery']['series'] * TT_Sim['battery']['parallel']
    # print('Battery has', str(n_cells), 'cells')
    energy_remaining = TT_Sim['battery']['E_init'] - energy_used
    energy_remaining_cell = energy_remaining / n_cells
    v_now = TT_Sim['battery']['series'] * np.interp(energy_remaining_cell, ref_cell_energy, ref_cell_voltage, 0, 5)

    status = []
    if v_now < 2.7 * TT_Sim['battery']['series']:
        print('Error: Battery overdischarged, cell V=', str(v_now))
        status = 'LVC'
    if v_now > 4.35 * TT_Sim['battery']['series']:
        print('Error: Battery probably on fire, cell V=', str(v_now))
        status = 'fire'
    else:
        if v_now > 4.2 * TT_Sim['battery']['series']:
            if verbosity > 1:
                print('Warning: Battery overcharged to', str(v_now))
    return [v_now, status]


def battery_heat(IR_cell, m_cell, I_cell, t, verbosity):
    shc = 1014
    P_cell = I_cell ** 2 * IR_cell
    E_cell = integrate.cumtrapz(P_cell, t, initial=0)
    T_cell = E_cell / shc / m_cell
    if verbosity > 0:
        print('Battery temperature rise (end)=%.1f°K' % T_cell[-1])
    return T_cell


def energy_losses(file, name, IR):
    data_file = sio.loadmat(file, struct_as_record=False, squeeze_me=True)
    Data = data_file[name]
    rhoCdA = Data.constants.rho * Data.constants.cd * Data.constants.area
    [Data.P.air, Data.P.roll] = power_mech_losses_no_grad(Data.v, rhoCdA, Data.constants.m, Data.constants.p_tyre)

    try:
        Data.Idc = Data.P.Drive / Data.Vdc
    except:
        Data.Idc = Data.P.Drive / Data.Vdc2
    Data.P.DC_loss = Data.Idc ** 2 * IR
    Data.P.Mech_brake = Data.P.Mech - (Data.torque * Data.Rpm / 30 * np.pi)  # .torque includes system braking torque

    Data.Energy.air = np.trapz(Data.P.air, Data.t)
    # Data.Energy.gradient = np.trapz(Data.P.gradient, Data.t)
    Data.Energy.roll = np.trapz(Data.P.roll, Data.t)
    Data.Energy.brake = np.trapz(Data.P.Mech_brake, Data.t)
    Data.Energy.DriveLosses = np.trapz(Data.P.DriveLosses, Data.t)
    Data.Energy.MotorLosses = np.trapz(Data.P.MotorLosses, Data.t)
    Data.Energy.DC_loss = np.trapz(Data.P.DC_loss, Data.t)
    Data.Energy.Battery = Data.Energy.Drive + Data.Energy.DC_loss

    torque_motor = Data.Iq * Data.constants.Km + 0.01
    e = chain_eff(Data.N[1], Data.N[0], 12.7, 1.21 / 1.27, 1, Data.Rpm / 30 * np.pi,
                        torque_motor, 0.11, 0.00508 * 1.2)[0]
    Data.P.Chain_loss = abs(Data.Rpm / 30 * np.pi * torque_motor * (1 - e))
    Data.Energy.Chain = np.trapz(Data.P.Chain_loss, Data.t)

    return Data


def wheel_forces(x, p, h, b, r, m, force_air, torque_wheel, show_figure):
    # wheel_forces(1.415, 0.6, 0.7, r, cd, area, np.square(v) * rho * cd * area / 2.0, torque_wheel, 1)
    # p=1.415;  % Horizontal distance between wheel centers
    # h=0.6;  % Vertical distance from ground to center of mass
    # b=0.7;  % Horizontal distance from REAR wheel to center of mass
    # force_air = 1 / 2.0 * np.square(v) * rho * cd * area
    g = 9.81
    force_wheel = torque_wheel / r
    force_rear = (m * g * (p - b) + h * (force_wheel + force_air)) / p
    force_front = (m * g * b - h * (force_wheel + force_air)) / p

    mu_tyre = 1.2  # Friction co-efficient of tyre
    torque_MaxFront = force_front * mu_tyre * r
    torque_MaxRear = force_rear * mu_tyre * r

    n = x # np.linspace(1,force_rear.size,force_rear.size)
    if show_figure == 1:
        fig, ax1 = plt.subplots()
        ax1.plot(n, force_front)
        ax1.set_xlabel('samples')
        # Make the y-axis label and tick labels match the line color.
        ax1.set_ylabel('Force (N)', color='b')
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        ax2 = ax1.twinx()
        ax2.plot(n, torque_MaxFront, 'r')
        ax2.set_ylabel('Max torque (Nm)', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        plt.title('Tyre contact forces')
        # plt.show()
        fig.show()

    return [force_rear, force_front, torque_MaxRear, torque_MaxFront]


def wheelie_torque_limit(TT_Sim, R, V):
    return TT_Sim['N'][1] / TT_Sim['N'][0] * R * (TT_Sim['constants']['m'] * 9.81 *
                                                    TT_Sim['constants']['b'] / TT_Sim['constants']['h']
                                                    - V ** 2 * TT_Sim['constants']['rho'] * TT_Sim['constants']['cd'] *
                                                    TT_Sim['constants']['area'] / 2)


def set_speed_limit(TT_Sim, course_speed_limit):
    w_max = course_speed_limit * TT_Sim['N'][0] / TT_Sim['N'][1] / TT_Sim['constants']['r']
    TT_Sim['motor']['W_speed_lim'] = 0.98 * w_max
    TT_Sim['motor']['W_lim'] = TT_Sim['motor']['W_speed_lim'] * 1.02  # 8400 / 30 * np.pi
    return TT_Sim


def v_tyre_load_sens(v_ref, m_ref, k, m):
    v = v_ref * np.power(m_ref, (1 - k) / 2) * np.power(m, (k - 1) / 2)
    return v


def motor_saturation_coeff(kt, i_pk, t_pk, x2, y2):
    # Call with rms values, returns peak values
    # T = a*I.^3+b*I.^2+kt*I+0;
    # i_pk = 987.1
    # t_pk = 308.64
    # kt = 0.366
    # x2 = 493.6  # stall current continuous
    # y2 = 178.47  # stall torque continuous

    b = (y2 + x2 ** 3 / i_pk ** 3 * (kt * i_pk - t_pk) - kt * x2) / (x2 ** 2 - x2 ** 3 / i_pk)
    a = 1.0 / i_pk ** 3 * (t_pk - kt * i_pk - b * i_pk ** 2)

    a /= 2 * 2 ** 0.5
    b /= 2
    c = kt / (2 ** 0.5)

    return [a, b, c, 0]


def motor_saturation_interp(sim):
    sim['motor']['I_interp'] = np.linspace(sim['motor']['I_con'], 1.3 * sim['motor']['I_pk'], 25)
    sim['motor']['T_interp'] = np.polyval(sim['motor']['co'], sim['motor']['I_interp'])
    # sim['motor']['T_interp'] = np.array([motor_torque(sim['motor']['co'], x) for x in sim['motor']['I_interp']])
    # sim['motor']['I_interp'] = np.append(sim['motor']['I_interp'], 10 * sim['motor']['I_pk'])
    # sim['motor']['T_interp'] = np.append(sim['motor']['T_interp'], 1.1 * sim['motor']['T_interp'][-1])
    A = sim['motor']['I_interp'][::-1]  # [:0:-1]
    B = sim['motor']['T_interp'][::-1]  # [:0:-1]
    sim['motor']['I_interp'] = np.concatenate([-1.0 * A, sim['motor']['I_interp']])
    sim['motor']['T_interp'] = np.concatenate([-1.0 * B, sim['motor']['T_interp']])
    return sim


def motor_torque(co, i):
    return co[0] * i ** 3 + co[1] * i ** 2 + co[2] * i + co[3]


def motor_current(p, torque):
    print('motor_current() may not work in some cases')
    if torque > 0:
        roots = sorted(np.roots(p[:3] + [-1.0 * torque]))[1]
    else:
        roots = sorted(np.roots(p[:3] + [1.0 * torque]))[1]
    return roots


def motor_currents(sim, torque):
    return np.interp(torque, sim['motor']['T_interp'], sim['motor']['I_interp'])


def motor_current_newton(co, torque, rated_torque, error):
    i = torque / co[2]
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
            if count > 1000:
                print('NEWTON ITERATION STUCK ' + str(torque) + ' ' + str(rated_torque))
                iterate = False
    return i


def motor_sizing(sim):
    if sim['motor']['manufacturer'] == 'Parker':
        turns = sim['motor']['N']
        l_core = sim['motor']['L_core']
        sim['motor']['poles'] = 12

        ref_cores = np.array([50, 100, 150, 175, 200, 300, 400])
        ref_ipk = [196, 422, 697, 978.1, 614, 651, 572]
        ref_icon = [98, 211, 348.5, 493.6, 307, 326, 286]
        ref_tcon = np.array([44, 95, 146.667, 178.47, 173, 298, 404])
        ref_tpk = np.array([79, 168, 256.96, 308.64, 314, 520, 700])
        ti = (ref_cores - 175) * np.pi * (7850 - 4450) * 0.0304 ** 2 * 0.001 + (7.85/4.45*3.618-3.618)
        ref_m = np.array([25, 35, 46, 51, 56, 77, 97]) - ti
        ref_j = [0.00824, 0.01578, 0.02332, 0.0253, 0.03086, 0.04609, 0.06117]
        ref_kt = ref_tcon / ref_icon
        ref_ke = ref_kt * (2/3.) ** 0.5
        ref_turns = np.round((ref_ke / 0.3469 * 150 / ref_cores * 11.5) * 2) / 2
        ref_turn = np.interp(l_core, ref_cores, ref_turns)

        sim['constants']['Km'] = np.interp(l_core, ref_cores, ref_kt) * turns / ref_turn

        sim['motor']['i_rms_pk'] = np.interp(l_core, ref_cores, ref_ipk) * ref_turn / turns
        sim['motor']['i_rms_con'] = np.interp(l_core, ref_cores, ref_icon) * ref_turn / turns
        sim['motor']['T_pk'] = np.interp(l_core, ref_cores, ref_tpk)
        sim['motor']['T_con'] = np.interp(l_core, ref_cores, ref_tcon)

        sim['motor']['m'] = np.interp(l_core, ref_cores, ref_m)
        sim['motor']['J'] = np.interp(l_core, ref_cores, ref_j)

        sim['motor']['Ld'] = 53e-6 * (turns / 11.5) ** 2 * l_core / 150
        sim['motor']['Lq'] = 61e-6 * (turns / 11.5) ** 2 * l_core / 150
        r_hot = 0.007313 * turns / 11.5  # seems 2 be 'cold' (60C) resistance at rated winding temperature (140C), V6
        r_end_turns = r_hot * 43.688 / (150 + 43.688)
        r_mid_turns = r_hot * 150.0 / (150 + 43.688)
        sim['motor']['Rs'] = r_end_turns + r_mid_turns * l_core / 150

        sim['motor']['k_rpm'] = [0, 0.01528 / 150 * l_core, 2.6853e-5 / 150 * l_core]

        sim['motor']['W_max'] = 11000 * np.pi / 30.0  # 5500 rpm no load, 4700rpm full load 680V

    if sim['motor']['manufacturer'] == 'me':
        sim['motor']['i_rms_pk'] = 490.8
        sim['motor']['i_rms_con'] = 245.4
        sim['motor']['T_pk'] = 300
        sim['motor']['T_con'] = 164.3
        sim['motor']['m'] = 22.26
        sim['motor']['J'] = 0.0285
        sim['motor']['N'] = 11
        # sim['motor']['L_core'] = 100  # USED FOR moving LOSS CALCULATION SCALED FROM PARKER 150mm
        sim['motor']['Ld'] = 74e-6
        sim['motor']['Lq'] = 77e-6
        sim['motor']['Rs'] = 0.0092  # check temperature
        sim['motor']['poles'] = 14
        sim['motor']['W_max'] = 11000 * np.pi / 30.0  # 5500 rpm no load, 4700rpm full load 680V
        sim['motor']['k_rpm'] = [10.37, 0.08, 5.378e-5]
        sim['constants']['Km'] = sim['motor']['T_con'] / sim['motor']['i_rms_con']  # Nm/Arms cont
    if sim['motor']['manufacturer'] == 'Emrax':
        sim['motor']['i_rms_pk'] = 400.0
        sim['motor']['i_rms_con'] = 190.0
        sim['motor']['T_pk'] = 500.0
        sim['motor']['T_con'] = 250.0
        sim['motor']['m'] = 20.3
        sim['motor']['J'] = 0.093
        sim['motor']['N'] = 4
        sim['motor']['L_core'] = 150  # USED FOR moving LOSS CALCULATION SCALED FROM PARKER 150mm
        sim['motor']['Ld'] = 126e-6
        sim['motor']['Lq'] = 118e-6
        sim['motor']['Rs'] = 0.0115*1.373  # Internal phase resistance at 25C = 0.0115, 120C = 0.393%*95 more
        sim['motor']['poles'] = 20.0
        sim['motor']['W_max'] = 5500 * np.pi / 30.0  # 5500 rpm no load, 4700rpm full load 680V
        sim['constants']['Km'] = sim['motor']['T_con'] / sim['motor']['i_rms_con']  # Nm/Arms cont
    if sim['motor']['manufacturer'] == 'import':
        # load from m file
        # Rs at rated temperature
        restore = False
        if 'T_max' in sim['motor']:
            T_max_temp = sim['motor']['T_max']
            P_max_temp = sim['motor']['P_max']
            W_max_temp = sim['motor']['W_max']
            W_lim_temp = sim['motor']['W_lim']
            W_speed_lim = sim['motor']['W_speed_lim']
            restore = True
        motor_import = loadmat('data_import/' + sim['file']['motorimport'])

        imp_turns = motor_import['motor']['N']
        act_turns = sim['motor']['N']
        # imp Ld,Lq,co

        sim['motor'] = motor_import['motor']
        sim['motor']['N'] = act_turns
        print('NEED TO ADJUST Rs, for imported motor')
        sim['motor']['Ld'] *= (act_turns / imp_turns) ** 2
        sim['motor']['Lq'] *= (act_turns / imp_turns) ** 2
        sim['motor']['co'] *= act_turns / imp_turns
        sim['motor']['i_rms_con'] *= imp_turns / act_turns
        sim['motor']['i_rms_pk'] *= imp_turns / act_turns
        sim['motor']['Rs'] *= imp_turns / act_turns # r_end_turns + r_mid_turns * l_core / 150

        print('check int division?')

        sim['motor']['manufacturer'] = 'import'  # just to make sure, else successive calls will fail
        if 'co' not in sim['motor']:
            print('ERROR, motor saturation co-efficients missing')
        if restore:
            sim['motor']['T_max'] = T_max_temp
            sim['motor']['P_max'] = P_max_temp
            sim['motor']['W_max'] = W_max_temp
            sim['motor']['W_lim'] = W_lim_temp
            sim['motor']['W_speed_lim'] = W_speed_lim
        else:
            sim['motor']['W_max'] = 11e3

        sim['constants']['Km'] = sim['motor']['T_con'] / sim['motor']['i_rms_con']  # Nm/Arms cont
        # sim['IGBT']['Fsw'] = 1000 * sim['motor']['poles']

    sim['motor']['Ke'] = sim['constants']['Km'] * (2.0 / 3.) ** 0.5 / (sim['motor']['poles'] / 2.0)
    sim['J']['motor'] = sim['motor']['J']

    if sim['motor']['manufacturer'] != 'import':
        sim['motor']['co'] = motor_saturation_coeff(sim['constants']['Km'], sim['motor']['i_rms_pk'],
                                                    sim['motor']['T_pk'], sim['motor']['i_rms_con'],
                                                    sim['motor']['T_con'])

    sim['motor']['I_pk'] = sim['motor']['i_rms_pk'] * 2 ** 0.5
    sim['motor']['I_con'] = sim['motor']['i_rms_con'] * 2 ** 0.5

    v_s_max = v_dq_pmsm(sim['motor']['Ke'], sim['motor']['poles'], sim['motor']['Rs'], sim['motor']['Ld'],
                    sim['motor']['Lq'], 0, sim['motor']['I_pk'], sim['motor']['W_max'])[0]
    #  print('rpmmax ' + str(sim['motor']['W_max']*30.0/3.14159))
    #  print('vsmax ' + str(v_s_max))
    #  sim['battery']['V_max'] = 680
    if v_s_max > 0.9 * sim['battery']['V_max']:
        wfw = np.linspace(0, sim['motor']['W_max'], 1000)
        v_s = v_dq_pmsm(sim['motor']['Ke'], sim['motor']['poles'], sim['motor']['Rs'], sim['motor']['Ld'],
                        sim['motor']['Lq'], 0, sim['motor']['I_pk'], wfw)[0]
        wfwlimited = wfw[v_s > 0.9 * sim['battery']['V_max']]
        sim['motor']['P_pk'] = sim['motor']['T_pk'] * wfwlimited[0]
    else:
        sim['motor']['P_pk'] = sim['motor']['T_pk'] * sim['motor']['W_max']

    sim['motor']['name'] = str(sim['motor']['manufacturer']) + '_GVM210-' + str(sim['motor']['L_core']) + '-N' + str(sim['motor']['N'])
    return sim


def charge_battery(sim, ratio):
    sim['battery']['E_rated'] = sim['battery']['series'] * sim['battery']['parallel'] * \
                                   sim['battery']['cellAh'] * sim['battery']['cellVnom']
    sim['battery']['E_init'] = sim['battery']['E_rated'] * ratio
    return sim


def weigh(sim):
    m_orig_bike = 235           # 235kg at Assen 2016 UoN02
    m_orig_battery = 4.8 * 18   # 18 of 4.8kg batteries at Assen
    m_orig_tray = 4.75          # 4.75kg for 02 main tray ONLY
    m_orig_motor = 47
    m_orig_drive = 5
    m_orig_parts = m_orig_bike - m_orig_battery - m_orig_tray - m_orig_motor - m_orig_drive

    m_cell = sim['battery']['cellAh'] * sim['battery']['cellVnom'] / sim['battery']['E_density']
    sim['mass'] = {'rider': 90}
    sim['mass']['motor'] = sim['motor']['m']
    sim['mass']['battery'] = m_cell * sim['battery']['series'] * sim['battery']['parallel']
    sim['mass']['drives'] = sim['drive']['n'] * sim['drive']['m']
    sim['mass']['tray'] = m_orig_tray * sim['mass']['battery'] / m_orig_battery

    sim['mass']['bike'] = m_orig_parts + sim['mass']['motor'] + sim['mass']['battery'] + sim['mass']['drives'] + sim['mass']['tray']
    sim['constants']['m'] = sim['mass']['rider'] + sim['mass']['bike']
    return sim


def scrutineering(sim, charge_ratio):
    sim['scrutineering']['score'] = 255
    sim = charge_battery(sim, charge_ratio)
    sim = weigh(sim)

    # Electrical scrutineering
    if motor_current_newton(sim['motor']['co'], sim['motor']['T_max'], 0, 0.01) / sim['drive']['n'] < sim['drive']['I_max']:
        #print(motor_current_newton(sim['motor']['co'], sim['motor']['T_max'], 0, 0.01))
        # Check if motor current required for specified torque is within drive rating
        #print('N= ' + str(sim['motor']['N']) + str(sim['drive']['n']) + ' Idrive=' + str(motor_current(sim['motor']['co'], sim['motor']['T_max']) / sim['drive']['n']))
        sim['scrutineering']['score'] -= 1
    if sim['motor']['T_max'] < (sim['motor']['T_pk'] * 1.05):
        sim['scrutineering']['score'] -= 2
    if sim['motor']['P_max'] <= sim['motor']['P_pk']:
        # print(sim['motor']['P_max'], sim['motor']['P_pk'])
        sim['scrutineering']['score'] -= 4
    if (sim['battery']['V_max']) < sim['scrutineering']['volt_limit']:
        sim['scrutineering']['score'] -= 8

    to_max = C_torque_limits(sim['constants']['r'], 0, sim['N'][1], sim['N'][0], sim['constants']['m'],
                             sim['constants']['b'], sim['constants']['h'], sim['constants']['rho'],
                             sim['constants']['cd'], sim['constants']['area'],
                             sim['constants']['mu_tyre'])[0]
    # if (sim['motor']['T_max'] < 1.1 * to_max) and (sim['motor']['T_max'] > 0.7 * to_max):
    if sim['motor']['T_max'] > 0.7 * to_max:
        sim['scrutineering']['score'] -= 32
    # print(sim['motor']['T_max'], 1.1 * to_max, 0.7 * to_max)
    to_max = C_torque_limits(sim['constants']['r'], sim['v_max'], sim['N'][1], sim['N'][0], sim['constants']['m'],
                             sim['constants']['b'], sim['constants']['h'], sim['constants']['rho'],
                             sim['constants']['cd'], sim['constants']['area'],
                             sim['constants']['mu_tyre'])[0]
    w_max = sim['v_max'] / sim['constants']['r'] / sim['N'][1] * sim['N'][0]
    p_max = w_max * to_max  # (sim['motor']['W_lim'] + sim['motor']['W_speed_lim']) * 0.5 * to_max
    if sim['motor']['P_max'] < p_max:  # Power = [100e3, wheelie lim]
        sim['scrutineering']['score'] -= 64
    # print(sim['motor']['P_max'], p_max, to_max)


    # Mechanical scrutineering
    if (sim['constants']['m'] - sim['mass']['rider']) < sim['scrutineering']['weight_limit']:
        sim['scrutineering']['score'] -= 16

    if sim['motor']['W_lim'] < sim['motor']['W_max']:
        sim['scrutineering']['score'] -= 128

    if sim['scrutineering']['score'] == 0:
        sim['scrutineering']['passed'] = True
    else:
        sim['scrutineering']['passed'] = False
    return sim


def w_fw(w, TT_Sim, v_dc, i_q):
    return vs_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'], TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'],
                   TT_Sim['motor']['Lq'], 0, i_q, w) - (v_dc * 0.9 * 0.866)


def id_fw(i_d, TT_Sim, v_dc, w, i_q):
    return vs_pmsm(TT_Sim['motor']['Ke'], TT_Sim['motor']['poles'], TT_Sim['motor']['Rs'], TT_Sim['motor']['Ld'],
                   TT_Sim['motor']['Lq'], i_d, i_q, w) - (v_dc * 0.9 * 0.866)


def iq_pmsm(ke, poles, r_s, l_d, l_q, i_d, w_m, v_s_max):
    w_e = w_m * (poles / 2)
    a = w_e ** 2 * l_q ** 2 + r_s ** 2
    b = 2 * r_s * w_e * (i_d * (l_d - l_d) + ke)
    c = r_s ** 2 * i_d ** 2 + w_e ** 2 * (l_d * i_d + ke) ** 2 - v_s_max ** 2
    return (-b + (b ** 2 - 4 * a * c) ** 0.5) / 2 / a


def id_pmsm(ke, poles, r_s, l_d, l_q, i_q, w_m, v_s_max):
    w_e = w_m * (poles / 2)
    a = w_e ** 2 * l_d ** 2 + r_s ** 2
    b = 2 * w_e * (w_e * l_d * ke + r_s * i_q * (l_d - l_q))
    c = (r_s * i_q + w_e * ke) ** 2 + w_e ** 2 * l_q ** 2 * i_q ** 2 - v_s_max ** 2
    return (-b + (b ** 2 - 4 * a * c) ** 0.5) / 2 / a


def w_pmsm(ke, poles, r_s, l_d, l_q, i_d, i_q, v_s_max):
    a = l_q ** 2 * i_q ** 2 + (l_d * i_d + ke) ** 2
    b = 2 * r_s * (i_d * i_q * (l_d - l_q) + i_q * ke)
    c = r_s ** 2 * (i_d ** 2 + i_q ** 2) - v_s_max ** 2
    w_e = (-b + (b ** 2 - 4 * a * c) ** 0.5) / 2 / a
    return 2.0 / poles * w_e


def vs_pmsm(ke, poles, r_s, l_d, l_q, i_d, i_q, w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    w_e = w_m * (poles / 2)
    v_d = r_s * i_d - w_e * l_q * i_q
    v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    v_s = (v_d ** 2 + v_q ** 2) ** 0.5
    return v_s


def torque_fw(TT_Sim, t_motor, t_motor_prev, w_m_now, vdc, voc):
    # t_motor_prev = t_motor
    Ke = TT_Sim['motor']['Ke']
    Poles = TT_Sim['motor']['poles']
    Rs = TT_Sim['motor']['Rs']
    Ld = TT_Sim['motor']['Ld']
    Lq = TT_Sim['motor']['Lq']
    IR = TT_Sim['battery']['IR']
    L_core = TT_Sim['motor']['L_core']
    T_con = TT_Sim['motor']['T_con']
    driven = TT_Sim['drive']['n']
    co = TT_Sim['motor']['co']
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
        motor_loss_now = C_motor_losses(is_now, w_m_now * 30 / np.pi, Rs, TT_Sim['motor']['k_rpm'][2],
                                        TT_Sim['motor']['k_rpm'][1])[0]
        [vs, vd, vq, PF] = C_v_dq_pmsm_array(Ke, Poles, Rs, Ld, Lq, id_now, iq_now, w_m_now)
        p_drive_loss_now = driven * C_inverter_loss(vdc, vs, is_now / np.sqrt(2) / driven, PF, 82e-6, 13e3, 0.8, 1,
                                                    0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
        # print(str(p_batt_now / (t_motor * w_m_now + motor_loss_now + p_drive_loss_now)))
        p_batt_now = t_motor * w_m_now + motor_loss_now + p_drive_loss_now
        vdc = voc - p_batt_now / vdc * IR
        id = C_id_pmsm(Ke, Poles, Rs,
                       Ld, Lq, iq_now, w_m_now, vdc * 0.9 * 0.866)
    else:
        id = 0
    return [t_motor, id, vdc]


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


