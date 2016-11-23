import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from os import remove
from scipy.integrate import odeint
from scipy import integrate


def motorbike_mech(v, t, r, rho, cd, jr, area, m, p_tyre, rpm, t_motor, n2, n1, gradient, gradient_t):
    motor_rpm = v / (1 / 60 * n2 / n1 * 2 * np.pi * r)
    motor_torque = np.interp(motor_rpm, rpm, t_motor, 100000.0, 100000.0)
    # if(empty(e_chain))
    e_chain = chain_eff(n2, n1, 12.7, 1.21 / 1.27, 1, motor_rpm / 30 * np.pi, motor_torque, 0.11, 0.00508 * 1.2)[
        0]  # Bush diameter guessed
    # end
    torque = e_chain * n1 / n2 * motor_torque  # v=rpm*1/60*n2/n1*2*np.pi*r
    g = 9.81
    # Losses *add bearing and transmission losses*
    torque_air = r / 2.0 * np.square(v) * rho * cd * area
    if v < (165.0 / 3.6):  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = r * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v * 3.6))
    else:
        torque_roll = r * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v * 3.6))

    torque_gradient = r * m * g * np.interp(t, gradient_t, gradient, -1, -1)  # -1 so bike falls off cliff
    # dydt = (T(t)-Tr-Ta)/Jr
    #  Evaluate ODE at time t
    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt


def motorbike_mech2(v, t, r, rho, cd, jr, area, m, p_tyre, t_mot, t_mott, n2, n1, e_chain, gradient, gradient_t):
    torque = e_chain * n1 / n2 * np.interp(t, t_mott, t_mot, 0, 0)  # v=rpm*1/60*n2/n1*2*np.pi*r
    g = 9.81
    # Losses *add bearing and transmission losses*
    torque_air = r / 2.0 * np.square(v) * rho * cd * area
    if v < (165.0 / 3.6):  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = r * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v * 3.6))
    else:
        torque_roll = r * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v * 3.6))

    torque_gradient = r * m * g * np.interp(t, gradient_t, gradient, -1, -1)  # -1 so bike falls off cliff
    # dydt = (T(t)-Tr-Ta)/Jr
    #  Evaluate ODE at time t
    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt


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
    if f_chaintension.max() > 20000.0:
        print('Warning, chain tension very high')

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
    return [e, speed, f_chaintension]


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


def lap_analyse(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                structure_map, var_name_brake, rated_energy, initial_energy, n_series):
    # Import course map and reference lap data
    mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
    TT_map = mat_contents[structure_map]
    corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
    locsmin = 1 + corners[var_name_brake]  # -1 as matlab indexing starts at 1

    # Preallocate model bike data structure
    TT_Sim['t'] = np.array([])
    TT_Sim['v'] = np.array([])
    TT_Sim['Iq'] = np.array([])
    TT_Sim['torque'] = np.array([])
    TT_Sim['Distance'] = np.array([])
    TT_Sim['Vdc_sim'] = np.zeros(locsmin.size + 1)  # Will store battery V at START of each corner
    TT_Sim['Vdc_sim'][first_corner] = battery_simple(0, rated_energy, initial_energy, n_series)
    TT_Sim['Vdc'] = np.array([])
    TT_Sim['P'] = {}
    TT_Sim['P']['Mech'] = np.array([])
    TT_Sim['P']['Motor'] = np.array([])
    TT_Sim['P']['Drive'] = np.array([])
    TT_Sim['P']['MotorLosses'] = np.array([])
    TT_Sim['P']['DriveLosses'] = np.array([])
    TT_Sim['P']['motRLoss'] = np.array([])
    TT_Sim['P']['motwLoss'] = np.array([])
    TT_Sim['PF'] = np.array([])
    TT_Sim['Energy'] = np.array([])
    TT_Sim['Vd'] = np.array([])
    TT_Sim['Vq'] = np.array([])
    TT_Sim['Vs'] = np.array([])

    TT_Sim['J']['linear'] = TT_Sim['constants']['m'] * TT_Sim['constants']['r'] ** 2
    TT_Sim['J']['r'] = 2 * TT_Sim['J']['wheel'] + TT_Sim['J']['linear'] + np.square(TT_Sim['N'][0] / TT_Sim['N'][1]) * \
                                                                          TT_Sim['J']['motor']  # J referred to wheel

    # Bit of a kludge to convert dict type to mat_structure
    # THIS SHOULD BE ELIMINATED BY USING DICTS THROUGHOUT
    sio.savemat('temp', {'TT_Sim': TT_Sim}, oned_as='column')
    TT_Sim = sio.loadmat('temp', struct_as_record=False, squeeze_me=True)['TT_Sim']
    remove('temp.mat')

    accumulated_energy = 0
    v_flag = -1

    for c in range(first_corner, min(locsmin.size - 1, last_corner)):
        corner_index = np.arange(locsmin[c], locsmin[c + 1] + 100)
        corner_index_end = locsmin[c + 1]

        if v_flag == -1:  # Set start vel to corner vel iff bike made that speed on previous
            v0 = v[corner_index[0]]
        else:
            v0 = v_flag

        # t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]], corner_index.size)
        [t, dt_a] = np.linspace(Ref_Race.t[corner_index_end], Ref_Race.t[corner_index[0]], corner_index.size, retstep=True)
        print('dtA = ', str(dt_a))
        t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[0]]-corner_index.size*dt_a,
                        corner_index.size)
        print('dtA = ', str(t[1] - t[0]))
        gradientt = t
        gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
                             0)  # Initially assume gradient(t) = same as lap data

        wfw = np.linspace(0, TT_Sim.motor.W_max, 1000)
        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(0.34688 / 6, 12, 0.007313, 53e-6, 61e-6, 0,
                                                  TT_Sim.motor.T_max / TT_Sim.constants.Km, wfw)
        wfwlimited = wfw[v_s > 0.9 * TT_Sim.Vdc_sim[c]]
        if wfwlimited.size == 0:
            p_weak = TT_Sim.motor.P_max
        else:
            p_weak = TT_Sim.motor.T_max * wfwlimited[0]

        if p_weak >= TT_Sim.motor.P_max:
            p_weak = TT_Sim.motor.P_max
        else:
            print('Power limited to', str(p_weak), 'by field weakening')

        TT_Sim.motor.w, TT_Sim.motor.t, TT_Sim.motor.p = motor_torque_speed(TT_Sim.motor.T_max, TT_Sim.motor.W_max,
                                                                            p_weak, TT_Sim.motor.W_lim, 50, 0)

        for a in range(0, 2):
            V = odeint(motorbike_mech, v0, t,
                       args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                             TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre,
                             TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t, TT_Sim.N[1], TT_Sim.N[0], gradient,
                             gradientt), hmin=1e-9)
            V = np.squeeze(V)
            D = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(V, t, initial=0))
            gradient = np.interp(D, TT_map.dist, TT_map.gradient, 0.0, 0.0)
            # print(rep)
            # print(D.shape)
            # print(V.shape)
            # print(TT_Race.Distance.shape)
            # print(v.shape)
            # print(TT_map.dist.shape)
            # print(TT_map.gradient.shape)
        #fig2 = plt.figure(2)
        #ax = fig2.add_subplot(1,1,1)
        #ax.plot(D, V, Ref_Race.Distance, v, TT_map.dist, TT_map.gradient * 100)
        #plt.plot(gradientt, gradient)  # t[corner_index], v[corner_index], t, Ref_Race.v,
        #plt.xlabel('Distance')
        #plt.ylabel('v')
        #plt.xlim(D[0], D[-1])
        #fig2.show()
        #plt.show()

        MOTORSPEED = V / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * TT_Sim.constants.r)  # in rpm
        MOTORTORQUE = np.interp(MOTORSPEED, TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t)

        T = t
        t = np.linspace(Ref_Race.t[corner_index_end], Ref_Race.t[corner_index[0]],
                        corner_index.size)  # as previous but flipped
        print('dtB = ', str(t[1] - t[0]))

        TBrake_t = t
        TBrake = -TT_Sim.brake.PeakTorque * TT_Sim.N[1] / TT_Sim.N[0] * np.ones(t.shape)
        rampIndex = TBrake_t > (TBrake_t[0] - TT_Sim.brake.RampTime)
        TBrake[rampIndex] = np.linspace(0, TBrake[-1], sum(rampIndex))
        TBrake_t = np.flipud(TBrake_t)
        TBrake = np.flipud(TBrake)

        # plt.close()
        # plt.plot(Tbraket, Tbrake)
        # plt.show()

        v0 = v[corner_index_end]

        gradientt = t
        gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
                             0)  # Initially gradient(t) = same as lap data
        gradientt = np.flipud(gradientt)

        for a in range(0, 2):
            e_chain = 1
            V2 = odeint(motorbike_mech2, v0, t,
                        args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                              TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
                              TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt))
            V2 = np.squeeze(V2)
            D2 = np.squeeze(Ref_Race.Distance[corner_index_end] + integrate.cumtrapz(V2, t, initial=0))
            gradient = np.interp(D2, TT_map.dist, TT_map.gradient, 0, 0)  # '-5.0, -5.0)

            # graidient goes out of range at end of track!!! should make circular array

        # T2 = np.flipud(t)
        D2 = np.flipud(D2)
        V2 = np.flipud(V2)
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
            print('No. of intersections = %d' % dout.size)

        if dout.size == 0:
            print('Bike too slow on corner %d' % c, ', perhaps you used an Agni motor?')
            print('Higher motor torque required to achieve desired corner speed')
            plt.close()
            fig4 = plt.figure(4)
            ax = fig4.add_subplot(1, 1, 1)
            #ax.plot(T, V, '-o', Ref_Race.t, v, T2, V2, '-o', gradientt, gradient * 100)
            ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o')
            plt.xlim(D[0], D[-1])
            plt.ylim(0, V.max() * 1.2)
            # plt.xlim(TT_Race.t[corner_index[0]], TT_Race.t[corner_index[-1]])
            ax.plot(D[dout], V[dout], 'ro')
            fig4.show()
            plt.show()
            dout = Ref_Race.Distance[corner_index_end]  # BIT BAD THIS - causes jump in Vel
            v_flag = V[(D < dout)]
            v_flag = v_flag[-1]
            print('End vel = %.1f' % v_flag)
        else:
            dout = np.median([D[dout]])
            v_flag = -1

        # print('dout =', str(dout))
        # print('Braking point = %d' % dout, 'm')
        D2i = np.squeeze((D2 > dout) & (D2 < 60725))
        Di = np.squeeze(D <= dout)

        #  Vboth = np.hstack((V[:dout], V2i[dout:]))

        Vboth = np.hstack((V[Di], V2[D2i]))
        Dboth = np.hstack((D[Di], D2[D2i]))
        ToBoth = np.hstack((MOTORTORQUE[Di], TBrake[D2i]))

        # vals = V2i != -1
        # Dboth = D[vals]
        # Vboth = Vboth[vals]

        dt = T[1] - T[0]
        Tboth = T[0] + np.arange(0, Vboth.size * dt, dt)

        # fig5 = plt.figure(5)
        # ax = fig5.add_subplot(1, 1, 1)
        # ax.plot(TT_Race.Distance, v, Dboth, Vboth)
        # plt.xlim(Dboth[0], Dboth[-1])
        # fig5.show()

        # PLOT THIS ONE
        # fig6 = plt.figure(6)
        # ax = fig6.add_subplot(1, 1, 1)
        # ax.plot(Ref_Race.t, v, Tboth, Vboth, TT_Sim.t, TT_Sim.v) # , Tboth, ToBoth
        # plt.xlim(Tboth[0], Tboth[-1])
        # fig6.show()
        # plt.show()

        tgained = Ref_Race.t[corner_index_end] - Tboth[-1]
        print('Time gained = %.1f s' % tgained, ' on corner %d' % c)

        if TT_Sim.t.size == 0:  # if not [] == true
            TT_Sim.t = np.hstack((TT_Sim.t, Tboth))
        else:
            TT_Sim.t = np.hstack((TT_Sim.t, Tboth - Tboth[0] + TT_Sim.t[-1]))
        TT_Sim.v = np.hstack((TT_Sim.v, Vboth))
        TT_Sim.torque = np.hstack((TT_Sim.torque, ToBoth))
        TT_Sim.Distance = np.hstack((TT_Sim.Distance, Dboth))

        # Now postprocess to find energy
        TT_Sim.t, indices = np.unique(TT_Sim.t, return_index=True)  # Remove stops in time
        TT_Sim.v = TT_Sim.v[indices]
        TT_Sim.torque = TT_Sim.torque[indices]
        TT_Sim.Distance = TT_Sim.Distance[indices]

        # Limit braking torque to rear wheel regenerative
        torque = braking_regen(TT_Sim.v / TT_Sim.constants.r, TT_Sim.torque * TT_Sim.N[0] / TT_Sim.N[1],
                               TT_Sim.brake.LimitTorque, TT_Sim.brake.k_wt)  # limits for traction motor braking
        # TT_Sim.Iq = np.hstack((TT_Sim.Iq, torque / TT_Sim.N[0] * TT_Sim.N[1] / TT_Sim.constants.Km))
        TT_Sim.Iq = torque / TT_Sim.N[0] * TT_Sim.N[1] / TT_Sim.constants.Km
        [total_loss, resistive, moving] = motor_losses(TT_Sim.Iq, TT_Sim.v / TT_Sim.constants.r / np.pi * 30 *
                                                       TT_Sim.N[0] / TT_Sim.N[1], 6.333e-3, 2.6853e-5, 0.01528)
        TT_Sim.P.MotorLosses = moving  # np.hstack((TT_Sim.P.MotorLosses, moving))
        TT_Sim.P.motRLoss = resistive  # np.hstack((TT_Sim.P.motRLoss, resistive))
        TT_Sim.P.motwLoss = total_loss  # np.hstack((TT_Sim.P.motwLoss, total_loss))
        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(0.34688 / 6, 12, 0.007313, 53e-6, 61e-6, 0, TT_Sim.Iq,
                                                  TT_Sim.v / TT_Sim.constants.r * TT_Sim.N[0] /
                                                  TT_Sim.N[1])
        TT_Sim.Vs = v_s  # np.hstack((TT_Sim.Vs, v_s))
        TT_Sim.Vd = v_d  # np.hstack((TT_Sim.Vd, v_d))
        TT_Sim.Vq = v_q  # np.hstack((TT_Sim.Vq, v_q))
        TT_Sim.PF = power_factor  # .hstack((TT_Sim.PF, power_factor))
        TT_Sim.Vdc = np.hstack((TT_Sim.Vdc, TT_Sim.Vdc_sim[c] * np.ones(TT_Sim.t.size - TT_Sim.Vdc.size)))
        p_drive_loss = inverter_losses(TT_Sim.Vdc, TT_Sim.Vs, abs(TT_Sim.Iq) / np.sqrt(2), TT_Sim.PF, 82e-6, 13e3, 0.8,
                                       1, 0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
        # TT_Sim.P.DriveLosses = np.hstack((TT_Sim.P.DriveLosses, p_drive_loss))
        TT_Sim.P.DriveLosses = p_drive_loss
        TT_Sim.P.Mech = TT_Sim.v / TT_Sim.constants.r * TT_Sim.Iq * TT_Sim.constants.Km * TT_Sim.N[0] / TT_Sim.N[1]
        TT_Sim.P.Motor = TT_Sim.P.Mech + TT_Sim.P.MotorLosses
        TT_Sim.P.Drive = TT_Sim.P.Motor + TT_Sim.P.DriveLosses
        # TT_Sim.Idc = TT_Sim.P.Drive / TT_Sim.Vdc
        TT_Sim.Energy = np.trapz(TT_Sim.P.Drive, TT_Sim.t)
        TT_Sim.Vdc_sim[c + 1] = battery_simple(TT_Sim.Energy / 3600, rated_energy, initial_energy, n_series)
        print('Energy=', str(TT_Sim.Energy / 3600), 'Wh')
        print('Voltage=', str(TT_Sim.Vdc_sim[c + 1]), 'V')

    # Limit to one lap (and remove some crazy points)
    # indices = TT_Sim_Distance < 60725

    # TT_Sim_t = TT_Sim_t[indices]
    # TT_Sim_v = TT_Sim_v[indices]
    # TT_Sim_torque = TT_Sim_torque[indices]
    # TT_Sim_Distance = TT_Sim_Distance[indices]

    return TT_Sim


def inverter_losses(v_bus, v_oll, i_o_rms, power_factor, l, f_sw, u_ce0=0.8, u_d0=1, r_c=0.95e-3, r_d=0.54e-3,
                    e_ton=12e-3, e_toff=25e-3, e_d=9.5e-3):
    # e.g. [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]
    # = inverter_losses(450,230,350,0.95,75e-6,1e4,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3)
    i_ripple = (v_bus - np.sqrt(2) * v_oll) * v_oll / (2 * l * v_bus * f_sw)
    i_opk = np.sqrt(2) * i_o_rms
    m = np.sqrt(1.0 / 3.0) * v_oll * np.sqrt(2) / v_bus
    p_ct = u_ce0 * i_opk * (1 / (2 * np.pi) + m * np.cos(power_factor) / 8) + r_c * i_opk ** 2 * (
        1 / 8 + m * np.cos(power_factor) / (3 * np.pi))
    p_cd = u_d0 * i_opk * (1 / (2 * np.pi) - m * np.cos(power_factor) / 8) + r_d * i_opk ** 2 * (
        1 / 8 - m * np.cos(power_factor) / (3 * np.pi))
    i_dc = i_opk / np.pi  # DC equivalent to actual AC output current
    i_c_on = i_dc - i_ripple / 2
    i_c_off = i_dc + i_ripple / 2
    p_st = (e_ton * i_c_on + e_toff * i_c_off) * f_sw / 300.0 * v_bus / 550.0  # Test at Vce=300V, Ic=550A
    p_sd = e_d * f_sw / 300.0 * v_bus / 550.0 * i_dc
    p_total = 6.0 * (p_ct + p_cd + p_st + p_sd)
    return [p_total, p_ct, p_cd, p_st, p_sd, i_ripple]


def v_dq_pmsm(ke, poles, r_s, l_d, l_q, i_d, i_q, w_m):
    # Calculate pmsm voltages with dq model, assuming dI/dt=0
    # Using peak stator convention - magnitude of 2-phase quantity = peak value of stator phase
    # e.g. = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,Id,Iq,Rpm/30*pi)
    # ke in V/electrical rad/s
    # poles = pole pair_s
    w_e = w_m * (poles / 2)
    v_d = r_s * i_d - w_e * l_q * i_q
    v_q = r_s * i_q + w_e * (l_d * i_d + ke)
    v_s = np.sqrt(v_d ** 2 + v_q ** 2)
    #
    v_alpha = np.arctan(v_d / (v_q + 0.0001))
    i_alpha = np.arctan(i_d / (i_q + 0.0001))
    power_factor = i_alpha - v_alpha
    power_factor[np.isinf(power_factor)] = 1
    return [v_s, v_d, v_q, power_factor]


def motor_losses(stator_current, speed_rpm, resistance=6.333e-3, k_rpm_2=2.6853e-5, k_rpm_1=0.01528):
    resistive = stator_current ** 2 * resistance  # check this, 3*?
    moving = k_rpm_2 * speed_rpm ** 2 + k_rpm_1 * speed_rpm
    total = resistive + moving
    return [total, resistive, moving]


def battery_simple(energy_used, rated_energy, initial_energy, n_series):
    cell_ve = sio.loadmat('cell_V_E.mat', struct_as_record=False, squeeze_me=True)['cell_ve']
    ref_cell_rated_energy = 38  # Cell energy at 4.2V
    ref_cell_energy = cell_ve.e * ref_cell_rated_energy
    ref_cell_voltage = cell_ve.v

    n_cells = rated_energy / ref_cell_rated_energy
    # print('Battery has', str(n_cells), 'cells')
    energy_remaining = initial_energy - energy_used
    energy_remaining_cell = energy_remaining / n_cells
    v_now = n_series * np.interp(energy_remaining_cell, ref_cell_energy, ref_cell_voltage, 0, 5)

    if v_now < 2.8 * n_series:
        print('Warning: Battery overdischarged, cell V=', str(v_now))
    if v_now > 4.35 * n_series:
        print('Warning: Battery probably on fire')
    else:
        if v_now > 4.2 * n_series:
            print('Warning: Battery overcharged to', str(v_now))
    return v_now
