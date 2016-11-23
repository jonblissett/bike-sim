import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from os import remove
from scipy.integrate import odeint
from scipy import integrate
from scipy.integrate import ode
# from numba import jit

#@jit
def motorbike_mech(t, v, r, rho, cd, jr, area, m, p_tyre, rpm, t_motor, n2, n1, gradient, gradient_t):
    motor_rpm = v / (1 / 60 * n2 / n1 * 2 * np.pi * r)
    motor_torque = np.interp(motor_rpm, rpm, t_motor, 0.001, 0.001)
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
    # if t > 1204.9:
    # print('t', str(t), 'v', str(v), 'torque', str(torque), 't_air', str(torque_air), 'gradient',
    #      str(np.interp(t, gradient_t, gradient, -1, -1)), 't_roll', str(torque_roll), 'e', str(e_chain), 'dvdt', str(r / jr * (torque - torque_air - torque_gradient - torque_roll)))
    # print('t', str(t), 'v', str(v), 'r', str(r))
    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt

def motorbike_mech3(t, v, r, rho, cd, jr, area, m, p_tyre, t_motor_t, t_motor, n2, n1, gradient, gradient_t):
    motor_rpm = v / (1 / 60 * n2 / n1 * 2 * np.pi * r)
    # print(motor_rpm)
    #if time - t[0] < 1:
    #    tramp = time - t[0]
    #else:
    #    tramp = 1
    #print(r)

    motor_torque = np.interp(t, t_motor_t, t_motor, -1, -1)
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
    # if t > 1204.9:
    # print('t', str(t), 'v', str(v), 'torque', str(torque), 't_air', str(torque_air), 'gradient',
    #      str(np.interp(t, gradient_t, gradient, -1, -1)), 't_roll', str(torque_roll), 'e', str(e_chain), 'dvdt', str(r / jr * (torque - torque_air - torque_gradient - torque_roll)))
    # print('t', str(t), 'v', str(v), 'r', str(r))

    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt

def motorbike_mech4(t, v, r, rho, cd, jr, area, m, p_tyre, motor_torque, n2, n1, gradient):
    motor_rpm = v / (1 / 60 * n2 / n1 * 2 * np.pi * r)

    e_chain = chain_eff(n2, n1, 12.7, 1.21 / 1.27, 1, motor_rpm / 30 * np.pi, motor_torque, 0.11, 0.00508 * 1.2)[
        0]  # Bush diameter guessed
    torque = e_chain * n1 / n2 * motor_torque  # v=rpm*1/60*n2/n1*2*np.pi*r
    g = 9.81
    # Losses *add bearing and transmission losses*
    torque_air = r / 2.0 * np.square(v) * rho * cd * area
    if v < (165.0 / 3.6):  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = r * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v * 3.6))
    else:
        torque_roll = r * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v * 3.6))

    torque_gradient = r * m * g * gradient  # -1 so bike falls off cliff
    # dydt = (T(t)-Tr-Ta)/Jr
    # print('t', str(t), 'v', str(v), 'torque', str(torque), 't_air', str(torque_air), 'gradient',
    #      str(np.interp(t, gradient_t, gradient, -1, -1)), 't_roll', str(torque_roll), 'e', str(e_chain), 'dvdt', str(r / jr * (torque - torque_air - torque_gradient - torque_roll)))
    # print('t', str(t), 'v', str(v), 'r', str(r))
    return r / jr * (torque - torque_air - torque_gradient - torque_roll)  # w=v/r; so dv/dt=r*dw/dt

#@jit
def motorbike_mech2(t, v, r, rho, cd, jr, area, m, p_tyre, t_mot, t_mott, n2, n1, e_chain, gradient, gradient_t):
    torque = e_chain * n1 / n2 * np.interp(t, t_mott, t_mot, 0, 0)  # v=rpm*1/60*n2/n1*2*np.pi*r
    g = 9.81
    # Losses *add bearing and transmission losses*
    torque_air = r / 2.0 * np.square(v) * rho * cd * area
    if v < (165.0 / 3.6):  # If vel < 165 kph from J. Bradley, 1996
        torque_roll = r * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v * 3.6))
    else:
        # print(t, v)
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
                structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings):
    # Import course map and reference lap data
    mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
    TT_map = mat_contents[structure_map]
    corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
    locsmin = corners[var_name_brake]  # -1 as matlab indexing starts at 1

    # Preallocate model bike data structure
    TT_Sim['t'] = np.array([])
    TT_Sim['v'] = np.array([])
    TT_Sim['Iq'] = np.array([])
    TT_Sim['torque'] = np.array([])
    TT_Sim['Rpm'] = np.array([])
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
    TT_Sim['Energy'] = {}
    TT_Sim['Energy']['Drive'] = np.array([])
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
        corner_index = np.arange(locsmin[c], locsmin[c + 1] - 1)
        # corner_index_end = locsmin[c + 1]

        if v_flag == -1:  # Set start vel to corner vel iff bike made that speed on previous
            v0 = v[corner_index[0]]
        else:
            v0 = v_flag

        dt_a = 0.05  #  0.008
        # t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]], corner_index.size)
        # [t, dt_a] = np.linspace(Ref_Race.t[corner_index[1]], Ref_Race.t[corner_index[0]], corner_index.size, retstep=True)
        # print('dtA = ', str(dt_a))
        t = Ref_Race.t[corner_index[0]] + np.linspace(0, corner_index.size * dt_a - dt_a, corner_index.size)
        # print('dtA = ', str(t[1] - t[0]))
        gradientt = t
        gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
                             0)  # Initially assume gradient(t) = same as lap data

        wfw = np.linspace(0, TT_Sim.motor.W_max, 1000)
        v_s = v_dq_pmsm(TT_Sim.motor.Ke, TT_Sim.motor.poles, TT_Sim.motor.Rs, TT_Sim.motor.Ld, TT_Sim.motor.Lq, 0,
                        TT_Sim.motor.T_max / TT_Sim.constants.Km, wfw)[0]

        wfwlimited = wfw[v_s > 0.9 * TT_Sim.Vdc_sim[c]]
        if wfwlimited.size == 0:
            p_weak = TT_Sim.motor.P_max
        else:
            p_weak = TT_Sim.motor.T_max * wfwlimited[0]

        if p_weak >= TT_Sim.motor.P_max:
            p_weak = TT_Sim.motor.P_max
        else:
            print('Power limited to', str(p_weak), 'by field weakening')

				# DOES THIS REALLY NEED TO BE CALCULATED EVERY TIME?
        TT_Sim.motor.w, TT_Sim.motor.t, TT_Sim.motor.p = motor_torque_speed(TT_Sim.motor.T_max, TT_Sim.motor.W_max,
                                                                            p_weak, TT_Sim.motor.W_lim, 50, 0)

        for a in range(0, 2):
            # print('a=', a)
            V = []  # list to hold solutions
            V.append(v0)  # Put y0 into solution
            # print(str(a))
            solver = ode(motorbike_mech)
            # solver.set_integrator('dopri5')     #   5th order runge-kutta
            solver.set_integrator('lsoda', with_jacobian=False)
            solver.set_initial_value(v0, t[0])
            solver.set_f_params(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                                TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre,
                                TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t, TT_Sim.N[1], TT_Sim.N[0], gradient,
                                gradientt)
            for time in t[1:]:

                V.append(solver.integrate(time))
                # position.append = ...
                # todo put field weakening calculation in here? include IR drop at least. Could just calc. motor curve for two adjecent points
                if not solver.successful():
                    print('Warning: integration not successful')

            # V = odeint(motorbike_mech, v0, t,
            #           args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
            #                 TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre,
            #                 TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t, TT_Sim.N[1], TT_Sim.N[0], gradient,
            #                 gradientt))#, hmin=1e-4,mxstep=5000
            V = np.squeeze(V)
            D = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(V, t, initial=0))
            gradient_save = gradient
            gradient = np.interp(D, TT_map.dist, TT_map.gradient, 0.0, 0.0)
            # print(rep)
            # print(D.shape)
            # print(V.shape)
            # print(TT_Race.Distance.shape)
            # print(v.shape)
            # print(TT_map.dist.shape)
            # print(TT_map.gradient.shape)
        # fig2 = plt.figure(2)
        # ax = fig2.add_subplot(1,1,1)
        # ax.plot(D, V, Ref_Race.Distance, v, TT_map.dist, TT_map.gradient * 100, D, gradient_save*100)
        # plt.xlabel('Distance')
        ##ax.plot(t, V, Ref_Race.t, v)
        # plt.ylabel('v')
        # plt.xlim(D[0], D[-1])
        ##plt.xlim(t[0], t[-1])
        # fig2.show()
        # plt.show()

        MOTORSPEED = V / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * TT_Sim.constants.r)  # in rpm
        MOTORTORQUE = np.interp(MOTORSPEED, TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t)

        T = t
        # t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[0]],
        #                 corner_index.size)  # as previous but flipped

        # print('dtB = ', str(t[1] - t[0]))
        dt_b = dt_a
        t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[-1]] - (corner_index.size * dt_b) + dt_b,
                        corner_index.size)
        # print('dtB = ', str(t[1] - t[0]))

        TBrake_t = t
        TBrake = -TT_Sim.brake.PeakTorque * TT_Sim.N[1] / TT_Sim.N[0] * np.ones(t.shape)
        rampIndex = TBrake_t > (TBrake_t[0] - TT_Sim.brake.RampTime)
        TBrake[rampIndex] = np.linspace(0, TBrake[-1], sum(rampIndex))
        TBrake_t = np.flipud(TBrake_t)
        TBrake = np.flipud(TBrake)

        # plt.close()
        # plt.plot(Tbraket, Tbrake)
        # plt.show()

        v0 = v[corner_index[-1]]

        gradientt = t
        gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
                             0)  # Initially gradient(t) = same as lap data
        gradientt = np.flipud(gradientt)

        for a in range(0, 2):
            # print('a=', a)
            e_chain = 1
            # V2 = odeint(motorbike_mech2, v0, t,
            #            args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
            #                  TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
            #                  TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt))

            V2 = []  # list to hold solutions
            V2.append(v0)  # Put y0 into solution
            # print(str(a))
            solver = ode(motorbike_mech2)
            # solver.set_integrator('dopri5')     #   5th order runge-kutta
            solver.set_integrator('lsoda', with_jacobian=False)
            solver.set_initial_value(v0, t[0])
            solver.set_f_params(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                                TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
                                TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt)

            v_max = max(V) * 1.1  # Need to stop simulating a bit after intersection, for gradient delta V
            # print('v_max = ', v_max, 't0 = ', t[0])
            for time in t[1:]:
                if V2[-1] >= v_max:
                    V2.append(v_max)  # would be neater to break the loop here, but messes up the array lengths
                else:
                    V2.append(solver.integrate(time))
                    if not solver.successful():
                        print('Warning: integration not successful')
                # position.append = ...
                # print('V2 = ', V2[j], 't = ', t[j])


            V2 = np.squeeze(V2)
            D2 = np.squeeze(Ref_Race.Distance[corner_index[-1]] + integrate.cumtrapz(V2, t, initial=0))
            # gradient = np.interp(D2, TT_map.dist, TT_map.gradient, 0, 0)  # '-5.0, -5.0)
            # gradientt = gradientt[:j]

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
            print('############################################################')
            print('Bike too slow on corner %d' % c, ', perhaps you used an Agni motor?')
            print('Higher motor torque required to achieve desired corner speed')
            if enable_warnings:
                plt.close()
                fig4 = plt.figure(4)
                ax = fig4.add_subplot(1, 1, 1)
                # ax.plot(Ref_Race.t, v, T, V, '-o', gradientt, gradient * 100)#, Ref_Race.t, 100*np.interp(Ref_Race.Distance,TT_map.dist,TT_map.gradient))#, np.flipud(t), V2, '-o'
                ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o', TT_map.dist,
                        TT_map.gradient * 100)  # , D2, gradient * 100)
                plt.xlim(D[0], D[-1])
                plt.ylim(0, V.max() * 1.2)
                # plt.xlim(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]])
                ax.plot(D[dout], V[dout], 'ro')
                fig4.show()
                plt.show()
            dout = Ref_Race.Distance[corner_index[-1]]  # BIT BAD THIS - causes jump in Vel
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
        Tboth = T[0] + np.arange(0, Vboth.size * dt - dt, dt)

        # fig5 = plt.figure(5)
        # ax = fig5.add_subplot(1, 1, 1)
        # ax.plot(TT_Race.Distance, v, Dboth, Vboth)
        # plt.xlim(Dboth[0], Dboth[-1])
        # fig5.show()

        # PLOT THIS ONE
        if enable_warnings:
            fig6 = plt.figure(6)
            ax = fig6.add_subplot(1, 1, 1)
            ax.plot(Ref_Race.t, v, Tboth, Vboth, TT_Sim.t, TT_Sim.v)  # , Tboth, ToBoth
            plt.xlim(Tboth[0], Tboth[-1])
            fig6.show()
            plt.show()

        tgained = Ref_Race.t[corner_index[-1]] - Tboth[-1]
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
        TT_Sim.P.MotorLosses = total_loss  # np.hstack((TT_Sim.P.MotorLosses, moving))
        TT_Sim.P.motRLoss = resistive  # np.hstack((TT_Sim.P.motRLoss, resistive))
        TT_Sim.P.motwLoss = moving  # np.hstack((TT_Sim.P.motwLoss, total_loss))
        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(TT_Sim.motor.Ke, TT_Sim.motor.poles, TT_Sim.motor.Rs, TT_Sim.motor.Ld,
                                                  TT_Sim.motor.Lq, 0, TT_Sim.Iq, TT_Sim.v / TT_Sim.constants.r *
                                                  TT_Sim.N[0] / TT_Sim.N[1])
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
        TT_Sim.Energy.Drive = np.trapz(TT_Sim.P.Drive, TT_Sim.t)
        TT_Sim.Vdc_sim[c + 1] = battery_simple(TT_Sim.Energy.Drive / 3600, rated_energy, initial_energy, n_series)
        print('Energy=', str(TT_Sim.Energy.Drive / 3600), 'Wh')
        print('Voltage=', str(TT_Sim.Vdc_sim[c + 1]), 'V')

    # Limit to one lap (and remove some crazy points)
    # indices = TT_Sim_Distance < 60725

    TT_Sim.Rpm = TT_Sim.v / TT_Sim.constants.r / np.pi * 30 * TT_Sim.N[0] / TT_Sim.N[1]
    # TT_Sim_t = TT_Sim_t[indices]
    # TT_Sim_v = TT_Sim_v[indices]
    # TT_Sim_torque = TT_Sim_torque[indices]
    # TT_Sim_Distance = TT_Sim_Distance[indices]

    return TT_Sim

def lap_analyse2(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings):
    ramp_start = 0.3;
    ramp_time = 1.5;    # Throttle ramp settings
    # Import course map and reference lap data
    mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
    Course_map = mat_contents[structure_map]
    corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
    locsmin = corners[var_name_brake]  # -1 as matlab indexing starts at 1

    # Preallocate model bike data structure
    TT_Sim['t'] = np.array([])
    TT_Sim['v'] = np.array([])
    TT_Sim['Iq'] = np.array([])
    TT_Sim['torque'] = np.array([])
    TT_Sim['Rpm'] = np.array([])
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
    TT_Sim['Energy'] = {}
    TT_Sim['Energy']['Drive'] = np.array([])
    TT_Sim['Vd'] = np.array([])
    TT_Sim['Vq'] = np.array([])
    TT_Sim['Vs'] = np.array([])
    TT_Sim['lean'] = np.array([])
    TT_Sim['Temperature'] = {}
    TT_Sim['Temperature']['Battery'] = np.array([])

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
        corner_index = np.arange(locsmin[c], locsmin[c + 1] - 1)
        # corner_index_end = locsmin[c + 1]

        if v_flag == -1:  # Set start vel to corner vel iff bike made that speed on previous
            v0 = v[corner_index[0]]
        else:
            v0 = v_flag
        d0 = Ref_Race.Distance[corner_index[0]]

        dt_a = 0.05  #  0.008
        # t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]], corner_index.size)
        # [t, dt_a] = np.linspace(Ref_Race.t[corner_index[1]], Ref_Race.t[corner_index[0]], corner_index.size, retstep=True)
        t = Ref_Race.t[corner_index[0]] + np.linspace(0, corner_index.size * dt_a - dt_a, corner_index.size)
        # print('dtA = ', str(t[0]), 'to ', str(t[-1]), 'step = ', str(t[1] - t[0]))
        #gradientt = t
        #gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
        #                     0)  # Initially assume gradient(t) = same as lap data

        wfw = np.linspace(0, TT_Sim.motor.W_max, 1000)
        v_s = v_dq_pmsm(TT_Sim.motor.Ke, TT_Sim.motor.poles, TT_Sim.motor.Rs, TT_Sim.motor.Ld, TT_Sim.motor.Lq, 0,
                        TT_Sim.motor.T_max / TT_Sim.constants.Km, wfw)[0]

        wfwlimited = wfw[v_s > 0.9 * TT_Sim.Vdc_sim[c]]
        if wfwlimited.size == 0:
            p_weak = TT_Sim.motor.P_max
        else:
            p_weak = TT_Sim.motor.T_max * wfwlimited[0]

        if p_weak >= TT_Sim.motor.P_max:
            p_weak = TT_Sim.motor.P_max
        else:
            print('Power limited to', str(p_weak), 'by field weakening')

				# DOES THIS REALLY NEED TO BE CALCULATED EVERY TIME?
        TT_Sim.motor.w, TT_Sim.motor.t, TT_Sim.motor.p = motor_torque_speed(TT_Sim.motor.T_max, TT_Sim.motor.W_max,
                                                                            p_weak, TT_Sim.motor.W_lim, 50, 0)
        # TT_Sim.motor.w = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
        # TT_Sim.motor.t = np.array([106, 106, 65.791, 45.949])
        # TT_Sim.motor.p = TT_Sim.motor.w * TT_Sim.motor.t

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
        R.append(TT_Sim.constants.r - 0.12 * (1 - np.cos(A[-1])))
        H.append(np.interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
        T_motor.append(ramp_start*np.interp(V[-1] / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * R[-1]),
                                 TT_Sim.motor.w / np.pi * 30, TT_Sim.motor.t, 0.001, 0.001))
        J_l.append(TT_Sim.constants.m * R[-1] ** 2)
        J_r.append(2 * TT_Sim.J.wheel + J_l[-1] +
                   np.square(TT_Sim.N[0] / TT_Sim.N[1]) * TT_Sim.J.motor)  # J referred to wheel

        solver = ode(motorbike_mech4)
        solver.set_integrator('lsoda', with_jacobian=False)
        solver.set_initial_value(V[0], t[0])
        solver.set_f_params(R[-1], TT_Sim.constants.rho, TT_Sim.constants.cd, J_r[-1], TT_Sim.constants.area,
                            TT_Sim.constants.m, TT_Sim.constants.p_tyre, T_motor[-1], TT_Sim.N[1], TT_Sim.N[0], G[-1])

        # print(t[0], V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
        # corner_index2 = np.arange(locsmin[c], locsmin[c + 3] - 1)
        wheelie = False

        for time in t[1:]:
            V.append(solver.integrate(time))
            G.append(np.interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
            A.append(np.interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
            R.append(TT_Sim.constants.r - 0.12 * (1 - np.cos(A[-1])))
            #  R[-1]=TT_Sim.constants.r
            T_motor.append(np.interp(V[-1] / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * R[-1]),
                                     TT_Sim.motor.w / np.pi * 30, TT_Sim.motor.t, 0.001, 0.001))
            if (time-t[0] < ramp_time) & (c != 0):   # ramp of 1.5 second
                tramp = ramp_start+(1-ramp_start)*(time-t[0])/ramp_time
                T_motor[-1] *= tramp
                #  print(time, time-t[0], tramp)

            T_max = TT_Sim.N[1] / TT_Sim.N[0] * R[-1] * (TT_Sim.constants.m * 9.81 * TT_Sim.constants.b /
                                                         TT_Sim.constants.h - np.square(V[-1]) * TT_Sim.constants.rho *
                                                         TT_Sim.constants.cd * TT_Sim.constants.area / 2.0)
            # print('T_max', T_max, T_motor[-1])
            if (T_motor[-1] > T_max):
                T_motor[-1] = T_max
                wheelie = True

            #  T_motor[-1] = Ref_Race.constants.Km * np.interp(time, Ref_Race.t, Ref_Race.Iq, -1, -1)

            J_l.append(TT_Sim.constants.m * R[-1] ** 2)
            J_r.append(2 * TT_Sim.J.wheel + J_l[-1] +
                       np.square(TT_Sim.N[0] / TT_Sim.N[1]) * TT_Sim.J.motor)  # J ref to wheel
            solver.set_f_params(R[-1], TT_Sim.constants.rho, TT_Sim.constants.cd, J_r[-1], TT_Sim.constants.area,
                                TT_Sim.constants.m, TT_Sim.constants.p_tyre, T_motor[-1], TT_Sim.N[1], TT_Sim.N[0],
                                G[-1])
            dD = np.squeeze(V[-1] + V[-2]) * dt_a / 2.0
            D.append(D[-1] + dD)

            # dx = np.interp(D[-1], Course_map.dist, Course_map.x) - np.interp(D[-2], Course_map.dist, Course_map.x)
            # dy = np.interp(D[-1], Course_map.dist, Course_map.y) - np.interp(D[-2], Course_map.dist, Course_map.y)
            # H.append(np.arctan2(dx,dy))
            H.append(np.interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
            dH = H[-1] - H[-2]
            if dH > np.pi / 2: dH -= np.pi
            if dH < -np.pi / 2: dH += np.pi
            w = dH / dt_a
            a_lateral = V[-1] * w
            lean.append(np.arctan(a_lateral / 9.81))

            # print(D[-1], H[-1], V[-1], lean[-1]*180/np.pi, A[-1]*180/np.pi)
            # print(time, V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
            # todo put field weakening calculation in here? include IR drop at least. Could just calc. motor curve for two adjecent points
            if not solver.successful():
                print('Warning: integration not successful')

        if (wheelie == True):
            print('Wheelie alert!')
        #V = np.squeeze(V)
        #D = np.squeeze(Ref_Race.Distance[corner_index[0]] + integrate.cumtrapz(V, t, initial=0))
        V = np.squeeze(V)
        R = np.squeeze(R)
        T_motor = np.squeeze(T_motor)
        A = np.squeeze(A)
        lean = np.squeeze(lean)
        D = np.squeeze(D)
        #gradient_save = gradient
        #gradient = np.interp(D, TT_map.dist, TT_map.gradient, 0.0, 0.0)

        MOTORSPEED = V / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * R)  # in rpm
        #MOTORTORQUE = np.interp(MOTORSPEED, TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t)
        MOTORTORQUE = T_motor

        T = t
        # t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[0]],
        #                 corner_index.size)  # as previous but flipped

        # print('dtB = ', str(t[1] - t[0]))
        dt_b = 0.05  # dt_a
        t = Ref_Race.t[corner_index[-1]] + np.linspace(0, dt_b * (1 - corner_index.size), corner_index.size)
        # print('dtB = ', str(t[1] - t[0]))

        TBrake_t = t
        TBrake = -TT_Sim.brake.PeakTorque * TT_Sim.N[1] / TT_Sim.N[0] * np.ones(t.shape)
        rampIndex = TBrake_t > (TBrake_t[0] - TT_Sim.brake.RampTime)
        TBrake[rampIndex] = np.linspace(0, TBrake[-1], sum(rampIndex))
        TBrake_t = np.flipud(TBrake_t)
        TBrake = np.flipud(TBrake)

        # plt.close()
        # plt.plot(Tbraket, Tbrake)
        # plt.show()

        v0 = v[corner_index[-1]]
        d0 = Ref_Race.Distance[corner_index[-1]]

        gradientt = t
        gradient = np.interp(Ref_Race.Distance[corner_index], Course_map.dist, Course_map.gradient, 0,
                             0)  # Initially gradient(t) = same as lap data
        gradientt = np.flipud(gradientt)

        for a in range(0, 2):
            # print('a=', a)
            e_chain = 1
            # V2 = odeint(motorbike_mech2, v0, t,
            #            args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
            #                  TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
            #                  TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt))

            V2 = []  # list to hold solutions
            D2 = []
            H2 = []
            lean2 = []
            V2.append(v0)  # Put y0 into solution
            D2.append(d0)
            H2.append(np.interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
            lean2.append(0.0)
            # print(str(a))
            solver = ode(motorbike_mech2)
            # solver.set_integrator('dopri5')     #   5th order runge-kutta
            solver.set_integrator('lsoda', with_jacobian=False)
            solver.set_initial_value(v0, t[0])
            solver.set_f_params(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                                TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
                                TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt)

            v_max = max(V) * 1.1  # Need to stop simulating a bit after intersection, for gradient delta V
            # print('v_max = ', v_max, 't0 = ', t[0])
            for time in t[1:]:
                if V2[-1] >= v_max:
                    V2.append(v_max)  # would be neater to break the loop here, but messes up the array lengths
                else:
                    V2.append(solver.integrate(time))
                    if not solver.successful():
                        print('Warning: integration not successful')

                D2.append(D2[-1] + np.squeeze(V2[-1] + V2[-2]) * -dt_b / 2.0)
                H2.append(np.interp(D2[-1], Course_map.dist, Course_map.heading, 0, 0))
                dH = -1*(H2[-1] - H2[-2])
                if dH > np.pi / 2: dH -= np.pi
                if dH < -np.pi / 2: dH += np.pi
                w = dH / dt_a
                a_lateral = V2[-1] * w
                lean2.append(np.arctan(a_lateral / 9.81))
                # position.append = ...
                # print('V2 = ', V2[j], 't = ', t[j])


        V2 = np.squeeze(V2)
        D2 = np.squeeze(D2)
        #D2 = np.squeeze(Ref_Race.Distance[corner_index[-1]] + integrate.cumtrapz(V2, t, initial=0))
        lean2 = np.squeeze(lean2)
        T2 = np.flipud(t)
        D2 = np.flipud(D2)
        V2 = np.flipud(V2)
        MOTORSPEED2 = V2 / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * TT_Sim.constants.r)  # in rpm
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
            print('No. of intersections = %d' % dout.size)

        if dout.size == 0:
            print('############################################################')
            print('Bike too slow on corner %d' % c, ', perhaps you used an Agni motor?')
            print('Higher motor torque required to achieve desired corner speed')
            if enable_warnings:
                plt.close()
                fig4 = plt.figure(4)
                ax = fig4.add_subplot(1, 1, 1)
                # ax.plot(Ref_Race.t, v, T, V, '-o', gradientt, gradient * 100)#, Ref_Race.t, 100*np.interp(Ref_Race.Distance,TT_map.dist,TT_map.gradient))#, np.flipud(t), V2, '-o'
                ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o', Course_map.dist,
                        Course_map.gradient * 100)  # , D2, gradient * 100)
                plt.xlim(D[0], D[-1])
                plt.ylim(0, V.max() * 1.2)
                # plt.xlim(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]])
                ax.plot(D[dout], V[dout], 'ro')
                fig4.show()
                plt.show()
            dout = Ref_Race.Distance[corner_index[-1]]  # BIT BAD THIS - causes jump in Vel
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
        RPMBoth = np.hstack((MOTORSPEED[Di], MOTORSPEED2[D2i]))
        leanBoth = np.hstack((lean[Di], lean2[D2i]))

        # vals = V2i != -1
        # Dboth = D[vals]
        # Vboth = Vboth[vals]

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
            ax.plot(Ref_Race.t, v, Tboth, Vboth, TT_Sim.t, TT_Sim.v)  # , Tboth, ToBoth
            plt.xlim(Tboth[0], Tboth[-1])
            fig6.show()
            plt.show()

        tgained = Ref_Race.t[corner_index[-1]] - Tboth[-1]
        print('Time gained = %.1f s' % tgained, ' on corner %d' % c)

        if TT_Sim.t.size == 0:  # if not [] == true
            TT_Sim.t = np.hstack((TT_Sim.t, Tboth))
        else:
            TT_Sim.t = np.hstack((TT_Sim.t, Tboth - Tboth[0] + TT_Sim.t[-1]))
        TT_Sim.v = np.hstack((TT_Sim.v, Vboth))
        TT_Sim.Rpm = np.hstack((TT_Sim.Rpm, RPMBoth))
        TT_Sim.torque = np.hstack((TT_Sim.torque, ToBoth))
        TT_Sim.Distance = np.hstack((TT_Sim.Distance, Dboth))
        TT_Sim.lean = np.hstack((TT_Sim.lean, leanBoth))

        # Now postprocess to find energy
        TT_Sim.t, indices = np.unique(TT_Sim.t, return_index=True)  # Remove stops in time
        TT_Sim.v = TT_Sim.v[indices]
        TT_Sim.Rpm = TT_Sim.Rpm[indices]
        TT_Sim.torque = TT_Sim.torque[indices]
        TT_Sim.Distance = TT_Sim.Distance[indices]
        TT_Sim.lean = TT_Sim.lean[indices]

        # Limit braking torque to rear wheel regenerative
        torque = braking_regen(TT_Sim.Rpm / 30 * np.pi * TT_Sim.N[1] / TT_Sim.N[0], TT_Sim.torque * TT_Sim.N[0] / TT_Sim.N[1],
                               TT_Sim.brake.LimitTorque, TT_Sim.brake.k_wt)  # limits for traction motor braking
        # TT_Sim.Iq = np.hstack((TT_Sim.Iq, torque / TT_Sim.N[0] * TT_Sim.N[1] / TT_Sim.constants.Km))
        TT_Sim.Iq = torque / TT_Sim.N[0] * TT_Sim.N[1] / TT_Sim.constants.Km
        [total_loss, resistive, moving] = motor_losses(TT_Sim.Iq, TT_Sim.Rpm, 6.333e-3, 2.6853e-5, 0.01528)
        TT_Sim.P.MotorLosses = total_loss  # np.hstack((TT_Sim.P.MotorLosses, moving))
        TT_Sim.P.motRLoss = resistive  # np.hstack((TT_Sim.P.motRLoss, resistive))
        TT_Sim.P.motwLoss = moving  # np.hstack((TT_Sim.P.motwLoss, total_loss))
        [v_s, v_d, v_q, power_factor] = v_dq_pmsm(TT_Sim.motor.Ke, TT_Sim.motor.poles, TT_Sim.motor.Rs, TT_Sim.motor.Ld,
                                                  TT_Sim.motor.Lq, 0, TT_Sim.Iq, TT_Sim.Rpm / 30 * np.pi)
        TT_Sim.Vs = v_s  # np.hstack((TT_Sim.Vs, v_s))
        TT_Sim.Vd = v_d  # np.hstack((TT_Sim.Vd, v_d))
        TT_Sim.Vq = v_q  # np.hstack((TT_Sim.Vq, v_q))
        TT_Sim.PF = power_factor  # .hstack((TT_Sim.PF, power_factor))
        TT_Sim.Vdc = np.hstack((TT_Sim.Vdc, TT_Sim.Vdc_sim[c] * np.ones(TT_Sim.t.size - TT_Sim.Vdc.size)))
        p_drive_loss = inverter_losses(TT_Sim.Vdc, TT_Sim.Vs, abs(TT_Sim.Iq) / np.sqrt(2), TT_Sim.PF, 82e-6, 13e3, 0.8,
                                       1, 0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
        # TT_Sim.P.DriveLosses = np.hstack((TT_Sim.P.DriveLosses, p_drive_loss))
        TT_Sim.P.DriveLosses = p_drive_loss
        TT_Sim.P.Mech = TT_Sim.Rpm / 30 * np.pi * TT_Sim.Iq * TT_Sim.constants.Km
        TT_Sim.P.Motor = TT_Sim.P.Mech + TT_Sim.P.MotorLosses
        TT_Sim.P.Drive = TT_Sim.P.Motor + TT_Sim.P.DriveLosses
        # TT_Sim.Idc = TT_Sim.P.Drive / TT_Sim.Vdc
        TT_Sim.Energy.Drive = np.trapz(TT_Sim.P.Drive, TT_Sim.t)
        TT_Sim.Vdc_sim[c + 1] = battery_simple(TT_Sim.Energy.Drive / 3600, rated_energy, initial_energy, n_series)
        print('Energy=%dWh' % (TT_Sim.Energy.Drive / 3600))
        print('Voltage=%.1fV' % (TT_Sim.Vdc_sim[c + 1]))

    n_cells = rated_energy / 38
    I_cell = TT_Sim.P.Drive / TT_Sim.Vdc / n_cells * n_series
    TT_Sim.Temperature.Battery = battery_heat(n_cells, TT_Sim.constants.IRcell, 0.19, I_cell, TT_Sim.t)
    TT_Sim.gradient = np.interp(TT_Sim.Distance, Course_map.dist, Course_map.gradient, 0, 0)
    TT_Sim.constants.R = TT_Sim.v * 30 / np.pi / TT_Sim.Rpm * TT_Sim.N[0] / TT_Sim.N[1]  # v/w = r

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
    v_s = np.sqrt(v_d ** 2 + v_q ** 2)
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

def battery_heat(n_cells, IR_cell, m_cell, I_cell, t):
    shc = 1014
    P_cell = I_cell ** 2 * IR_cell
    E_cell = integrate.cumtrapz(P_cell, t, initial=0)
    T_cell = E_cell / shc / m_cell
    print('Battery temperature rise (end)=%.1f°K' % T_cell[-1])
    return T_cell

def energy_losses(file, name, IR):
    data_file = sio.loadmat(file, struct_as_record=False, squeeze_me=True)
    Data = data_file[name]
    rhoCdA = Data.constants.rho * Data.constants.cd * Data.constants.area
    [Data.P.air, Data.P.roll] = power_mech_losses_no_grad(Data.v, rhoCdA, Data.constants.m,
                                                                              Data.constants.p_tyre)

    Data.Idc = Data.P.Drive / Data.Vdc
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

def wheel_forces(p, h, b, r, m, force_air, torque_wheel, show_figure):
    # wheel_forces(1.415, 0.6, 0.7, r, cd, area, np.square(v) * rho * cd * area / 2.0, torque_wheel, 1)
    # p=1.415;  % Horizontal distance between wheel centers
    # h=0.6;  % Vertical distance from ground to center of mass
    # b=0.7;  % Horizontal distance from REAR wheel to center of mass
    # force_air = 1 / 2.0 * np.square(v) * rho * cd * area
    g = 9.81
    force_wheel = torque_wheel / r
    force_rear = (m * g * (p - b)+ h * (force_wheel + force_air)) / p
    force_front = (m * g * b - h * (force_wheel + force_air)) / p

    mu_tyre = 0.8; # Friction co-efficient of tyre
    torque_MaxFront = force_front * mu_tyre * r
    torque_MaxRear = force_rear * mu_tyre * r

    n = np.linspace(1,force_rear.size,force_rear.size)
    if show_figure == 1:
        fig, ax1 = plt.subplots()
        ax1.plot(n, force_rear, n, force_front)
        ax1.set_xlabel('samples')
        # Make the y-axis label and tick labels match the line color.
        ax1.set_ylabel('Force (N)', color='b')
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        ax2 = ax1.twinx()
        ax2.plot(n, torque_MaxRear, n, torque_MaxFront)
        ax2.set_ylabel('Max torque (Nm)', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        plt.title('Tyre contact forces')
        # plt.show()
        fig.show()

    return [force_rear, force_front, torque_MaxRear, torque_MaxFront]

