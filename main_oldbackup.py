import numpy as np
import matplotlib.pyplot as plt

import time
import os
import scipy.io as sio
from motorbike_functions import *  # motorbike_mech, chain_eff, motor_torque_speed, motorbike_mech2
from scipy.integrate import odeint
from scipy import integrate

# Select which corners to analyse
first_corner = 3  # 0 = first
last_corner = 7  # set to large number to ignore

# Export parameters
filename_exp = 'Python_sim.mat'
structure_exp = 'TT_sim'
# Import filenames
filename_ref_lap = 'TT_Laps_2016.mat'
filename_ref_map = 'TT_map.mat'
filename_ref_brake = 'TT_Race_2016_manual_braking_pts.mat'
structure_map = 'TT_map'
structure_lap = 'TT_Race'
var_name_brake = 'locsMin'

# Reference data mechanical specifications
N1 = 83.0
N2 = 19.0
r = 2.003 / 2 / np.pi
wMotor_ref = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
TMotor_ref = np.array([106, 106, 65.791, 45.949])

# Model bike mechanical specifications
TT_Sim = {'N': ([83.0, 19.0]), 'Km': 257 / 977,
          'constants': {'cd': 0.39, 'area': 1, 'rho': 1.204, 'm': 279.0 + 93, 'p_tyre': 1.9, 'r': 2.003 / 2 / np.pi},
          'J': {'wheel': 1.35, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 200.0, 'k_wt': 1.15}}

# Model Brake characteristic
# Let'go' time of braking instance
# Regenerative torque limit
# Braking co-efficient, relating braking torque to w

# Model motor torque vs speed characteristic
MotPMax = 2 * 50500
MotTMax = 2 * 106
TT_Sim['motor'] = {}
[TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = motor_torque_speed(MotTMax, 10500 / 30 * np.pi,
                                                                                        MotPMax, 10600 / 30 * np.pi, 20)

# NOW COMMENCE ACTUAL SIMULATION
timer = time.time()

# Preallocate model bike data structure and input constants
TT_Sim['t'] = np.array([])
TT_Sim['v'] = np.array([])
TT_Sim['Iq'] = np.array([])
TT_Sim['torque'] = np.array([])
TT_Sim['Distance'] = np.array([])
TT_Sim['J']['linear'] = TT_Sim['constants']['m'] * TT_Sim['constants']['r'] ** 2
TT_Sim['J']['r'] = 2 * TT_Sim['J']['wheel'] + TT_Sim['J']['linear'] + np.square(TT_Sim['N'][0] / TT_Sim['N'][1]) * \
                                                                      TT_Sim['J']['motor']  # J referred to wheel

# Bit of a kludge to convert dict type to mat_structure
sio.savemat('temp', {'TT_Sim': TT_Sim}, oned_as='column')
s = sio.loadmat('temp', struct_as_record=False, squeeze_me=True)
TT_Sim = s['TT_Sim']
# os.remove('temp.mat')
del s

# Import course map and reference lap data
mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
TT_map = mat_contents[structure_map]
mat_contents = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)
Ref_Race = mat_contents[structure_lap]
corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
locsmin = corners[var_name_brake]
v = r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
del r
del N2
del N1  # Just to ensure they are not used accidentally

for c in range(first_corner, min(locsmin.size - 1, last_corner)):
    corner_index = np.arange(locsmin[c], 1 + locsmin[c + 1])

    v0 = v[corner_index[0]]
    t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]], corner_index.size)

    gradientt = t
    gradient = np.interp(Ref_Race.Distance[corner_index], TT_map.dist, TT_map.gradient, 0,
                         0)  # Initially assume gradient(t) = same as lap data

    for a in range(0, 3):
        V = odeint(motorbike_mech, v0, t,
                   args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                         TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre,
                         TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t, TT_Sim.N[1], TT_Sim.N[0], gradient, gradientt))
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
    # fig2 = plt.figure(2)
    # ax = fig2.add_subplot(1,1,1)
    # ax.plot(D, V, TT_Race.Distance, v, TT_map.dist, TT_map.gradient * 100)
    # plt.plot(time, v, time[corner_index], v[corner_index], gradientt, gradient)
    # plt.xlabel('Distance')
    # plt.ylabel('v')
    # plt.xlim(D[0], D[-1])
    # fig2.show()

    MOTORSPEED = V / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * TT_Sim.constants.r)  # in rpm
    MOTORTORQUE = np.interp(MOTORSPEED, TT_Sim.motor.w * 30 / np.pi, TT_Sim.motor.t)

    T = t
    t = np.linspace(Ref_Race.t[corner_index[-1]], Ref_Race.t[corner_index[0]],
                    corner_index.size)  # as previous but flipped
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
    graidient = np.flipud(gradient)

    for a in range(0, 3):
        e_chain = 1
        V2 = odeint(motorbike_mech2, v0, t,
                    args=(TT_Sim.constants.r, TT_Sim.constants.rho, TT_Sim.constants.cd, TT_Sim.J.r,
                          TT_Sim.constants.area, TT_Sim.constants.m, TT_Sim.constants.p_tyre, TBrake, TBrake_t,
                          TT_Sim.N[1], TT_Sim.N[0], e_chain, gradient, gradientt))
        V2 = np.squeeze(V2)
        D2 = np.squeeze(Ref_Race.Distance[corner_index[-1]] + integrate.cumtrapz(V2, t, initial=0))
        gradient = np.interp(D2, TT_map.dist, TT_map.gradient, 0, 0)  # '-5.0, -5.0)

        # graidient goes out of range at end of track!!! should make circular array

    T2 = np.flipud(t)
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
        # ax.plot(T, V, '-o', TT_Race.t, v, T2, V2, '-o', gradientt, gradient * 100)
        ax.plot(Ref_Race.Distance, v, D, V, D, V2i, D2, V2, 'o')
        plt.xlim(D[0], D[-1])
        plt.ylim(0, V.max() * 1.2)
        # plt.xlim(TT_Race.t[corner_index[0]], TT_Race.t[corner_index[-1]])
        ax.plot(D[dout], V[dout], 'ro')
        fig4.show()
        plt.show()
        dout = D[-1]  # BIT BAD THIS - causes jump in Vel
    else:
        dout = np.median([D[dout]])

    # print('Braking point = %d' % dout, 'm')
    D2i = np.squeeze((D2 > dout) & (D2 < 60725))
    Di = np.squeeze(D < dout)

    #  Vboth = np.hstack((V[:dout], V2i[dout:]))

    Vboth = np.hstack((V[Di], V2[D2i]))
    Dboth = np.hstack((D[Di], D2[D2i]))

    # vals = V2i != -1
    # Dboth = D[vals]
    # Vboth = Vboth[vals]

    dt = T[1] - T[0]
    Tboth = T[0] + np.arange(0, Vboth.size * dt, dt)

    ToBoth = np.hstack((MOTORTORQUE[Di], TBrake[D2i]))

    # fig5 = plt.figure(5)
    # ax = fig5.add_subplot(1, 1, 1)
    # ax.plot(TT_Race.Distance, v, Dboth, Vboth)
    # plt.xlim(Dboth[0], Dboth[-1])
    # fig5.show()

    # PLOT THIS ONE
    # fig6 = plt.figure(6)
    # ax = fig6.add_subplot(1, 1, 1)
    # ax.plot(time, v, Tboth, ToBoth, Tboth, Vboth)
    # plt.xlim(Tboth[0], Tboth[-1])
    # fig6.show()

    tgained = Ref_Race.t[corner_index[-1]] - Tboth[-1]
    print('Time gained = %.1f s' % tgained, ' on corner %d' % c)

    if TT_Sim.t.size == 0:  # if not [] == true
        TT_Sim.t = np.hstack((TT_Sim.t, Tboth))
    else:
        TT_Sim.t = np.hstack((TT_Sim.t, Tboth - Tboth[0] + TT_Sim.t[-1]))
    TT_Sim.v = np.hstack((TT_Sim.v, Vboth))
    TT_Sim.torque = np.hstack((TT_Sim.torque, ToBoth))
    TT_Sim.Distance = np.hstack((TT_Sim.Distance, Dboth))

# Limit to one lap (and remove some crazy points)
# indices = TT_Sim_Distance < 60725

# TT_Sim_t = TT_Sim_t[indices]
# TT_Sim_v = TT_Sim_v[indices]
# TT_Sim_torque = TT_Sim_torque[indices]
# TT_Sim_Distance = TT_Sim_Distance[indices]

TT_Sim.t, indices = np.unique(TT_Sim.t, return_index=True)  # Remove stops in time
TT_Sim.v = TT_Sim.v[indices]
TT_Sim.torque = TT_Sim.torque[indices]
TT_Sim.Distance = TT_Sim.Distance[indices]

# Limit braking torque to rear wheel regenerative
torque = braking_regen(TT_Sim.v / TT_Sim.constants.r, TT_Sim.torque * TT_Sim.N[0] / TT_Sim.N[1],
                       TT_Sim.brake.LimitTorque, TT_Sim.brake.k_wt)  # limits for traction motor braking
TT_Sim.Iq = torque / TT_Sim.N[0] * TT_Sim.N[1] / TT_Sim.Km
# END OF SIMULATION
print('Simulation duration = ', str(time.time() - timer), ' seconds')

fig7 = plt.figure(7)
ax = fig7.add_subplot(2, 2, 1)
ax.plot(Ref_Race.Distance, v, TT_Sim.Distance, TT_Sim.v)
# plt.xlim(Dboth[0], Dboth[-1])

ax = fig7.add_subplot(2, 2, 2)
ax.plot(Ref_Race.t, v, TT_Sim.t, TT_Sim.v)
# plt.xlim(Dboth[0], Dboth[-1])

ax = fig7.add_subplot(2, 2, 3)
ax.plot(Ref_Race.Distance, Ref_Race.Iq, TT_Sim.Distance, TT_Sim.Iq)
# plt.xlim(Dboth[0], Dboth[-1])

ax = fig7.add_subplot(2, 2, 4)
ax.plot(Ref_Race.t, Ref_Race.Iq, TT_Sim.t, TT_Sim.Iq)
# plt.xlim(Dboth[0], Dboth[-1])
fig7.show()

plt.show()

sio.savemat(filename_exp, {structure_exp: TT_Sim}, oned_as='column')
aaa = 111
