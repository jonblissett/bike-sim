#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import motorbike_functions as bike
from pylab import *

filename_ref_map = 'TT_map.mat'
structure_map = 'TT_map'

TT_map = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)[structure_map]

TT = sio.loadmat('TT_Laps_2016.mat', struct_as_record=False, squeeze_me=True)
TT['TT_Race'].gradient = np.interp(TT['TT_Race'].Distance, TT_map.dist, TT_map.gradient, 0, 0)
TT['TT_Qualy'].gradient = np.interp(TT['TT_Qualy'].Distance, TT_map.dist, TT_map.gradient, 0, 0)


Data_A = TT['TT_Race']
Data_B = TT['TT_Qualy']

TT_Sim = {'N': ([83.0, 19.0]),
          'constants': {'Km': 257 / 977, 'cd': 0.38, 'area': 1, 'rho': 1.204, 'm': 270.0 + 90, 'p_tyre': 1.9,
                        'r': 2.003 / 2 / np.pi},
          'J': {'wheel': 1.35, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 200.0, 'k_wt': 1.15}}
turns = 11.5
L_core = 150
TT_Sim['motor'] = {}
# TT_Sim['motor']['T_max'] = 106
# TT_Sim['motor']['W_max'] = 10500 / 30 * np.pi
# TT_Sim['motor']['W_lim'] = 10600 / 30 * np.pi
# TT_Sim['motor']['P_max'] = 50500
TT_Sim['motor']['poles'] = 12
TT_Sim['motor']['Ld'] = 53e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Lq'] = 61e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Rs'] = 0.007313 * (turns/11.5) * L_core / 150
TT_Sim['motor']['Ke'] = 0.3467/6 * (turns/11.5) * L_core / 150


r = TT_Sim['constants']['r']
m = TT_Sim['constants']['m']
rho = TT_Sim['constants']['rho']
cd = TT_Sim['constants']['cd']
area = TT_Sim['constants']['area']
p_tyre = TT_Sim['constants']['p_tyre']
n1 = TT_Sim['N'][0]
n2 = TT_Sim['N'][1]
Km = TT_Sim['constants']['Km']
poles = TT_Sim['motor']['poles']
Ld = TT_Sim['motor']['Ld']
Lq = TT_Sim['motor']['Lq']
Rs = TT_Sim['motor']['Rs']
Ke = TT_Sim['motor']['Ke']

t = Data_A.t
w_m = Data_A.Rpm / 30 * np.pi
w = w_m / n1 * n2
v = w_m / n1 * n2 * r
Iq = Data_A.Iq
gradient = Data_A.gradient
Vdc = Data_A.Vdc

# know motor torque from iq, k*Iq*N1/N2-torque_roll -torque_air = torque_gradient
# include wheel forces
#    e_chain = chain_eff(n2, n1, 12.7, 1.21 / 1.27, 1, motor_rpm / 30 * np.pi, motor_torque, 0.11, 0.00508 * 1.2)[0]
# wheel_torque = e_chain * n1 / n2 * motor_torque  # v=rpm*1/60*n2/n1*2*np.pi*r
g = 9.81
# Losses *add bearing and transmission losses*
p_air = w * r / 2.0 * np.square(v) * rho * cd * area
p_roll = v * 0
th = 165.0 / 3.6
p_roll[v < th] = w[v < th] * r * m * g * (0.0085 + 0.018 / p_tyre + 1.59e-06 / p_tyre * np.square(v[v < th] * 3.6))
p_roll[v > th] = w[v > th] * r * m * g * (0.018 / p_tyre + 2.91e-06 / p_tyre * np.square(v[v > th] * 3.6))
# torque_gradient = torque_motor - torque_roll - torque_air
p_gradient = w * r * m * g * gradient

# MOTOR LOSSES
# [total_loss, resistive, moving] = motor_losses(TT_Sim.Iq, TT_Sim.v / TT_Sim.constants.r / np.pi * 30 *
 #                                               TT_Sim.N[0] / TT_Sim.N[1], 6.333e-3, 2.6853e-5, 0.01528)
[p_MotorLosses, p_motRLoss, p_motwLoss] = bike.motor_losses(Iq, w_m * 30 / np.pi)

# Cable, drive losses
[v_s, PF] = bike.v_dq_pmsm(Ke, poles, Rs, Ld, Lq, 0, Iq, w_m)[::3]

p_DriveLosses = bike.inverter_losses(Vdc, v_s, abs(Iq) / np.sqrt(2), PF, 82e-6, 13e3, 0.8,
                               1, 0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
# TT_Sim.P.DriveLosses = np.hstack((TT_Sim.P.DriveLosses, p_drive_loss))
p_Mech = w_m * Km * Iq
p_Motor = p_Mech + p_MotorLosses
p_Drive = p_Motor + p_DriveLosses

IR = 0.1
p_DCLosses = (p_Drive / Vdc)**2 * IR

print(np.trapz(p_air, t),np.trapz(p_roll, t))

N = 5
menMeans = (20, 35, 30, 35, 27)
menStd = (2, 3, 4, 1, 2)

ind = np.arange(N)  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

womenMeans = (25, 32, 34, 20, 25)
womenStd = (3, 5, 2, 3, 3)
rects2 = ax.bar(ind + width, womenMeans, width, color='y', yerr=womenStd)

# add some text for labels, title and axes ticks
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(ind + width)
ax.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))

ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))


def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                '%d' % int(height),
                ha='center', va='bottom')


autolabel(rects1)
autolabel(rects2)

plt.show()

#fig7 = plt.figure(7)
#ax = fig7.add_subplot(2, 1, 1)
#ax.plot(t, v, t, Iq/10)
#plt.xlim(0, 70)
#ax = fig7.add_subplot(2, 1, 2)
#ax.plot(t, p_air, t, p_roll, t, p_MotorLosses, t, p_DriveLosses, t, p_DCLosses)
#plt.xlim(0, 70)
#plt.ylim(0, 50e3)
#plt.xlabel('Time (s)')
#plt.ylabel('Loss power (W)')
# fig7.show()

tot = p_air + p_roll + p_MotorLosses + p_DriveLosses + p_DCLosses

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(t, p_air, t, p_roll, t, p_MotorLosses, t, p_DriveLosses, t, p_DCLosses, t, tot, 'g-')
ax2.plot(t, v, 'm-')

ax1.set_xlabel('Time (s)')
ax1.set_xlim(0, 70)
ax1.set_ylim(0, 70e3)
ax1.set_ylabel('Power in system losses (W)', color='g')
ax2.set_ylabel('Velocity (m/s)', color='m')

plt.show()

ax = axes([0.1, 0.1, 0.8, 0.8])

fracs = [15, 30, 45, 10]
pie(fracs, labels=['Frogs', 'Hogs', 'Dogs', 'Logs'], autopct='%1.1f%%')
title('Test pie')
plt.show()