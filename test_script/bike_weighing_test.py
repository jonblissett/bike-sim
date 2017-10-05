from motorbike_functions import motor_torque, motor_torque_speed, motor_sizing, battery_simple, motor_currents, \
    motor_saturation_interp, motor_current, charge_battery, weigh, motor_losses, inverter_losses, vs_pmsm
import numpy as np
import scipy.io as sio
try:
    import matplotlib.pyplot as plt
except:
    pass
sim = {'N': ([83.0, 19.0]),
          'constants': {'cd': 0.32, 'area': 1, 'rho': 1.204, 'p_tyre': 1.9,
                        'r': 2.16 / 2 / np.pi, 'b': 0.7, 'h': 0.6, 'k_tyre': 0.7},
          'J': {'wheel': 1.35 - 0.445, 'motor': 0.0233},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 1.615}
          }

sim['J'] = {}
sim['battery'] = {}
sim['battery']['series'] = 120
sim['battery']['parallel'] = 3
sim['battery']['cellAh'] = 10
sim['battery']['cellVnom'] = 3.8
sim['battery']['E_density'] = 3.8 * 6 * 40 / 4.8
sim['battery']['cellIR'] = 0.00438
cell_ve = sio.loadmat('../data_import/cell_V_E.mat', struct_as_record=False, squeeze_me=True)['cell_ve']

sim['battery']['cell_ve'] = {'v': np.array(cell_ve.v), 'e': np.array(cell_ve.e)}
sim['battery']['IR'] = sim['battery']['cellIR'] * sim['battery']['series'] / sim['battery']['parallel']
sim = charge_battery(sim, 1.096)

sim['motor'] = {}
sim['motor']['N'] = 11.5
sim['motor']['L_core'] = 150.0
sim['motor']['W_max'] = 10500 / 30 * np.pi
sim['motor']['W_lim'] = 10600 / 30 * np.pi
sim['motor']['poles'] = 12
sim['motor']['manufacturer'] = 'Parker'

sim['drive'] = {}
sim['drive']['n'] = 1.0
sim['drive']['m'] = 5.0
sim['drive']['I_max'] = 900.0

v_max = battery_simple(sim, 0, 3)[0]
sim['battery']['V_max'] = battery_simple(sim, 0, 3)[0]
sim = motor_sizing(sim)
sim = weigh(sim)

krpm2 = 2.6853e-5 / 150 * sim['motor']['L_core']
krpm1 = 0.01528 / 150 * sim['motor']['L_core']

[total_loss, resistive, moving] = motor_losses(sim['motor']['I_pk'] * 0.72, 0.4 * 30 / np.pi * sim['motor']['W_max'],
                                               sim['motor']['Rs'], krpm2, krpm1)

print('Motor torque:' + str(motor_torque(sim['motor']['co'], sim['motor']['I_pk'] * 0.75)))
print('Motor speed:' + str(0.4 * 30 / np.pi * sim['motor']['W_max']))
print('Bike mass estimate: ' + str(sim['mass']['bike']))
print('Bike with rider mass: ' + str(sim['constants']['m']))

print(sim['motor'])
print(sim['mass'])
print('Krpm1 = ' + str(krpm1))
print('Krpm2 = ' + str(krpm2))
print('Losses (resistive,moving): (' + str(resistive) + ', ' + str(moving) + ')')

w = np.linspace(0, sim['motor']['W_max'])
i = np.linspace(0, sim['motor']['I_pk'])
[total_loss, resistive, moving] = motor_losses(i, 30 / np.pi * w, sim['motor']['Rs'], krpm2, krpm1)

Vs = 400  # vs_pmsm
PF = 0.9
DriveLosses_to = sim['drive']['n'] * inverter_losses(sim['battery']['cellVnom'] * sim['battery']['series'], Vs,
                                                                    i / np.sqrt(2) /
                                                                    sim['drive']['n'], PF, 82e-6,
                                                                    13e3, 0.8, 1, 0.95e-3, 0.54e-3, 12e-3,
                                                                    25e-3, 9.5e-3)[0]
Vs = vs_pmsm(sim['motor']['Ke'], sim['motor']['poles'], sim['motor']['Rs'], sim['motor']['Ld'],
             sim['motor']['Lq'], 0, sim['motor']['I_pk'], w)
DriveLosses_w = sim['drive']['n'] * inverter_losses(sim['battery']['cellVnom'] * sim['battery']['series'], Vs,
                                                    sim['motor']['I_pk'] / np.sqrt(2) /
                                                                    sim['drive']['n'], PF, 82e-6,
                                                                    13e3, 0.8, 1, 0.95e-3, 0.54e-3, 12e-3,
                                                                    25e-3, 9.5e-3)[0]
fig, ax1 = plt.subplots()
ax1.plot([],[],color='b', label='Motor core losses from speed', linewidth=5)
ax1.plot([],[],color='g', label='Motor resistive losses from torque', linewidth=5)
ax1.plot([],[],color='r', label='Drive losses at peak torque from speed', linewidth=5)
ax1.plot([],[],color='c', label='Drive losses from torque', linewidth=5)

ax1.plot(motor_torque(sim['motor']['co'], i) / sim['motor']['T_pk'], resistive, 'b-', w / w[-1], moving, 'g-',
         motor_torque(sim['motor']['co'], i) / sim['motor']['T_pk'], DriveLosses_w, 'r-',
         motor_torque(sim['motor']['co'], i) / sim['motor']['T_pk'], DriveLosses_to, 'c-')
ax1.set_xlabel('Normalised motor speed and torque')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Power loss (W)')#, color='b')
#for tl in ax1.get_yticklabels():
#    tl.set_color('b')
ax1.set_ylim(0, 7000)
#ax2 = ax1.twinx()
#ax2.plot(w / w[-1], moving, 'g-')
#ax2.set_ylabel('Power (kW)', color='r')
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')
#plt.title('Motor Torque and Power vs Speed')

ax1.legend(loc='upper left')

plt.show()
fig.show()


v = w * sim['N'][1] / sim['N'][0] * sim['constants']['r']
# Limit motor torque for wheelies
T_max = sim['N'][1] / sim['N'][0] * sim['constants']['r'] * (sim['constants']['m'] * 9.81 * sim['constants']['b'] / sim['constants']['h']
                                                   - v ** 2 * sim['constants']['rho'] *
                                                   sim['constants']['cd'] *
                                                   sim['constants']['area'] / 2.0)
P_max = T_max * w

fig, ax1 = plt.subplots()
ax1.plot([],[],color='b', label='Maximum motor torque', linewidth=5)
ax1.plot([],[],color='g', label='Maximum power', linewidth=5)

ax1.plot(v, T_max, 'b-')
ax1.set_xlabel('Bike speed (m/s)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Motor torque limit (Nm)', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax1.set_ylim(0, 300)

ax2 = ax1.twinx()
ax2.plot(v, P_max / 1e3, 'g-')
ax2.set_ylabel('Power (kW)', color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')
#plt.title('Motor Torque and Power vs Speed')

ax1.legend(loc='lower center')
ax1.set_xlim(0, 86)
ax2.set_xlim(0, 86)

plt.show()
fig.show()