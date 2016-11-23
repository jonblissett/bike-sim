import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import integrate

import time
from motorbike_functions import lap_analyse2, lap_analyse, motor_torque_speed

# done: bus voltage estimation during sim,
# done: add field weakening torque limiting - make electrical_functions.py
# can call torque-speed function with new P lim = f(Vbus)
# still need correct value of id for loss calculation
# TODO postprocessing - losses etc
# TODO assess gains of front wheel brake
# TODO find optimal N1/N2 for ???
# done: use correct speed value if unreachable
# TODO check if regen chain efficiency is correct?
# TODO try torque vs speed as a ramp? as in like the power requirement

enable_warnings = False

# Find list of N1&N2 to test with

# Select which corners to analyse
first_corner = 0    # 6#62  # 0 = first
last_corner = 15#16     # 7#63  # set to large number to use all

end_dist = 4600  # Distance at lap end for timing

# Export parameters
filename_exp = 'Python_sim_Portimao_DM1.mat'
structure_exp = 'Portimao_sim_DM1'
# Import filenames
filename_ref_lap = 'Portimao_RaceDM1.mat'
filename_ref_map = 'Portimao_map.mat'
filename_ref_brake = 'Portimao_RaceDM1.mat'
structure_map = 'Portimao_map'
structure_lap = 'RaceDM1'
var_name_brake = 'locsMin'

# Reference data mechanical specifications
# N1 = 83.0
# N2 = 19.0
# r = 2.003 / 2 / np.pi
# wMotor_ref = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TMotor_ref = np.array([106, 106, 65.791, 45.949])

# Model bike mechanical specifications
#  Model Brake characteristic
#  - Release time of braking instance
#  - Regenerative torque limit
#  - Braking co-efficient, relating braking torque to w
TT_Sim = {'N': ([83.0, 18.0]),
          'constants': {'Km': 257/977, 'cd': 0.43, 'area': 1, 'rho': 1.204, 'm': 235 + 90, 'p_tyre': 1.9,
                        'r': 2.18 / 2 / np.pi},
          'J': {'wheel': 1.1, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 1.5}}
# UoN02, zolder, 108s4p 10Ah, 235kg, 320kg with JJ & no leathers


# CHECK WHEEL INERTIA - IS THIS FRONT AND REAR OR TWO FRONT WHEELS?

# Model bike electrical specifications
#  motor torque vs speed characteristic
n_series = 108
rated_energy = n_series*3.7*40
initial_energy = rated_energy
turns = 11.5
L_core = 150
TT_Sim['motor'] = {}
TT_Sim['motor']['T_max'] = 600 * TT_Sim['constants']['Km']
TT_Sim['motor']['W_max'] = 8600 / 30 * np.pi
TT_Sim['motor']['W_lim'] = 9000 / 30 * np.pi
TT_Sim['motor']['P_max'] = 350000  # 50500
TT_Sim['motor']['poles'] = 12
TT_Sim['motor']['Ld'] = 53e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Lq'] = 61e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Rs'] = 0.007313 * (turns/11.5) * L_core / 150
TT_Sim['motor']['Ke'] = 0.3467/6 * (turns/11.5) * L_core / 150  # Back-emf constant, V/erad/s
[TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = motor_torque_speed(TT_Sim['motor']['T_max'],
                                                                                        TT_Sim['motor']['W_max'],
                                                                                        TT_Sim['motor']['P_max'],
                                                                                        TT_Sim['motor']['W_lim'], 50, 1)

# TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
# TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']
# SETTING LIKE THIS DOESN@T WORK AS THEY ARE RECALCULATED

# END OF USER PARAMETERS

timer = time.time()

Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
# v = Ref_Race.vGPS
v = Ref_Race.vGPS# *0.8
# v = r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
# v = v * 0.9
# del r
# del N2
# del N1  # Just to ensure they are not used accidentally

TT_Sim = lap_analyse2(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                     structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings)

print('Simulation duration =', str(time.time() - timer), 'seconds')
sio.savemat(filename_exp, {structure_exp: TT_Sim}, oned_as='column')
print('Simulation saved to', filename_exp, 'as structure named', structure_exp)

print('Motor max speed', str(max(TT_Sim.v) / TT_Sim.constants.r * 30 / np.pi * TT_Sim.N[0] / TT_Sim.N[1]))
print('Bike max speed (mph) = ', str(max(TT_Sim.v) * 2.23))
print('Simulated lap time (s) = ', str((TT_Sim.t[TT_Sim.Distance < end_dist])[-1]))
print('Simulated lap speed (mph) = ', str(4.545*0.621/TT_Sim.t[-1]*3600))
# print(integrate.cumtrapz(TT_Sim.v, TT_Sim.t))

fig7 = plt.figure(7)
ax = fig7.add_subplot(3, 1, 1)
ax.plot(Ref_Race.Distance, v, TT_Sim.Distance, TT_Sim.v, 'o')
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
ax = fig7.add_subplot(3, 1, 3)
ax.plot(Ref_Race.Distance, Ref_Race.Iq, TT_Sim.Distance, TT_Sim.Iq)
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
ax = fig7.add_subplot(3, 1, 2)
ax.plot(Ref_Race.Distance, Ref_Race.lean*180/np.pi, TT_Sim.Distance, TT_Sim.lean*180/np.pi)
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
fig7.show()
#ax = fig7.add_subplot(2, 2, 4)
#ax.plot(Ref_Race.t, Ref_Race.Iq, TT_Sim.t, TT_Sim.Iq)
#plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
#fig7.show()
#ax = fig7.add_subplot(4, 1, 2)
#ax.plot(Ref_Race.t, v, TT_Sim.t, TT_Sim.v, 'o')
#plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])

plt.show()
