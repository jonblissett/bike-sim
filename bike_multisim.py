import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
# from scipy import integrate

import time
from motorbike_functions import lap_analyse, motor_torque_speed

# 1) Set paramater to change
# 2) Find best gear ratio for new parameters - probably simulate for only limited corners
# 3) Simulate full lap
# 4) Store results to analyse in MATLAB


enable_warnings = False
enable_messages = False

# Select which corners to analyse
first_corner = 0    # 6#62  # 0 = first
last_corner = 99    # 7#63  # set to large number to use all

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
# wMotor_ref = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TMotor_ref = np.array([106, 106, 65.791, 45.949])

# Model bike mechanical specifications
#  Model Brake characteristic
#  - Release time of braking instance
#  - Regenerative torque limit
#  - Braking co-efficient, relating braking torque to w
TT_Sim = {'N': ([83.0, 19.0]),
          'constants': {'Km': 257 / 977, 'cd': 0.38, 'area': 1, 'rho': 1.204, 'm': 270.0 + 90, 'p_tyre': 1.9,
                        'r': 2.003 / 2 / np.pi},
          'J': {'wheel': 1.35, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 200.0, 'k_wt': 1.15}}

# CHECK WHEEL INERTIA - IS THIS FRONT AND REAR OR TWO FRONT WHEELS?

# Model bike electrical specifications
#  motor torque vs speed characteristic
rated_energy = 18500
initial_energy = 19000
n_series = 120
turns = 11.5
L_core = 150
TT_Sim['motor'] = {}
TT_Sim['motor']['T_max'] = 106
TT_Sim['motor']['W_max'] = 10500 / 30 * np.pi
TT_Sim['motor']['W_lim'] = 10600 / 30 * np.pi
TT_Sim['motor']['P_max'] = 50500
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
v = r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
del r
del N2
del N1  # Just to ensure they are not used accidentally

ratios = [[83, 20], [83, 19], [79, 20], [73, 20], [67, 20]]

TT_Sim_INIT = TT_Sim
Lap_Times = []
Lap_Energy = []
m = TT_Sim['constants']['m']

for rat_in in range(0, len(ratios)):    # better, while f() > better than previous ratio
    for dothisdifferently in range(0, 2):
        TT_Sim = TT_Sim_INIT  # Reset variable to initial conditions
        # TT_Sim['constants']['m'] = m
        TT_Sim['N'][0] = ratios[rat_in][0]
        TT_Sim['N'][1] = ratios[rat_in][1]
        TT_Sim = lap_analyse(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                             structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings)
        print('Simulation duration =', str(time.time() - timer), 'seconds')
        print('Motor max speed', str(max(TT_Sim.v) / TT_Sim.constants.r * 30 / np.pi * TT_Sim.N[0] / TT_Sim.N[1]))
        print("Bike max speed (mph) = ", str(max(TT_Sim.v) * 2.23))
        print('Simulated lap time (s) = ', str(TT_Sim.t[-1]))
        print('Simulated lap speed (mph) = ', str(37.733/TT_Sim.t[-1]*3600))
        print('Energy=', str(TT_Sim.Energy / 3600), 'Wh')
        # Estimate weight from required battery capacity
        energy_spare = rated_energy - TT_Sim.Energy / 3600
        m = TT_Sim.constants.m - energy_spare / 185    # 185Wh/kg
        rated_energy -= energy_spare
        initial_energy -= energy_spare
        # resize battery here so fw and losses are accurate
        print('')
        print('')
        print('New mass =', m)
        print('New batt =', rated_energy)
        print('New batt =', initial_energy)
        print('')
        print('')
        # Adjust corner speeds for new bike mass - shouldn't change, but suspension has to settle (JJ)
        # v = v - (new_mass - previous_mass)*some_constant

    f_name = "%s_%s" % (filename_exp, rat_in)
    sio.savemat(f_name, {structure_exp: TT_Sim}, oned_as='column')
    print("Simulation saved to", f_name, "as structure named", structure_exp)
    print("N={0:d}:{1:d}".format(TT_Sim.N[0], TT_Sim.N[1]))
    # Do we pick lowest energy, fastest lap, or minimise lap*energy ...? cost function
    # Wh*s ? --> minimise 'action'
    # cost function should favour energy over lap time, as lap time will decrease and energy increase with more mass
    Lap_Energy.append(TT_Sim.Energy / 3600)
    Lap_Times.append(TT_Sim.t[-1])
    m = TT_Sim.constants.m


print('')
print(Lap_Energy)
print(Lap_Times)
print([a*b/1000 for a, b in zip(Lap_Energy, Lap_Times)])
print('')

best_ratio = min(range(len(Lap_Times)), key=Lap_Times.__getitem__)

print('Best ratio = ', best_ratio)

last_corner = 99
# Simulate
TT_Sim = TT_Sim_INIT
TT_Sim['N'][0] = ratios[best_ratio][0]
TT_Sim['N'][1] = ratios[best_ratio][1]
TT_Sim = lap_analyse(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                     structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings)
print('Simulation duration =', str(time.time() - timer), 'seconds')
print('Motor max speed', str(max(TT_Sim.v) / TT_Sim.constants.r * 30 / np.pi * TT_Sim.N[0] / TT_Sim.N[1]))
print('Bike max speed (mph) = ', str(max(TT_Sim.v) * 2.23))
print('Simulated lap time (s) = ', str(TT_Sim.t[-1]))
print('Simulated lap speed (mph) = ', str(37.733 / TT_Sim.t[-1] * 3600))
print('Energy=', str(TT_Sim.Energy / 3600), 'Wh')


TT_Sim = TT_Sim_INIT
TT_Sim['N'][0] = ratios[best_ratio][0]
TT_Sim['N'][1] = ratios[best_ratio][1]

# TODO this is no good as new weights will affect gear ratio choice

# Adjust corner speeds for new bike mass
# v = v - (new_mass - previous_mass)*some_constant

# Re-simulate

# Save

f_name = "%s_%s" % (filename_exp, 'best')
sio.savemat(f_name, {structure_exp: TT_Sim}, oned_as='column')
print("Simulation saved to", f_name, "as structure named", structure_exp)
print("N={0:d}:{1:d}".format(TT_Sim.N[0], TT_Sim.N[1]))

fig7 = plt.figure(7)
ax = fig7.add_subplot(2, 2, 1)
ax.plot(Ref_Race.Distance, v, TT_Sim.Distance, TT_Sim.v, 'o')
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
ax = fig7.add_subplot(2, 2, 2)
ax.plot(Ref_Race.t, v, TT_Sim.t, TT_Sim.v, 'o')
plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
ax = fig7.add_subplot(2, 2, 3)
ax.plot(Ref_Race.Distance, Ref_Race.Iq, TT_Sim.Distance, TT_Sim.Iq)
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
ax = fig7.add_subplot(2, 2, 4)
ax.plot(Ref_Race.t, Ref_Race.Iq, TT_Sim.t, TT_Sim.Iq)
plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
fig7.show()

plt.show()
