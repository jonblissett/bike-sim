import numpy as np
import scipy.io as sio
from scipy import signal

import time
from motorbike_functions import lap_analyse2, motor_torque_speed, wheel_forces

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
enable_plotting = False

if enable_plotting:
    try:
        import matplotlib.pyplot as plt
    except:
        pass

# Find list of N1&N2 to test with

# Select which corners to analyse
first_corner = 0    # 6#62  # 0 = first
last_corner = 99    # 7#63  # set to large number to use all

# Export parameters
filename_exp = 'Python_TT_sim_UoN01_max.mat'
structure_exp = 'TT_sim_UoN01_max'
# Import filenames
filename_ref_lap = 'TT_Laps_2016.mat'
filename_ref_map = 'TT_map_zerolean.mat'
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
TT_Sim = {'N': ([83.0, 18.0]),
          'constants': {'Km': 257 / 977, 'cd': 0.38, 'area': 1, 'rho': 1.204, 'm': 290.0 + 90 - 12.2, 'p_tyre': 1.9,
                        'r': 2.003 / 2 / np.pi, 'IRcell': 0.00438, 'b': 0.7, 'h': 0.6},
          'J': {'wheel': 1.35, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 1.615}}

# CHECK WHEEL INERTIA - IS THIS FRONT AND REAR OR TWO FRONT WHEELS?

# Model bike electrical specifications
#  motor torque vs speed characteristic
n_series = 120
rated_energy = n_series*3.8*60
initial_energy = rated_energy*1.096
turns = 11.5
L_core = 150
TT_Sim['motor'] = {}
TT_Sim['motor']['T_max'] = 977 * TT_Sim['constants']['Km']
TT_Sim['motor']['W_max'] = 10500 / 30 * np.pi
TT_Sim['motor']['W_lim'] = 10600 / 30 * np.pi
TT_Sim['motor']['P_max'] = 165e3  # 141000
TT_Sim['motor']['poles'] = 12
TT_Sim['motor']['Ld'] = 53e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Lq'] = 61e-6 * (turns/11.5)**2 * L_core / 150
TT_Sim['motor']['Rs'] = 0.007313 * (turns/11.5) * L_core / 150
TT_Sim['motor']['Ke'] = 0.3467/6 * (turns/11.5) * L_core / 150  # Back-emf constant, V/erad/s
[TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = motor_torque_speed(TT_Sim['motor']['T_max'],
                                                                                        TT_Sim['motor']['W_max'],
                                                                                        TT_Sim['motor']['P_max'],
                                                                                        TT_Sim['motor']['W_lim'],
                                                                                        50, enable_plotting)
# TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
# TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']
# SETTING LIKE THIS DOESN@T WORK AS THEY ARE RECALCULATED

# END OF USER PARAMETERS

timer = time.time()

Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
v = 1.0*r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
del r
del N2
del N1  # Just to ensure they are not used accidentally

TT_Sim = lap_analyse2(TT_Sim, Ref_Race, v, first_corner, last_corner, filename_ref_map, filename_ref_brake,
                      structure_map, var_name_brake, rated_energy, initial_energy, n_series, enable_warnings)

print('Simulation duration =', str(time.time() - timer), 'seconds')
sio.savemat(filename_exp, {structure_exp: TT_Sim}, oned_as='column')
print('Simulation saved to', filename_exp, 'as structure named', structure_exp)

print('Motor max speed', str(max(TT_Sim.v) / TT_Sim.constants.r * 30 / np.pi * TT_Sim.N[0] / TT_Sim.N[1]))
print('Bike max speed (mph) = ', str(max(TT_Sim.v) * 2.23))
print('Simulated lap time (s) = ', str(TT_Sim.t[-1]))
print('Simulated lap speed (mph) = ', str(37.733/TT_Sim.t[-1]*3600))

if enable_plotting:
    wheel_forces(1.415, 0.6, 0.7, TT_Sim.constants.r, TT_Sim.constants.m,
                 np.square(TT_Sim.v) * TT_Sim.constants.rho * TT_Sim.constants.cd * TT_Sim.constants.area / 2.0,
                 TT_Sim.torque * TT_Sim.N[0] / TT_Sim.N[1], 1)

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

    b, a = signal.butter(3, 0.1)


    fig, ax1 = plt.subplots()
    ax1.plot([], [], color='b', label='Speed (Race)', linewidth=5)
    ax1.plot([], [], color='g', label='Speed (Simulated)', linewidth=5)
    ax1.plot(Ref_Race.Distance, signal.filtfilt(b, a, v), TT_Sim.Distance, TT_Sim.v)
    ax1.set_xlabel('Distance (m)', fontsize=18)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('Speed (m/s)', fontsize=18)
    # for tl in ax1.get_yticklabels():
    #     tl.set_color('b')
    plt.xlim(TT_Sim.Distance[0]-25, TT_Sim.Distance[-1])

    ax2 = ax1.twinx()
    ax2.plot([], [], color='c', label='Torque (Race)', linewidth=5)
    ax2.plot([], [], color='r', label='Torque (Simulated)', linewidth=5)
    ax2.plot(Ref_Race.Distance, signal.filtfilt(b, a, Ref_Race.Iq)*TT_Sim.constants.Km, 'c-', TT_Sim.Distance, TT_Sim.Iq*TT_Sim.constants.Km, 'r-', zorder=1)
    ax2.set_ylabel('Motor Torque (Nm)', fontsize=18)
    # for tl in ax2.get_yticklabels():
    #    tl.set_color('r')
    # plt.title('Motor Torque and Power vs Speed')
    plt.xlim(TT_Sim.Distance[0]-25, TT_Sim.Distance[-1])
    plt.setp(ax1.get_xticklabels(), fontsize=16)
    plt.setp(ax2.get_xticklabels(), fontsize=16)
    plt.setp(ax1.get_yticklabels(), fontsize=16)
    plt.setp(ax2.get_yticklabels(), fontsize=16)
    ax1.legend(loc='lower left', fancybox=True)
    ax2.legend(loc='lower right', fancybox=True)
    # plt.show()
    fig.savefig('TT_sim_model.png')
    fig.show()

    plt.show()
