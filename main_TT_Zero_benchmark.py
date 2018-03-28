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
import scipy.io as sio
from scipy import signal
from os import remove, rename, path
from itertools import product
from itertools import izip_longest  # This is used to deal with variable length of lists
import time
import motorbike_functions as bike

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

verbosity = 0  # 0, no print, 1, final stats, 2, per corner stats and warnings, 3 everything
enable_warnings = False
enable_plotting = False
enable_parallel = True  # ipcluster start -n 4
save_data_files = True
dummy_run = False
calibration_mode = False
optimise_ratio = False
battery_fixed = False
fake_parallel = False
motor_manufacturer = 'Parker'  # 'me', 'Emrax'
igbt = 'SEMiX603_SiC'  # 'FF600'  # 'SEMiX603_SiC'

course_speed_limit = 200 / 2.23

track = 'TT'

if enable_plotting:
    try:
        import matplotlib.pyplot as plt
    except:
        pass


parallel_queue = 600

# Model bike mechanical specifications
#  Model Brake characteristic
#  - Release time of braking instance
#  - Regenerative torque limit
#  - Braking co-efficient, relating braking torque to w
TT_Sim = {'N': ([71.0, 18.0]),
          'constants': {'cd': 0.32, 'area': 1, 'rho': 1.204, 'm': 290.0 + 90, 'p_tyre': 1.9,
                        'r': 2.16 / 2 / np.pi, 'b': 0.725, 'h': 0.56, 'k_tyre': 0.7, 'mu_tyre': 1.2},
          'J': {'wheel': 1.35 - 0.445, 'motor': 0.0233},
          'brake': {'RampTime': 1.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 0},
          #'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 1.615},
          'battery': {},
          'motor': {'manufacturer': motor_manufacturer},
          'drive': {},
          'IGBT': bike.spec_igbt(igbt),
          'v_max': course_speed_limit,
          'file': {'motorimport': 'MotorLAB_export.mat'}#_Mr25
          }

if track == 'TT':
    # Select which corners to analyse
    first_corner = 0  # 6#62  # 0 = first
    last_corner = 999  # 7#63  # set to large number to use all
    corner_delete = [11, 60, 93, 106, 107]
    laps = 1

    end_dist = 60.7e3  # Distance at lap end for timing

    # Export parameters
    filename_exp = 'data_export/Python_Sims_FW_TT_2017.mat'
    # filename_exp = 'data_export/TT_SpeedLimits/Python_Sims_FW_TT_SpeedLimit_' + str(int(2.23*course_speed_limit)) + 'mph.mat'
    # filename_exp = 'data_export/18s16p_P_T_9_various_Mr25/Python_Sim_' + motor_manufacturer + '_motor_power_varied_mph_Mr25_regen.mat'

    structure_exp = 'TT_sims'

    # Import filenames
    if enable_parallel:
        #filename_ref_lap = 'data_import/TT_Race_2016_small.mat'
        filename_ref_lap = 'data_import/TT_Race_Louis.mat'
    else:
        filename_ref_lap = 'data_import/TT_Race_Louis.mat'
        # filename_ref_lap = 'data_import/TT_Laps_2016.mat'
    filename_ref_map = 'data_import/TT_map_zerolean.mat'
    # filename_ref_brake = 'data_import/TT_Race_2016_manual_braking_pts.mat'
    filename_ref_brake = 'data_import/TT_Race_Louis_manual_braking_pts.mat'
    filename_ref_bat = 'data_import/Python_Sims_FW_TT_guessed_noregen_120_vlim.mat'
    # filename_motor = 'data_import/MotorLAB_export.mat'
    structure_map = 'TT_map'
    structure_lap = 'TT_Race'
    var_name_brake = 'locsMin'
    structure_bat = 'TT_Sim'

    # Reference data mechanical specifications
    N1 = 83.0
    N2 = 19.0
    r = 2.16 / 2 / np.pi
    m_ref = 270.0 + 90
    # m_ref = 240
    # wMotor_ref = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
    # TMotor_ref = np.array([106, 106, 65.791, 45.949])
    Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
    # v_ref = 1.0 * r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
    v_ref = Ref_Race.vGPS
    TT_Sim['scrutineering'] = {'score': 0.0, 'weight_limit': 305.0, 'volt_limit': 800.0}
    TT_Sim['battery']['series'] = 162
    TT_Sim['battery']['parallel'] = 4
    TT_Sim['battery']['cellAh'] = 10  # -0.28
    TT_Sim['battery']['cellVnom'] = 3.7
    TT_Sim['battery']['cellIR'] = 0.00438
    TT_Sim['battery']['E_density'] = 3.7 * 6 * 40 / 4.8  # 3.7*8/0.175
    charge_ratio = 1.04  # 1.096
    motor_mass_kf = 1
    scrutineering_cheat = 4

    if battery_fixed:
        sim = sio.loadmat(filename_ref_bat, struct_as_record=False, squeeze_me=True)[structure_bat]
        TT_Sim['Vdc_sim'] = sim.Vdc_sim

    variables_list = {
        'P_max': np.arange(100, 210, 5) * 1000,
        'T_max': np.arange(180, 300, 10),
        # 'n0': range(42, 84, 41),
        'n1': range(17, 25, 1),
        'v_max': np.arange(159, 176, 1) / 2.23,
        # 'parallel': np.arange(4.5, 7, 0.5),
        # 'L_core': range(100, 525, 25),
        # 'L_core': np.arange(150, 200, 25),
        # 'turns': np.arange(8.5, 16.5, 2),
        # 'drives': range(1, 3),
        # CdA ??
        # regen?
    }
elif track == 'Portimao':
    # Select which corners to analyse
    first_corner = 0
    last_corner = 16
    corner_delete = []
    laps = 6

    end_dist = 4494 * laps - 250 - 230 # Distance at lap end for timing
    print(' CHECK RADIUS WITH TRACK LENGTH')

    # Export parameters
    filename_exp = 'data_export/Python_Sims_FW_P_HalfMmotor'
    structure_exp = 'Portimao_sims'

    # Import filenames
    if enable_parallel:
        filename_ref_lap = 'data_import/Portimao_RaceDM2_small.mat'
    else:
        filename_ref_lap = 'data_import/Portimao_RaceDM2.mat'
    filename_ref_map = 'data_import/Portimao_map.mat'
    filename_ref_brake = 'data_import/Portimao_RaceDM2.mat'
    structure_map = 'Portimao_map'
    structure_lap = 'RaceDM2'
    var_name_brake = 'locsMin'

    # Reference data mechanical specifications
    N1 = 83.0
    N2 = 16.0
    r = 2.16 / 2 / np.pi
    # m_ref = 270.0 + 90
    m_ref = 240-3 + 90
    Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
    v_ref = Ref_Race.vGPS
    TT_Sim['scrutineering'] = {'score': 0.0, 'weight_limit': 250.0 * 1.01, 'volt_limit': 600.0}
    TT_Sim['battery']['series'] = 168
    TT_Sim['battery']['parallel'] = 2
    TT_Sim['battery']['cellAh'] = 8
    TT_Sim['battery']['cellVnom'] = 3.7
    TT_Sim['battery']['cellIR'] = 0.00219 # 0.001 # 0.00438
    TT_Sim['battery']['E_density'] = 169  #169 148  # 3.8 * 6 * 40 / 4.8
    charge_ratio = 1.0
    motor_mass_kf = 1
    scrutineering_cheat = 4

    variables_list = {
        'P_max': np.arange(390, 400, 10) * 1000,
        'T_max': np.arange(220, 1000, 10),
        #'n0': range(42, 84, 41),
        #'n1': range(15, 43, 1),
        'parallel': np.arange(3.0, 3.5, 0.5),
        'L_core': range(100, 525, 50),
        'turns': np.arange(8.5, 22.5, 2),
        'drives': range(1, 3),
        # CdA ??
        # regen?
    }
elif track == 'Drag':
    # Select which corners to analyse
    first_corner = 0
    last_corner = 1
    laps = 1

    end_dist = 1609.3 # Distance at lap end for timing

    # Export parameters
    filename_exp = 'data_export/Python_sim_Drag.mat'
    structure_exp = 'Drag_sim'

    # Import filenames
    filename_ref_lap = 'data_import/Drag_Race_pretend_2.mat'
    filename_ref_map = 'data_import/Drag_Race_pretend_2.mat'
    filename_ref_brake = 'data_import/Drag_Race_pretend_2.mat'
    structure_map = 'Drag_Map'
    structure_lap = 'Drag_Race'
    var_name_brake = 'locsMin'

    TT_Sim['N'] = ([83.0, 21.0])

    # Reference data mechanical specifications
    N1 = 83.0
    N2 = 19.0
    r = 2.003 / 2 / np.pi
    # m_ref = 270.0 + 90
    m_ref = 240 - 3
    Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
    v_ref = Ref_Race.v
    TT_Sim['scrutineering'] = {'score': 0.0, 'weight_limit': 250.0 * 1.01, 'volt_limit': 600.0}
    TT_Sim['battery']['series'] = 120
    TT_Sim['battery']['parallel'] = 2
    TT_Sim['battery']['cellAh'] = 8
    TT_Sim['battery']['cellVnom'] = 3.8
    TT_Sim['battery']['cellIR'] = 0.00219  # 0.001 # 0.00438
    TT_Sim['battery']['E_density'] = 169  # 148  # 3.8 * 6 * 40 / 4.8
    charge_ratio = 1.0
    motor_mass_kf = 1.0
    scrutineering_cheat = 4

    variables_list = {
        'P_max': np.arange(390, 400, 10) * 1000,
        'T_max': np.arange(220, 1000, 10),
        'n0': range(42, 84, 41),
        'n1': range(15, 43, 1),
        'parallel': np.arange(3.0, 3.5, 0.5),
        'L_core': range(100, 525, 50),
        'turns': np.arange(8.5, 22.5, 2),
        'drives': range(1, 3),
        # CdA ??
        # regen?
    }

# Model bike electrical specifications
cell_ve = sio.loadmat('data_import/cell_V_E.mat', struct_as_record=False, squeeze_me=True)['cell_ve']

TT_Sim['battery']['cell_ve'] = {'v': np.array(cell_ve.v), 'e': np.array(cell_ve.e)}
TT_Sim['battery']['IR'] = TT_Sim['battery']['cellIR'] * TT_Sim['battery']['series'] / TT_Sim['battery']['parallel']
TT_Sim = bike.charge_battery(TT_Sim, charge_ratio)
TT_Sim['battery']['V_max'] = bike.battery_simple(TT_Sim, 0, verbosity)[0]
#################################

if TT_Sim['motor']['manufacturer'] == 'Parker':
    TT_Sim['motor']['N'] = 18.5  # p is 18.5
    TT_Sim['motor']['L_core'] = 150.0

TT_Sim = bike.motor_sizing(TT_Sim)
# TT_Sim = bike.set_speed_limit(TT_Sim, TT_Sim['v_max'])
w_max = 6400 /30 * np.pi # 130 / 2.23 * TT_Sim['N'][0] / TT_Sim['N'][1] / TT_Sim['constants']['r']
TT_Sim['motor']['W_speed_lim'] = 0.99 * w_max
TT_Sim['motor']['W_lim'] = 7100 / 30 * np.pi # TT_Sim['motor']['W_speed_lim'] * 1.15  # 8400 / 30 * np.pi
print('SPEED LIMIT SET GENTLY')

TT_Sim['motor']['T_max'] = bike.motor_torque(TT_Sim['motor']['co'], 452*0.6)  # TT_Sim['motor']['T_pk']  #
print('Motor torque = ' + str(TT_Sim['motor']['T_max']) + ' Nm')
TT_Sim['motor']['P_max'] = 150e3  # TT_Sim['motor']['P_pk'] * 0.99  # *0.66
[TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = bike.motor_torque_speed(TT_Sim['motor']['T_max'],
                                                                                             TT_Sim['motor']['W_speed_lim'],
                                                                                             TT_Sim['motor']['P_max'],
                                                                                             TT_Sim['motor']['W_lim'],
                                                                                             50, enable_plotting)

TT_Sim['drive']['n'] = 1.0
TT_Sim['drive']['m'] = 5.0
TT_Sim['drive']['I_max'] = 1001  # THIS IS ONLY FOR SCRUTINEERING FUNCTION

# TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
# TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']
# SETTING LIKE THIS DOESN@T WORK AS THEY ARE RECALCULATED

# END OF USER PARAMETERS

v = bike.v_tyre_load_sens(v_ref, m_ref, TT_Sim['constants']['k_tyre'], TT_Sim['constants']['m'])

del r
del N2
del N1  # Just to ensure they are not used accidentally


def save_parallel(filename, async):
    try:
        results = np.asarray([data.get() for data in async])
        Sims = {}
        Sims['ALL'] = np.array(results)
        Sims['names'] = 'Pmax,Tmax,Eused,Erated,Einit,tmax,N1,N2,Vend,m,dTbat,vmax,RPMmax,ERROR,Lcore,turns,Ns,Np,drives'
        Sims['Pmax'] = np.array(results[:, 1])
        Sims['Tmax'] = np.array(results[:, 2])
        Sims['Eused'] = np.array(results[:, 3])
        Sims['Erated'] = np.array(results[:, 4])
        Sims['Einit'] = np.array(results[:, 5])
        Sims['tmax'] = np.array(results[:, 6])
        Sims['N1'] = np.array(results[:, 7])
        Sims['N2'] = np.array(results[:, 8])
        Sims['Vend'] = np.array(results[:, 9])
        Sims['m'] = np.array(results[:, 10])
        Sims['dTbat'] = np.array(results[:, 11])
        Sims['vmax'] = np.array(results[:, 12])
        Sims['RPMmax'] = np.array(results[:, 13])
        Sims['ERROR'] = np.array(results[:, 14])
        Sims['Lcore'] = np.array(results[:, 15])
        Sims['turns'] = np.array(results[:, 16])
        Sims['Ns'] = np.array(results[:, 17])
        Sims['Np'] = np.array(results[:, 18])
        Sims['drives'] = np.array(results[:, 19])
        sio.savemat(filename, {'Sims': Sims}, oned_as='column')
        print('Saved results to filename: ' + filename)
        return 0
    except ZeroDivisionError:
        print('####################################')
        print('##### FAILED TO SAVE SOME DATA #####')
        print('######FILE:' + filename + '#####')
        print('####################################')
        return 1


if enable_parallel:
    n_sims = 1.0
    for a, b in variables_list.iteritems():
        n_sims *= len(b)
    print('Number of sims expected:', n_sims)

    dict1 = [dict(izip_longest(variables_list, va)) for va in product(*variables_list.values())]

    from ipyparallel import Client

    full_data_exp = False
    # rc = Client(profile='/home/jpb/.ipython/profile_ssh/pid/ipcluster.pid')
    # rc = Client('/home/jpb/.ipython/profile_ssh/security/ipcontroller-client.json', sshserver='jpb@80.229.16.48:9989')
    # rc = Client(profile='ssh')
    if not dummy_run and not fake_parallel:
        rc = Client()
        print('Cores assigned', rc.ids)
        with rc[:].sync_imports():
            import motorbike_functions
        view = rc.load_balanced_view()
        view.block = True
        if save_data_files:
            sio.savemat(filename_exp + '_variables_list.mat', {'list': variables_list}, oned_as='column')
            sio.savemat(filename_exp + '_constants_list.mat', {'settings': TT_Sim}, oned_as='column')

    async_results = []

    timer = time.time()

    count = 0
    count_all = 0
    save_ID = 0
    prev_m = 0
    for dict2 in dict1:
        count_all += 1
        prev_N = TT_Sim['motor']['N']
        prev_L = TT_Sim['motor']['L_core']
        prev_m = TT_Sim['constants']['m']
        prev_v = TT_Sim['battery']['V_max']
        # print('Progress %.1f%%' % (count / n_sims * 100))
        for key, value in dict2.items():
            if key in TT_Sim:
                TT_Sim[key] = value
            else:
                if key in TT_Sim['battery']:
                    TT_Sim['battery'][key] = value * 1.0
                if key in TT_Sim['motor']:
                    TT_Sim['motor'][key] = value * 1.0
                if key in TT_Sim['constants']:
                    TT_Sim['constants'][key] = value * 1.0
                if key in TT_Sim['brake']:
                    TT_Sim['brake'][key] = value * 1.0
                if key == 'n0':
                    TT_Sim['N'][0] = value * 1.0
                if key == 'n1':
                    TT_Sim['N'][1] = value * 1.0
                if key == 'turns':
                    TT_Sim['motor']['N'] = value * 1.0
                if key == 'drives':
                    TT_Sim['drive']['n'] = value * 1.0
                # If key not found??
        TT_Sim = bike.set_speed_limit(TT_Sim, TT_Sim['v_max'])
        if (TT_Sim['motor']['L_core'] != prev_L) or (TT_Sim['motor']['N'] != prev_N):
            TT_Sim = bike.charge_battery(TT_Sim, charge_ratio)
            TT_Sim['battery']['V_max'] = bike.battery_simple(TT_Sim, 0, 0)[0]
            TT_Sim = bike.motor_sizing(TT_Sim)
            TT_Sim['motor']['m'] *= motor_mass_kf
        TT_Sim = bike.scrutineering(TT_Sim, charge_ratio)
        if TT_Sim['scrutineering']['passed'] or TT_Sim['scrutineering']['score'] == scrutineering_cheat:
            count += 1
            if TT_Sim['constants']['m'] != prev_m:
                v = bike.v_tyre_load_sens(v_ref, m_ref, TT_Sim['constants']['k_tyre'], TT_Sim['constants']['m'])
            if not dummy_run:
                if not fake_parallel:
                    if optimise_ratio:
                        ar = view.apply_async(bike.gear_optimise, TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps,
                                              end_dist, filename_ref_map, filename_ref_brake, structure_map, var_name_brake,
                                              enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed)
                    else:
                        ar = view.apply_async(bike.lap_analyse3, TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps,
                                              end_dist, filename_ref_map, filename_ref_brake, structure_map, var_name_brake,
                                              enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed)
                else:
                    if optimise_ratio:
                        ar = bike.gear_optimise(TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps, end_dist,
                                                filename_ref_map, filename_ref_brake, structure_map, var_name_brake,
                                                enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed)
                    else:
                        ar = bike.lap_analyse3(TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps, end_dist,
                                               filename_ref_map, filename_ref_brake, structure_map, var_name_brake,
                                               enable_warnings, verbosity, calibration_mode, full_data_exp, battery_fixed)
                async_results.append(ar)
                if count % parallel_queue == 0:
                    print('Waiting for results, %.1f%% completed' % (count_all / n_sims * 100))
                    if not fake_parallel:
                        rc.wait_interactive(async_results)
                        if count % (parallel_queue * 1) == 0:
                            save_ID += 1
                    if save_data_files:
                        # save_parallel(filename_exp + 'teMP_' + str(count), async_results)
                        s = save_parallel(filename_exp + '_' + str(int(round((count_all / n_sims * 1000)))) +
                                          '_part_' + str(save_ID), async_results)
                        del async_results
                        async_results = []
            print('Race ' + str(count) + ' started: m = ' + str(round(TT_Sim['mass']['bike'])), TT_Sim['motor']['P_max'],
                  TT_Sim['motor']['T_max'], TT_Sim['N'],
                  str(TT_Sim['battery']['series']) + 's' + str(TT_Sim['battery']['parallel']) + 'p',
                  str(TT_Sim['motor']['name']) + '-X' + str(TT_Sim['drive']['n']), str(TT_Sim['v_max']*2.23) + ' mph')
        # else:
        #     print('Scrutineering FAILED, score: ' + str(TT_Sim['scrutineering']['score']))

    #print('ET=', 0.7 * count / 60, ' minutes')
    print('ET=', 3.97 * count / 3600, ' hours')

    if not dummy_run:
        print('Waiting for results, %.1f%% completed' % (count_all / n_sims * 100))
        if not fake_parallel:
            print(rc.wait_interactive(async_results))  # Wait until all tasks are done
        print('Completed', len(async_results), 'simulations in', time.time() - timer, 'seconds')
        save_ID += 1
        if save_data_files:
            # save_parallel(filename_exp, async_results)
            s = save_parallel(filename_exp + '_' + str(int(round((count_all / n_sims * 1000)))) + '_part_' +
                              str(save_ID), async_results)
            # Sims = sio.loadmat(filename_exp, struct_as_record=False, squeeze_me=True)['Sims']

    #for a in range(0, len(async_results)):
    #    print(async_results[a].stdout.encode('utf-8'))
    print('Completed', count_all, 'simulations in', time.time() - timer, 'seconds')

else:
    full_data_exp = True
    timer = time.time()
    TT_Sim['battery']['V_max'] = bike.battery_simple(TT_Sim, 0, verbosity)[0]
    TT_Sim = bike.motor_sizing(TT_Sim)
    TT_Sim['motor']['m'] *= motor_mass_kf
    TT_Sim = bike.scrutineering(TT_Sim, charge_ratio)
    if TT_Sim['scrutineering']['passed']:
        print('Scrutineering passed')
    else:
        print('Scrutineering FAILED, score: ' + str(TT_Sim['scrutineering']['score']))
    v = bike.v_tyre_load_sens(v_ref, m_ref, TT_Sim['constants']['k_tyre'], TT_Sim['constants']['m'])

    print('Race ' + 'started: m = %.1f' % (TT_Sim['mass']['bike']) + 'kg, %.0fkWp ' % (TT_Sim['motor']['P_max'] / 1000),
          TT_Sim['motor']['T_max'], TT_Sim['N'],
          str(TT_Sim['battery']['series']) + 's' + str(TT_Sim['battery']['parallel']) + 'p',
          str(TT_Sim['motor']['name']) + '-X' + str(TT_Sim['drive']['n']))

    if optimise_ratio:
        TT_Sim = bike.gear_optimise(TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map,
                                    filename_ref_brake, structure_map, var_name_brake, enable_warnings,
                                    verbosity, calibration_mode, full_data_exp, battery_fixed)
    else:
        TT_Sim = bike.lap_analyse3(TT_Sim, Ref_Race, v, first_corner, last_corner, corner_delete, laps, end_dist, filename_ref_map,
                                   filename_ref_brake, structure_map, var_name_brake, enable_warnings, verbosity,
                                   calibration_mode, full_data_exp, battery_fixed)
    print('Simulation duration = %.2f s' % (time.time() - timer))
    # Bit of a kludge to convert dict type to mat_structure
    # THIS COULD BE ELIMINATED BY USING DICTS THROUGHOUT
    sio.savemat('temp', {'TT_Sim': TT_Sim}, oned_as='column')
    TT_Sim = sio.loadmat('temp', struct_as_record=False, squeeze_me=True)['TT_Sim']

    if save_data_files:
        if path.exists(filename_exp):
            remove(filename_exp)
        rename('temp.mat', filename_exp)
        if verbosity > 0:
            print('Simulation saved to ' + str(filename_exp) + ' as structure named ' + str(structure_exp))
    else:
        remove('temp.mat')
        if verbosity > 0:
            print('Simulation NOT saved')

    if verbosity > 0:
        end_dista = TT_Sim.Distance_race < end_dist
        tmax = TT_Sim.t[sum(end_dista)-1]
        print('Estimated bike mass = %.1f kg' % (TT_Sim.constants.m - 90))
        print('Motor max speed = %d rpm' % max(TT_Sim.Rpm))
        print('Motor mass = %.1f kg' % TT_Sim.motor.m + ' + inertial mass of %.2f kg' %
              (((TT_Sim.N[0] / TT_Sim.N[1]) ** 2 * TT_Sim.J.motor)/(TT_Sim.constants.r ** 2)))
        print('Bike max speed = %.1f mph' % (max(TT_Sim.v) * 2.23))
        print('Simulated lap time = %.2f s' % tmax)
        print('Simulated lap speed = %.2f mph' % (37.733 / TT_Sim.t[-1] * 3600))

    if enable_plotting:
        bike.wheel_forces(TT_Sim.Distance,1.415, TT_Sim.constants.h, TT_Sim.constants.b, TT_Sim.constants.R,
                          TT_Sim.constants.m, np.square(TT_Sim.v) * TT_Sim.constants.rho * TT_Sim.constants.cd
                          * TT_Sim.constants.area / 2.0, TT_Sim.torque * TT_Sim.N[0] / TT_Sim.N[1], 1)

        fig7 = plt.figure(7)
        ax = fig7.add_subplot(2, 2, 1)
        ax.plot(Ref_Race.Distance, v, TT_Sim.Distance_race, TT_Sim.v)
        plt.xlim(TT_Sim.Distance_race[0], end_dist)
        ax = fig7.add_subplot(2, 2, 2)
        ax.plot(Ref_Race.t, v, TT_Sim.t, TT_Sim.v)
        plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
        ax = fig7.add_subplot(2, 2, 3)
        ax.plot(Ref_Race.Distance, Ref_Race.Iq, TT_Sim.Distance_race, TT_Sim.Iq)
        plt.xlim(TT_Sim.Distance_race[0], end_dist)
        ax = fig7.add_subplot(2, 2, 4)
        ax.plot(Ref_Race.t, Ref_Race.Iq, TT_Sim.t, TT_Sim.Iq)
        plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
        fig7.show()

        b, a = signal.butter(3, 0.1)

        fig, ax1 = plt.subplots()
        # ax1.plot([], [], color='b', label='Speed (Race)', linewidth=5)
        ax1.plot([], [], color='g', label='Speed (Simulated)', linewidth=5)
        ax1.plot(Ref_Race.Distance, signal.filtfilt(b, a, v), TT_Sim.Distance, TT_Sim.v)
        ax1.set_xlabel('Distance (m)', fontsize=18)
        # Make the y-axis label and tick labels match the line color.
        ax1.set_ylabel('Speed (m/s)', fontsize=18)
        # for tl in ax1.get_yticklabels():
        #     tl.set_color('b')
        plt.xlim(TT_Sim.Distance[0] - 25, TT_Sim.Distance[-1])

        ax2 = ax1.twinx()
        # ax2.plot([], [], color='c', label='Torque (Race)', linewidth=5)
        ax2.plot([], [], color='r', label='Torque (Simulated)', linewidth=5)
        #ax2.plot(Ref_Race.Distance, signal.filtfilt(b, a, Ref_Race.Iq) * Ref_Race.constants.Km, 'c-', TT_Sim.Distance,
        #         bike.motor_torque(TT_Sim.motor.co, TT_Sim.Iq), 'r-', zorder=1)
        ax2.plot(TT_Sim.Distance, bike.motor_torque(TT_Sim.motor.co, TT_Sim.Iq), 'r-', zorder=1)
        ax2.set_ylabel('Motor Torque (Nm)', fontsize=18)
        # for tl in ax2.get_yticklabels():
        #    tl.set_color('r')
        # plt.title('Motor Torque and Power vs Speed')
        plt.xlim(TT_Sim.Distance[0] - 25, TT_Sim.Distance[-1])
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
