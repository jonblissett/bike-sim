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
from motorbike_functions import motor_torque, motor_torque_speed, motor_sizing, battery_simple, motor_currents, \
    motor_saturation_interp, motor_current, charge_battery, motor_current_newton
import numpy as np
import scipy.io as sio

enable_plotting = True

sim = {}
sim['motor'] = {}
sim['motor']['L_core'] = 150
sim['motor']['N'] = 13.5
sim['motor']['W_max'] = 10500 / 30 * 3.14159
sim['motor']['W_lim'] = 10600 / 30 * 3.14159
sim['motor']['poles'] = 12

cell_ve = sio.loadmat('../data_import/cell_V_E.mat', struct_as_record=False, squeeze_me=True)['cell_ve']

sim['battery'] = {}
sim['battery']['cell_ve'] = {'v': np.array(cell_ve.v), 'e': np.array(cell_ve.e)}
sim['battery']['series'] = 120
sim['battery']['parallel'] = 6
sim['battery']['cellAh'] = 10
sim['battery']['cellVnom'] = 3.8

sim = charge_battery(sim, 1.096)
sim['constants'] = {}
sim['J'] = {}

sim['battery']['V_max'] = battery_simple(sim, 0, 0)[0]
sim = charge_battery(sim, 1.096)
sim = motor_sizing(sim)

sim = motor_saturation_interp(sim)
T100 = motor_current(sim['motor']['co'], 100)
T101 = motor_current_newton(sim['motor']['co'], 100, sim['motor']['T_con'], 0.01)
print('T=', T100, T101)
Ts = np.linspace(-300, 300, 500)
Is = motor_currents(sim, Ts)

print(sim['motor'])
print(sim['constants'])
print('co=', sim['motor']['co'])

if enable_plotting:
    import numpy as np
    import matplotlib.pyplot as plt
    I = np.linspace(0, sim['motor']['i_rms_pk'] * 2 ** 0.5, 100)
    T = motor_torque(sim['motor']['co'], I)
    fig, ax1 = plt.subplots()
    ax1.plot([], [], color='b', label='Torque with saturation', linewidth=5)
    ax1.plot([], [], color='g', label='Torque without saturation', linewidth=5)
    ax1.plot([], [], color='r', label='Torque interpolated', linewidth=5)
    ax1.plot(I, T, I, I*sim['motor']['co'][2], Is, Ts)
    ax1.set_xlabel('Current (A-pk)')
    ax1.set_ylabel('Torque (Nm)')
    plt.setp(ax1.get_xticklabels())
    plt.setp(ax1.get_yticklabels())
    ax1.legend(loc='lower right', fancybox=True)
    fig.show()

    # plt.savefig('../data_export/graphs/Python_saturation_150.pgf')

    [sim['motor']['w'], sim['motor']['t'], sim['motor']['p']] = motor_torque_speed(sim['motor']['T_pk'],
                                                                                   sim['motor']['W_max'],
                                                                                   sim['motor']['P_pk'],
                                                                                   sim['motor']['W_lim'], 50, True)
    plt.show()
