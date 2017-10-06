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
import motorbike_functions as bike
from C_pmsm import C_id_pmsm, C_w_pmsm, C_vs_pmsm, C_v_dq_pmsm, C_motor_current_newton
from C_losses import C_motor_losses, C_inverter_losses, C_inverter_loss


import matplotlib.pyplot as plt

motor_manufacturer = 'Parker'  # 'me', 'Emrax'
course_speed_limit = 200 / 2.23

filename = './data_export/Python_sim_motor_fw_test.mat'

TT_Sim = {'N': ([71.0, 18.0]),
          'constants': {'cd': 0.32, 'area': 1, 'rho': 1.204, 'm': 290.0 + 90, 'p_tyre': 1.9,
                        'r': 2.16 / 2 / np.pi, 'b': 0.725, 'h': 0.56, 'k_tyre': 0.7, 'mu_tyre': 1.2},
          'J': {'wheel': 1.35 - 0.445, 'motor': 0.0233},
          'brake': {'RampTime': 1.6, 'PeakTorque': 1100.0, 'LimitTorque': 300.0, 'k_wt': 0},
          #  'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 300.0, 'k_wt': 1.615},
          'battery': {},
          'motor': {'manufacturer': motor_manufacturer},
          'drive': {},
          'IGBT': bike.spec_igbt(igbt),
          'v_max': course_speed_limit,
          'file': {'motorimport': 'MotorLAB_export.mat'}  #_Mr25
          }

TT_Sim['motor']['N'] = 18.5  # p is 18.5
TT_Sim['motor']['L_core'] = 150.0

TT_Sim['battery']['IR'] = 0.0
TT_Sim['battery']['V_max'] = 687

TT_Sim['drive']['n'] = 1

TT_Sim = bike.motor_sizing(TT_Sim)


Ke = TT_Sim['motor']['Ke']
Poles = TT_Sim['motor']['poles']
Rs = TT_Sim['motor']['Rs']
Ld = TT_Sim['motor']['Ld']
Lq = TT_Sim['motor']['Lq']
IR = TT_Sim['battery']['IR']
L_core = TT_Sim['motor']['L_core']
T_con = TT_Sim['motor']['T_con']
driven = TT_Sim['drive']['n']
co = TT_Sim['motor']['co']

i_smax = 400*2**0.5
vdc = TT_Sim['battery']['V_max']
# w_m = 8000*np.pi/30
i_q0 = i_smax
i_d0 = 0
w_m_list = np.arange(0,11000,10)*np.pi/30
i_d_list = w_m_list * 0
i_q_list = w_m_list * 0
v_bus_list = w_m_list * 0
j = 0
m = 0.9    # modulation limit
IR = 0.18
i_d = i_d0

for w_m in w_m_list:
    iterate = True
    count = 0
    i_q = i_q0
    while iterate:
        count += 1
        P_mech = w_m * bike.motor_torque(TT_Sim['motor']['co'],i_q)
        motor_loss = C_motor_losses(np.sqrt(i_d ** 2 + i_q ** 2), w_m * 30 / np.pi, Rs, TT_Sim['motor']['k_rpm'][2],
                                        TT_Sim['motor']['k_rpm'][1])[0]
        [vs, vd, vq, PF] = C_v_dq_pmsm(Ke, Poles, Rs, Ld, Lq, i_d, i_q, w_m)
        drive_loss= driven * C_inverter_loss(vdc, vs, np.sqrt(i_d ** 2 + i_q ** 2) / np.sqrt(2) / driven, PF, 82e-6, 13e3, 0.8, 1,
                                                    0.95e-3, 0.54e-3, 12e-3, 25e-3, 9.5e-3)[0]
        P = P_mech + motor_loss + drive_loss
        # print(P_mech, motor_loss, drive_loss)

        idc = P / vdc
        vir = idc * IR
        vbus = vdc - vir

        i_d = C_id_pmsm(Ke, Poles, Rs, Ld, Lq, i_q, w_m, vbus * 0.866 * m)
        i_s = np.sqrt(i_d ** 2 + i_q ** 2)
        i_s_error = i_smax - i_s

        if i_d < 0:
            i_q = i_q + i_s_error
            iterate = abs(i_s_error) > 1
        #print(count, w_m, i_d, i_q, i_s, i_smax)
        else:
            i_d = 0
            if np.isnan(i_s):
                i_q = 0
            iterate = False
        if count > 100:
            print('ITERATION STUCK ' + str(i_s) + ' ' + str(i_smax))
            iterate = False
    i_d_list[j] = i_d
    i_q_list[j] = i_q
    v_bus_list[j] = vbus
    j += 1
    print(count, w_m, i_d, i_q, i_s, i_smax)


torque = bike.motor_torque(TT_Sim['motor']['co'],i_q_list)
power = torque * w_m_list

fwsim = {}
fwsim['torque'] = torque
fwsim['power'] = power
fwsim['w'] = w_m_list
fwsim['rpm'] = w_m_list * 30.0 / np.pi
fwsim['Id'] = i_d_list
fwsim['Iq'] = i_q_list
fwsim['Is'] = (i_d_list ** 2 + i_q_list ** 2) ** 0.5
fwsim['Vdc'] = v_bus_list

sio.savemat(filename, {'fwsim': fwsim}, oned_as='column')
print('Saved results to filename: ' + filename)


fig, ax1 = plt.subplots()
ax1.plot(w_m_list * 30.0 / np.pi, torque, 'b-')
ax1.set_xlabel('Angular Speed (rpm)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Torque (Nm)', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(w_m_list * 30.0 / np.pi, power / 1000, 'g-')
ax2.set_ylabel('Power (kW)', color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')
plt.title('Motor Torque and Power vs Speed')

fig.show()

print(TT_Sim['motor']['Ke'])

fig, ax1 = plt.subplots()
ax1.plot(w_m_list * 30.0 / np.pi, i_d_list*0.707, 'b-', w_m_list * 30.0 / np.pi, i_q_list*0.707, 'g-')
ax1.set_xlabel('Angular Speed (rpm)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Id (Arms)', color='b')
#for tl in ax1.get_yticklabels():
#    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(w_m_list * 30.0 / np.pi, v_bus_list, 'r-')
ax2.set_ylabel('Vbus (V)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.title('Motor Torque and Power vs Speed')

fig.show()


plt.show()

# xrdb -merge ~/.Xresources