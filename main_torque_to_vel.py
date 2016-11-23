import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import integrate
from os import remove

import time
from motorbike_functions import lap_analyse, motor_torque_speed, motorbike_mech4
from scipy.integrate import odeint
from scipy import integrate
from scipy.integrate import ode
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
c = 0
# Export parameters
filename_exp = 'Python_sim_torque_to_vel.mat'
structure_exp = 'Portimao_sim'
# Import filenames
filename_ref_lap = 'Portimao_RaceDM2.mat'
filename_ref_map = 'Portimao_map.mat'
filename_ref_brake = 'Portimao_RaceDM2.mat'
structure_map = 'Portimao_map'
structure_lap = 'RaceDM2'
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
TT_Sim = {'N': ([83.0, 16.0]),
          'constants': {'Km': 257/977, 'cd': 0.38, 'area': 1, 'rho': 1.204, 'm': 235 + 90, 'p_tyre': 1.9,
                        'r': 2.18 / 2 / np.pi},
          'J': {'wheel': 1.1, 'motor': 0},
          'brake': {'RampTime': 2.6, 'PeakTorque': 830.0, 'LimitTorque': 200.0, 'k_wt': 1.15}}
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
# [TT_Sim['motor']['w'], TT_Sim['motor']['t'], TT_Sim['motor']['p']] = motor_torque_speed(TT_Sim['motor']['T_max'],
#                                                                                        TT_Sim['motor']['W_max'],
#                                                                                        TT_Sim['motor']['P_max'],
#                                                                                        TT_Sim['motor']['W_lim'], 50, 1)

# TT_Sim['motor']['w'] = np.array([0, 4166.7, 7333.3, 10500]) / 30 * np.pi  # Daley TTZ 2016 limits
# TT_Sim['motor']['t'] = np.array([106, 106, 65.791, 45.949])
# TT_Sim['motor']['p'] = TT_Sim['motor']['w'] * TT_Sim['motor']['t']
# SETTING LIKE THIS DOESN@T WORK AS THEY ARE RECALCULATED

# END OF USER PARAMETERS
TT_Sim['t'] = np.array([])
TT_Sim['v'] = np.array([])
TT_Sim['Iq'] = np.array([])
TT_Sim['torque'] = np.array([])
TT_Sim['Distance'] = np.array([])

TT_Sim['J']['linear'] = TT_Sim['constants']['m'] * TT_Sim['constants']['r'] ** 2
TT_Sim['J']['r'] = 2 * TT_Sim['J']['wheel'] + TT_Sim['J']['linear'] + np.square(TT_Sim['N'][0] / TT_Sim['N'][1]) * \
                                                                      TT_Sim['J']['motor']  # J referred to wheel

# Bit of a kludge to convert dict type to mat_structure
# THIS SHOULD BE ELIMINATED BY USING DICTS THROUGHOUT
sio.savemat('temp', {'TT_Sim': TT_Sim}, oned_as='column')
TT_Sim = sio.loadmat('temp', struct_as_record=False, squeeze_me=True)['TT_Sim']
remove('temp.mat')


timer = time.time()

Ref_Race = sio.loadmat(filename_ref_lap, struct_as_record=False, squeeze_me=True)[structure_lap]
v = Ref_Race.vGPS
# v = Ref_Race.v
d = Ref_Race.Distance

# Import course map and reference lap data
mat_contents = sio.loadmat(filename_ref_map, struct_as_record=False, squeeze_me=True)
Course_map = mat_contents[structure_map]
corners = sio.loadmat(filename_ref_brake, squeeze_me=True)  # Get track corner locations
locsmin = corners[var_name_brake] - 1  # -1 as matlab indexing starts at 1
# v = r * Ref_Race.Rpm / 30 * np.pi * N2 / N1
# v = v * 0.9
# del r
# del N2
# del N1  # Just to ensure they are not used accidentally
corner_index = np.arange(locsmin[c], locsmin[c + 1] - 1)

dt_a = 0.05  # 0.008
# t = np.linspace(Ref_Race.t[corner_index[0]], Ref_Race.t[corner_index[-1]], corner_index.size)
# [t, dt_a] = np.linspace(Ref_Race.t[corner_index[1]], Ref_Race.t[corner_index[0]], corner_index.size, retstep=True)
# print('dtA = ', str(dt_a))

#t_length = corner_index.size
t_length = 20/dt_a
t = Ref_Race.t[corner_index[0]] + np.linspace(0, t_length * dt_a - dt_a, t_length)
# print('dtA = ', str(t[1] - t[0]))

V = []          # list to hold solutions
D = []          # distances
G = []          # course gradients
A = []          # lean angle f(D)
lean = []
T_motor = []    # motor torque f(t)
J_l = []
J_r = []
R = []
H = []          # Course heading f(D)
V.append(v[corner_index[0]])    # Put y0 into solution
D.append(d[corner_index[0]])
G.append(np.interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
A.append(np.interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
lean.append(0.0)
R.append(TT_Sim.constants.r - 0.12 * (1 - np.cos(A[-1])))
H.append(np.interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
T_motor.append(Ref_Race.constants.Km * np.interp(t[0], Ref_Race.t, Ref_Race.Iq, -1, -1))
J_l.append(TT_Sim.constants.m * R[-1] ** 2)
J_r.append(2 * TT_Sim.J.wheel + J_l[-1] + np.square(TT_Sim.N[0] / TT_Sim.N[1]) * TT_Sim.J.motor)  # J referred to wheel

solver = ode(motorbike_mech4)
solver.set_integrator('lsoda', with_jacobian=False)
solver.set_initial_value(V[0], t[0])
solver.set_f_params(R[-1], TT_Sim.constants.rho, TT_Sim.constants.cd, J_r[-1], TT_Sim.constants.area,
                    TT_Sim.constants.m, TT_Sim.constants.p_tyre, T_motor[-1], TT_Sim.N[1], TT_Sim.N[0], G[-1])

print(t[0],V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
corner_index2=np.arange(locsmin[c], locsmin[c + 3] - 1)

for time in t[1:]:
    V.append(solver.integrate(time))
    G.append(np.interp(D[-1], Course_map.dist, Course_map.gradient, 0, 0))
    A.append(np.interp(D[-1], Course_map.dist, Course_map.lean, 0, 0))
    R.append(TT_Sim.constants.r - 0.12 * (1 - np.cos(A[-1])))
    T_motor.append(Ref_Race.constants.Km * np.interp(time, Ref_Race.t, Ref_Race.Iq, -1, -1))
    J_l.append(TT_Sim.constants.m * R[-1] ** 2)
    J_r.append(2 * TT_Sim.J.wheel + J_l[-1] + np.square(TT_Sim.N[0] / TT_Sim.N[1]) * TT_Sim.J.motor)  # J ref to wheel
    solver.set_f_params(R[-1], TT_Sim.constants.rho, TT_Sim.constants.cd, J_r[-1], TT_Sim.constants.area,
                        TT_Sim.constants.m, TT_Sim.constants.p_tyre, T_motor[-1], TT_Sim.N[1], TT_Sim.N[0], G[-1])
    dD = np.squeeze(V[-1] + V[-2]) * dt_a / 2.0
    D.append(D[-1] + dD)
    #dx = np.interp(D[-1], Course_map.dist, Course_map.x) - np.interp(D[-2], Course_map.dist, Course_map.x)
    #dy = np.interp(D[-1], Course_map.dist, Course_map.y) - np.interp(D[-2], Course_map.dist, Course_map.y)
    #H.append(np.arctan2(dx,dy))
    H.append(np.interp(D[-1], Course_map.dist, Course_map.heading, 0, 0))
    #if dH > np.pi / 2:
    #    dH -= np.pi
    #if dH < -np.pi / 2:
    #    dH += np.pi
    w = (H[-1] - H[-2]) / dt_a
    #w = np.interp(D[-1],Course_map.dist, Course_map.w)
    a_lateral = V[-1] * w

    #vvv = np.interp(D[-1],Ref_Race.Distance[corner_index2], Ref_Race.vGPS[corner_index2])
    #vvv = np.interp(D[-1],Course_map.dist,Course_map.dDist*20)
    #a_lateral = vvv * W[-1]
    lean.append(np.arctan(a_lateral / 9.81))
    #print(D[-1], H[-1], V[-1], lean[-1]*180/np.pi, A[-1]*180/np.pi)
    #print(time, V[-1], D[-1], G[-1], T_motor[-1], J_r[-1], A[-1], R[-1])
    # todo put field weakening calculation in here? include IR drop at least. Could just calc. motor curve for two adjecent points
    if not solver.successful():
        print('Warning: integration not successful')

#for dist in D[1:]:


V = np.squeeze(V)
R = np.squeeze(R)
T_motor = np.squeeze(T_motor)
A = np.squeeze(A)
lean = np.squeeze(lean)

MOTORSPEED = V / (1 / 60 * TT_Sim.N[1] / TT_Sim.N[0] * 2 * np.pi * R)  # in rpm

TT_Sim.v = V
TT_Sim.t = t
TT_Sim.Distance = D
TT_Sim.Iq = T_motor / TT_Sim.constants.Km

#print('Simulation duration =', str(time.time() - timer), 'seconds')
sio.savemat(filename_exp, {structure_exp: TT_Sim}, oned_as='column')
print('Simulation saved to', filename_exp, 'as structure named', structure_exp)

print('Motor max speed', str(max(TT_Sim.v) / TT_Sim.constants.r * 30 / np.pi * TT_Sim.N[0] / TT_Sim.N[1]))
print('Bike max speed (mph) = ', str(max(TT_Sim.v) * 2.23))
print('Simulated lap time (s) = ', str(TT_Sim.t[-1]))
print('Simulated lap speed (mph) = ', str(4.545*0.621/TT_Sim.t[-1]*3600))
# print(integrate.cumtrapz(TT_Sim.v, TT_Sim.t))




fig, ax1 = plt.subplots()
ax1.plot(TT_Sim.t, TT_Sim.Iq*TT_Sim.constants.Km, Ref_Race.t,Ref_Race.Iq*Ref_Race.constants.Km, 'b-')
ax1.set_xlabel('Time (s)')
ax1.set_xlim(TT_Sim.t[0], TT_Sim.t[-1])
ax1.set_ylim(0, 180)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Torque (Nm)', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(TT_Sim.t, TT_Sim.v, 'r-', Ref_Race.t, Ref_Race.vGPS, 'g-')
ax2.set_ylabel('Bike speed (m/s)', color='r')
ax2.set_xlim(TT_Sim.t[0], TT_Sim.t[-1])
ax2.set_ylim(0, 80)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
# plt.title('Motor Torque and Power vs Speed')
# plt.show()

#fig.show()
#plt.figure(2)
#plt.plot(Ref_Race.t, Ref_Race.Rpm, TT_Sim.t, MOTORSPEED, 'o')
#plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
#plt.show()

fig.show()
plt.figure(3)
plt.plot(Course_map.dist, Course_map.lean*180/np.pi, TT_Sim.Distance, lean*180/np.pi)
plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
plt.show()

#fig.show()
#plt.figure(3)
#plt.plot(TT_Sim.t, A*180/np.pi, TT_Sim.t, lean*180/np.pi)
#plt.xlim(TT_Sim.t[0], TT_Sim.t[-1])
#plt.show()

#fig.show()
#plt.figure(3)
#plt.plot(Course_map.dist, Course_map.heading, TT_Sim.Distance, H)
#plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
#plt.show()

#fig.show()
#plt.figure(3)
#plt.plot(Course_map.dist,Course_map.dDist*20, Ref_Race.Distance[corner_index2], Ref_Race.vGPS[corner_index2])
#plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
#plt.show()

max_V = abs(9.81*np.tan(49/180*np.pi)*dt_a/(Course_map.w/20))
fig.show()
plt.figure(4)
plt.plot(Course_map.dist, max_V, Ref_Race.Distance, Ref_Race.vGPS)
# plt.xlim(TT_Sim.Distance[0], TT_Sim.Distance[-1])
plt.ylim(0, max(Course_map.dDist*20))
plt.show()

plt.show()
