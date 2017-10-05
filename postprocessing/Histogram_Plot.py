#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import motorbike_functions as bike
import seaborn as sns
import scipy.integrate as int
sns.set(style="ticks")

directory = '../data_export/'

filename_exp = directory + 'mymotor_300Nm'  # 'Python_TT_sim_UoN01_max.mat'
structure_exp = directory + 'TT_Sim'   # 'TT_sim_UoN01_max'

IR = 0.1
# Data_A = bike.energy_losses(directory + 'Python_Sims2_best120s6p.mat', 'TT_Sim', IR)
# Data_A = bike.energy_losses(directory + 'Portimao_FW_10C_150_4lap.mat', 'TT_Sim', IR)
Data_A = bike.energy_losses(directory + 'Python_Sim_import_motor_power_varied_mph_400Nm.mat', 'TT_Sim', IR)

print('Motor Losses (Wh): %d' % (Data_A.Energy.MotorLosses/3600))
print('Drive Losses (Wh): %d' % (Data_A.Energy.DriveLosses/3600))
print('Chain Losses (Wh): %d' % (Data_A.Energy.Chain/3600))
print('Total        (Wh): %d' % ((Data_A.Energy.MotorLosses+Data_A.Energy.DriveLosses+Data_A.Energy.Chain)/3600))
print('Air   Losses (Wh): %d' % (Data_A.Energy.air/3600))
print('Peak Motor loss (W): %d' % max(Data_A.P.MotorLosses))
print('Peak Drive loss (W): %d' % max(Data_A.P.DriveLosses))
print('Peak Chain loss (W): %d' % max(Data_A.P.Chain_loss))
print('Peak Aero  loss (W): %d' % max(Data_A.P.air))


# sio.savemat(filename_exp, {structure_exp: Data_A}, oned_as='column')
# Data_A = bike.energy_losses(directory + 'Python_TT_sim_UoN01_2015.mat', 'TT_sim_UoN01_2015', IR)
# Data_A = bike.energy_losses(directory + 'Python_sim_PP.mat', 'PP_sim', IR)
# sio.savemat(filename_exp, {structure_exp: Data_A}, oned_as='column')
# Data_A = bike.energy_losses('Python_sim_PP.mat', 'PP_sim', IR)
Data_B = bike.energy_losses(directory + 'Python_sim_Portimao_DM1.mat', 'Portimao_sim_DM1', IR)
# Data_B = bike.energy_losses('Python_TT_sim_UoN01_120.mat', 'TT_sim_UoN01_120', IR)
Data_C = bike.energy_losses(directory + 'Python_TT_sim_UoN01_2015.mat', 'TT_sim_UoN01_2015', IR)

# the histogram of the data
n, bins, patches = plt.hist(Data_A.v, 50, normed=1, facecolor='green', alpha=0.75)
n, bins, patches = plt.hist(Data_B.v, 50, normed=1, facecolor='blue', alpha=0.75)
n, bins, patches = plt.hist(Data_C.v, 50, normed=1, facecolor='red', alpha=0.75)

# add a 'best fit' line
# y = mlab.normpdf( bins, mu, sigma)
# l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Speed (m/s)', size=16)
plt.ylabel('Probability', size=16)
# plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
# plt.title('Histogram of Bike speed in two races')
# plt.axis([0, 80, 0, 0.14])
plt.grid(True)

plt.show()


t = Data_A.t
x = Data_A.Rpm
y = bike.motor_torque(Data_A.motor.co, Data_A.Iq)
# Z = int.cumtrapz(Data_A.P.DriveLosses, t, initial=0) + int.cumtrapz(Data_A.P.MotorLosses, t, initial=0)
# Z = int.trapz(Data_A.P.MotorLosses, t)
# Z = Data_A.P.DriveLosses  # int.cumtrapz(Data_A.P.air, t, initial=0)
Z = Data_A.P.MotorLosses


num = (t[-1] - t[0]) / min(np.diff(t))
if num / np.size(t) > 1:
    print('Data re-sampled to even timebase')
    t0 = np.linspace(t[0], t[-1], num=num)
    x = np.interp(t0, t, x)
    y = np.interp(t0, t, y)
    Z = np.interp(t0, t, Z)
    t = t0


x2 = x[y < 0]
y2 = abs(y[y < 0])
if x2.size == 0:
    regen = False
else:
    regen = True
x1 = x[y > 0]
y1 = abs(y[y > 0])

#  ax = sns.jointplot(x, abs(y), gridsize=50, kind="hex", space=0, stat_func=None, color="#4CB391", xlim=(0, 1.05*max(x)), ylim=(1.05 * min(y), 1.05 * max(y))).set_axis_labels('Motor Speed (RPM)', 'Motor Torque (Nm)', size=16).plot_joint(sns.kdeplot, zorder=1, n_levels=10, gridsize=30)

# plt.show()
# ax = sns.jointplot(x1, abs(y1), C=Z, gridsize=50, kind="kde", space=0, stat_func=None, xlim=(0, 1.05*max(x)), ylim=(1.05 * min(y), 1.05 * max(y)))
# plt.show()


g = sns.JointGrid(x1, abs(y1), space=0)
# g.plot_joint(plt.plot, linewidth=0.5, linestyle='dotted', color="Green")
# g.x = x1
# g.y = y1
g = g.plot_marginals(sns.distplot, kde=False, color="Blue")
g = g.plot_joint(sns.kdeplot, gridsize=50, cmap="Blues", shade=True, shade_lowest=False, cut=0, figsize=(8, 8))
g.set_axis_labels("Motor speed (Rpm)", "Torque (Nm)", size=16)
if regen:
    g.x = np.append(x2, np.linspace(max(x2), max(x1), 50))
    g.y = np.append(y2, np.linspace(min(y2), max(y1), 50))
    g = g.plot_marginals(sns.distplot, kde=False, color="Red")
    g.plot_joint(sns.kdeplot, cmap='Reds', gridsize=200, shade=False, shadelowest=False, cut=0)
else:
    g.x = np.append(x2, np.linspace(0, max(x1), 50))
    g.y = np.append(y2, np.linspace(0, max(y1), 50))

plt.show()

# Set up the figure
f, ax = plt.subplots(figsize=(8, 8))

# Draw the two density plots
if regen:
    ax = sns.kdeplot(x2, y2, gridsize=50, cmap="Reds", shade=True, shade_lowest=False, cut=0)  #, xlim=(0, 1.05*max(x)), ylim=(1.05 * min(y), 1.05 * max(y)))
ax = sns.kdeplot(x1, y1, gridsize=25, cmap="Blues", shade=True, shade_lowest=False, cut=0)  #, xlim=(0, 1.05*max(x)), ylim=(1.05 * min(y), 1.05 * max(y)))

# Add labels to the plot
red = sns.color_palette("Reds")[-2]
blue = sns.color_palette("Blues")[-2]
ax.text(6000, 200, "Acceleration", size=16, color=blue)
if regen:
    ax.text(2000, 50, "Braking", size=16, color=red)
else:
    ax.text(2000, 50, "No Braking", size=16, color=red)
plt.xlabel('Motor speed (Rpm)', size=16)
plt.ylabel('Torque (Nm)', size=16)

plt.show()

# plot
# ========================================
graph = sns.jointplot(x=x1, y=y1, kind='kde', color='b', gridsize=50, space=0, shade_lowest=False)
if regen:
    graph.x = x2
    graph.y = y2
    graph.plot_joint(sns.kdeplot, cmap='Reds', gridsize=50, shade_lowest=False)
    plt.xlabel('Motor speed (Rpm)', size=16)
    plt.ylabel('Torque (Nm)', size=16)
else:
    graph.set_axis_labels("Motor speed (Rpm)", "Torque (Nm)", size=16)

# graph.text(6000, 200, "Acceleration", size=16, color=blue)
# graph.text(2000, 50, "Braking", size=16, color=red)

plt.show()
