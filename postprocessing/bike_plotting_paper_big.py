import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import motorbike_functions as bike
from pylab import *

IR = 0.1
directory = '../data_export/'

Data_A = bike.energy_losses(directory + 'Python_sim_Portimao_A.mat', 'Portimao_sim_A', IR)
# Data_A = bike.energy_losses('Python_sim.mat', 'TT_sim', IR)
# Data_A = bike.energy_losses('Python_sim_PP.mat', 'PP_sim', IR)
Data_B = bike.energy_losses(directory + 'Python_sim_Portimao_DM1.mat', 'Portimao_sim_DM1', IR)
Data_C = bike.energy_losses(directory + 'Python_sim_Portimao_DM2.mat', 'Portimao_sim_DM2', IR)
Data_D = bike.energy_losses(directory + 'Python_sim_Portimao_D.mat', 'Portimao_sim_D', IR)

# TODO CHAIN LOSSES??

print(Data_A.Energy.Chain/3600,Data_B.Energy.Chain/3600,Data_C.Energy.Chain/3600,Data_D.Energy.Chain/3600)

labels = ['Air', 'Brake', 'Tyre', 'Battery', 'Motor', 'Drive', 'Chain']


def fracs(data):
    frac = [data.Energy.air, data.Energy.brake, data.Energy.roll, data.Energy.DC_loss,
         data.Energy.MotorLosses,  data.Energy.DriveLosses, data.Energy.Chain]
    print(data.Energy.Battery / 3600, sum(frac) / 3600)
    return frac


fig7 = plt.figure(7)
ax = fig7.add_subplot(2, 2, 1)
ax.pie(fracs(Data_A), labels=labels, autopct='%1.1f%%')
title('Data A Loss Distribution')
ax = fig7.add_subplot(2, 2, 2)
ax.pie(fracs(Data_B), labels=labels, autopct='%1.1f%%')
title('Data B Loss Distribution')
ax = fig7.add_subplot(2, 2, 3)
ax.pie(fracs(Data_C), labels=labels, autopct='%1.1f%%')
title('Data C Loss Distribution')
ax = fig7.add_subplot(2, 2, 4)
ax.pie(fracs(Data_D), labels=labels, autopct='%1.1f%%')
title('Data D Loss Distribution')
plt.show()
fig7.savefig('Python_pie_losses.pgf')


N = 7
# air = (Data_A.Energy.air, Data_B.Energy.air, Data_C.Energy.air, Data_D.Energy.air)
# brake = (Data_A.Energy.brake, Data_B.Energy.brake, Data_C.Energy.brake, Data_D.Energy.brake)
# tyre = (Data_A.Energy.roll, Data_B.Energy.roll, Data_C.Energy.roll, Data_D.Energy.roll)
# drive = (Data_A.Energy.DriveLosses, Data_B.Energy.DriveLosses, Data_C.Energy.DriveLosses, Data_D.Energy.DriveLosses)
# motor = (Data_A.Energy.MotorLosses, Data_B.Energy.MotorLosses, Data_C.Energy.MotorLosses, Data_D.Energy.MotorLosses)
# battery = (Data_A.Energy.DC_loss, Data_B.Energy.DC_loss, Data_C.Energy.DC_loss, Data_D.Energy.DC_loss)

ind = np.arange(N)  # the x locations for the groups
width = 1/(4+1)  # the width of the bars

fig, ax = plt.subplots()
# rects1 = ax.bar(ind, air, width, color='b')
# rects2 = ax.bar(ind + 1*width, brake, width, color='g')
# rects3 = ax.bar(ind + 2*width, tyre, width, color='r')
# rects4 = ax.bar(ind + 3*width, drive, width, color='c')
# rects5 = ax.bar(ind + 4*width, motor, width, color='m')
# rects6 = ax.bar(ind + 5*width, battery, width, color='y')
rects1 = ax.bar(ind + 0*width, [x / 1000 for x in fracs(Data_A)], width, color='b')
rects2 = ax.bar(ind + 1*width, [x / 1000 for x in fracs(Data_B)], width, color='g')
rects3 = ax.bar(ind + 2*width, [x / 1000 for x in fracs(Data_C)], width, color='r')
rects4 = ax.bar(ind + 3*width, [x / 1000 for x in fracs(Data_D)], width, color='c')

# add some text for labels, title and axes ticks
ax.set_ylabel('Energy (kJ)')
# ax.set_title('Scores by group and gender')
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(labels)

ax.legend((rects1[0], rects2[0], rects3[1], rects4[0]), ('A', 'B', 'C', 'D'))


def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                '%d' % int(height),
                ha='center', va='bottom')


#autolabel(rects1)
#autolabel(rects2)

plt.show()
fig.savefig('Python_bar_losses.pgf')


fig, ax1 = plt.subplots()
ax1.plot([], [], color='b', label='Air', linewidth=5)
ax1.plot([], [], color='g', label='Tyre', linewidth=5)
ax1.plot([], [], color='r', label='Battery', linewidth=5)
ax1.plot([], [], color='c', label='Motor', linewidth=5)
ax1.plot([], [], color='m', label='Drive', linewidth=5)
ax1.plot([], [], color='y', label='Chain', linewidth=5)
ax1.stackplot(Data_A.t, Data_A.P.air/1000, Data_A.P.roll/1000, Data_A.P.DC_loss/1000,
              Data_A.P.MotorLosses/1000, Data_A.P.DriveLosses/1000, Data_A.P.Chain_loss/1000)
ax1.set_xlim(60, Data_A.t[-1])
ax1.set_xlabel('Time (s)', fontsize=18)
ax1.set_ylabel('Power (kW)', fontsize=18)
ax1.legend(loc='upper left', fancybox=True, fontsize=16)
plt.setp(ax1.get_xticklabels(), fontsize=16)
plt.setp(ax1.get_yticklabels(), fontsize=16)

#plt.title('Motor Power vs Time')
plt.show()
#fig.savefig('Python_stack_losses.png')


fig, ax1 = plt.subplots()
ax1.plot([], [], color='b', label='Air', linewidth=5)
ax1.plot([], [], color='g', label='Tyre', linewidth=5)
ax1.plot([], [], color='r', label='Battery', linewidth=5)
ax1.plot([], [], color='c', label='Motor', linewidth=5)
ax1.plot([], [], color='m', label='Drive', linewidth=5)
ax1.plot([], [], color='y', label='Chain', linewidth=5)
P_total = Data_A.P.air + Data_A.P.roll + Data_A.P.DC_loss + Data_A.P.MotorLosses + \
          Data_A.P.DriveLosses + abs(Data_A.P.Chain_loss)
ax1.stackplot(Data_A.t, Data_A.P.air/P_total, Data_A.P.roll/P_total, Data_A.P.DC_loss/P_total,
              Data_A.P.MotorLosses/P_total, Data_A.P.DriveLosses/P_total, abs(Data_A.P.Chain_loss)/P_total)
ax1.set_xlim(0, Data_A.t[-1])
ax1.set_ylim(0, 1)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Power (%)')
ax1.legend(loc='lower right')

#plt.title('Motor Power vs Time')
plt.show()

fig.show()
