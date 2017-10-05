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
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import motorbike_functions as bike
import scipy.io as sio
import seaborn as sns
import time
sns.set(style="ticks")

directory = './data_export/' # '../data_export/'

filename_exp = directory + 'Python_TT_sim_UoN01_max.mat'
structure_exp = directory + 'TT_sim_UoN01_max'

IR = 0.1
Data_A = bike.energy_losses(directory + 'Python_Sims2_best120s6p.mat', 'TT_Sim', IR)
v = Data_A.v
t = Data_A.t
Idc = Data_A.Idc

num = (t[-1] - t[0]) / max(np.diff(t))
if num / np.size(t) > 1:
    print('Data re-sampled to even timebase')
    t0 = np.linspace(t[0], t[-1], num=num)
    v = np.interp(t0, t, v)
    t = t0

import matplotlib.animation as animation

fig, ax = plt.subplots()

#x = np.arange(0, 2*np.pi, 0.01)        # x-array
#line, = ax.plot(x, np.sin(x))
#x = t
#x = t[range(0, 500)]
#line, = ax.plot(x, v[range(0, 500)])
#line, = ax.plot(t[0],v[0])
# ax.plot(t, v, t, Idc*0.1)
line, = ax.plot([], [])
#plt.xlim(t[0],t[-1])


def animate(i):
    # line.set_ydata(np.sin(x+i/10.0))  # update the data
    #line.set_xdata(t[np.arange(0, i)])
    #line.set_ydata(v[range(0, i)])

    # Try plotting just the data in the interval and the x limits - maybe faster?
    if i % 5 == 0:
        plt.xlim(t[i]-10, t[i])
        print(t[i])
    if i == 999:
        print('Simulation duration =', str(time.time() - timer), 'seconds')
    return line,

#Init only required for blitting to give a clean slate.
def init():
    ax.set_ylim(0, max(v) * 1.1)
    ax.set_xlim(0, 10)
    ax.plot(t, v, t, Idc * 0.1)
    #line.set_ydata(np.ma.array(x, mask=True))
    return line,

timer = time.time()
ani = animation.FuncAnimation(fig, animate, np.arange(1, 1000),
    interval=1, repeat=False, init_func=init)#, blit=True)
plt.show()