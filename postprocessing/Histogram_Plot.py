#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import motorbike_functions as bike
import scipy.io as sio

filename_exp = 'Python_TT_sim_UoN01_max.mat'
structure_exp = 'TT_sim_UoN01_max'

IR = 0.1
#Data_A = bike.energy_losses('Python_TT_sim_UoN01_max.mat', 'TT_sim_UoN01_max', IR)
#sio.savemat(filename_exp, {structure_exp: Data_A}, oned_as='column')
Data_A = bike.energy_losses('Python_TT_sim_UoN01_2015.mat', 'TT_sim_UoN01_2015', IR)
#sio.savemat(filename_exp, {structure_exp: Data_A}, oned_as='column')
# Data_A = bike.energy_losses('Python_sim_PP.mat', 'PP_sim', IR)
#Data_B = bike.energy_losses('Python_sim_Portimao_DM1.mat', 'Portimao_sim_DM1', IR)
Data_B = bike.energy_losses('Python_TT_sim_UoN01_120.mat', 'TT_sim_UoN01_120', IR)
Data_C = bike.energy_losses('Python_TT_sim_UoN01_max.mat', 'TT_sim_UoN01_max', IR)

mu, sigma = 100, 15
x = mu + sigma*np.random.randn(10000)

# the histogram of the data
n, bins, patches = plt.hist(Data_A.v, 50, normed=1, facecolor='green', alpha=0.75)
n, bins, patches = plt.hist(Data_B.v, 50, normed=1, facecolor='blue', alpha=0.75)
n, bins, patches = plt.hist(Data_C.v, 50, normed=1, facecolor='red', alpha=0.75)

# add a 'best fit' line
# y = mlab.normpdf( bins, mu, sigma)
# l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Speed (m/s)')
plt.ylabel('Probability')
# plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
plt.title('Histogram of Bike speed at TT vs Race')
plt.axis([0, 80, 0, 0.14])
plt.grid(True)

plt.show()