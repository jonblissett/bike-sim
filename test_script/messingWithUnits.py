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
import matplotlib.pyplot as plt
import scipy.integrate as sci

import scipy.io as sio

from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

mat_contents = sio.loadmat('TT_Qualy_2016.mat', struct_as_record=False, squeeze_me=True)

TT_Qualy = mat_contents['TT_Qualy']



index = np.arange(1000,3000)

N = ([19, 67])
r = 2.001/np.pi/2 * ureg.meter      # wheel radius
time = Q_(TT_Qualy.t[index],'second')
# time = TT_Qualy.t * ureg.second
# time = time.astype(float)
w = TT_Qualy.Rpm[index]*N[0]/N[1]* ureg.rpm
v = w*r
v.ito_base_units()

plt.plot(TT_Qualy.t[index], v)
plt.xlabel('t')
plt.ylabel('v')
plt.show()

