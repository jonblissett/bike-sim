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

