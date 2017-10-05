# -*- coding: utf-8 -*-
# cython: profile=False
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
cimport numpy as np

ctypedef np.double_t DTYPE_t

# @cython.boundscheck(False) # turn off bounds checking for this func.
cpdef double fast_interp(double val, np.ndarray[DTYPE_t, ndim=1, negative_indices=False] x, np.ndarray[DTYPE_t, ndim=1, negative_indices=False] y,double left,double right):
    cdef unsigned int index = 0
    # cdef double _xrange, xdiff, modolo, ydiff
    cdef double y_interp = left
    while x[index] <= val:
        index += 1
        if index < len(x):
            # _xrange = x[index] - x[index-1]
            # xdiff = val - x[index-1]
            # modolo = xdiff / _xrange
            # ydiff = y[index] - y[index-1]
            # y_interp = y[index-1] + modolo*ydiff
            y_interp = y[index-1] + (val - x[index-1]) / (x[index] - x[index-1]) * (y[index] - y[index-1])
        else:
            return right
    return y_interp

"""
import numpy as np
from numpy.core.numerictypes import typecodes, number
from numpy.core.multiarray import digitize, bincount, interp as compiled_interp


def interp(x, xp, fp, left=None, right=None, period=None):
    #  One-dimensional linear interpolation.
    #  Cython implementation of  numpy.interp functionbase.py, v1.11.0
    # if period is None:
    #    if isinstance(x, (float, int, number)):
    #        return compiled_interp([x], xp, fp, left, right).item()
    #    elif isinstance(x, np.ndarray) and x.ndim == 0:
    #        return compiled_interp([x], xp, fp, left, right).item()
    #    else:
    #        return compiled_interp(x, xp, fp, left, right)
    # else:

    if period == 0:
        raise ValueError("period must be a non-zero value")
    period = abs(period)
    left = None
    right = None
    return_array = True
    if isinstance(x, (float, int, number)):
        return_array = False
        x = [x]
    x = np.asarray(x, dtype=np.float64)
    xp = np.asarray(xp, dtype=np.float64)
    fp = np.asarray(fp, dtype=np.float64)
    if xp.ndim != 1 or fp.ndim != 1:
        raise ValueError("Data points must be 1-D sequences")
    if xp.shape[0] != fp.shape[0]:
        raise ValueError("fp and xp are not of the same length")
    # normalizing periodic boundaries
    x = x % period
    xp = xp % period
    asort_xp = np.argsort(xp)
    xp = xp[asort_xp]
    fp = fp[asort_xp]
    xp = np.concatenate((xp[-1:]-period, xp, xp[0:1]+period))
    fp = np.concatenate((fp[-1:], fp, fp[0:1]))
    if return_array:
        return compiled_interp(x, xp, fp, left, right)
    else:
        return compiled_interp(x, xp, fp, left, right).item()
from bisect import bisect_left


def fast_interp(np.ndarray[np.double_t, ndim=1] val_array, np.ndarray[np.double_t, ndim=1] x, np.ndarray[np.double_t, ndim=1] y,left=None,right=None):
    cdef unsigned int index, trap
    cdef unsigned int ntraps=val_array.size
    cdef long double val, _xrange, xdiff, modolo, ydiff
    cdef np.ndarray[np.double_t, ndim=1] y_interp = np.zeros(ntraps, dtype=np.double_t)
    for trap in range(ntraps):
        index = 0
        val = val_array[trap]
    while x[index] <= val:
        index += 1
        _xrange = x[index] - x[index-1]
        xdiff = val - x[index-1]
        modolo = xdiff/_xrange
        ydiff = y[trap,index] - y[trap,index-1]
        y_interp[trap] = y[trap,index-1] + modolo*ydiff
    return y_interp


def fast_interp(val, x, y, left=None, right=None):
    index = bisect_left(x, val)
    _xrange = x[index] - x[index-1]
    xdiff = val - x[index-1]
    modolo = xdiff/_xrange
    ydiff = y[index] - y[index-1]
    return y[index-1] + modolo*ydiff

cpdef double fast_interp_iffy(double x, np.ndarray[np.double_t] xp, np.ndarray[np.double_t] fp,double left,double right):
    # if not extrapolate:
    #     x = np.clip(x, xp[0], xp[-1])
    # i = np.clip(np.searchsorted(xp, x), 1, len(xp) - 1)
    cdef int i # , j
    # cdef double xpj, fpj
    i = np.searchsorted(xp, x)
    if i < len(xp):
        if i > 0:
            # j = i - 1
            # xpj, fpj = xp[i-1], fp[j]
            # return   fpj + (x -     xpj) * (fp[i] -     fpj) / (xp[i] -     xpj)
            return fp[i-1] + (x - xp[i-1]) * (fp[i] - fp[i-1]) / (xp[i] - xp[i-1])
        else:
            return left
    else:
        return right


def _interp_linear(np.ndarray[np.double_t, ndim=1] yin not None,
                   np.ndarray[np.double_t, ndim=1] a not None,
                   np.ndarray[np.double_t, ndim=1] v not None):

    #yout = (xout - x0) / (x1 - x0) * (y1 - y0) + y0
    #Parameters
    #----------
    #yin : 1-D array_like
    #     Input y values
    #a : 1-D array_like
    #     Input x values corresponding to yin
    #v : 1-D array_like
    #     Output x values corresponding to yout

    cdef double vi
    cdef int nv = v.shape[0]
    cdef np.ndarray[np.double_t, ndim=1] yout = np.empty(nv, dtype=np.float64)
    cdef unsigned int na = a.shape[0]
    cdef unsigned int ia = 0
    cdef unsigned int iv
    cdef unsigned int na1 = na - 1
    cdef unsigned int na2 = na - 2
    cdef unsigned int ia1

    for iv in range(nv):
        vi = v[iv]
        while True:
            if ia < na:
                if vi <= a[ia]:
                    if ia == 0:
                        yout[iv] = (vi - a[0]) / (a[1] - a[0]) * (yin[1] - yin[0]) + yin[0]
                    else:
                        ia1 = ia - 1
                        yout[iv] = (vi - a[ia1]) / (a[ia] - a[ia1]) * (yin[ia] - yin[ia1]) + yin[ia1]
                    break
                else:
                    ia += 1
            else:  # ia == na without vi ever being less than a[ia]
                   # Thus vi > all values of a
                yout[iv] = (vi - a[na1]) / (a[na1] - a[na2]) * (yin[na1] - yin[na2]) + yin[na1]
                break

    return yout
"""


