from motorbike_functions import v_dq_pmsm

Iq = 600
Id = 0
N = 8.5
c = 175
w = 10500/30*3.14159

pp = 12
Ke = 0.34688 / (pp/2) / 11.5 * N * c/150
R = 0.007313 * N / 11.5  # Windings in parallel
Ld = 53e-6 * (N*N)/(11.5*11.5) * c/150
Lq = 61e-6 * (N*N)/(11.5*11.5) * c/150

[v_s, v_d, v_q, power_factor] = v_dq_pmsm(Ke, pp, R, Ld, Lq, 0, Iq, w)

print('R =', R)
print('Ke =', Ke)
print('Ld, Lq = ', Ld, Lq)
print('vd, vq = ', v_d, v_q)
print('vs =', v_s)
print('w =', w)
