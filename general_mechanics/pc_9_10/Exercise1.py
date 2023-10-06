import matplotlib.pyplot as plt
import numpy as np


a = 1
nu = 0.3
p = 1000

N_zr = lambda r: -p*r/2

# case: simply supported

Mrr_ss = lambda r: p/16*(3+nu)*(a**2 - r**2)
Mtt_ss = lambda r: p/16*((1+3*nu)*r**2 - (3+nu)*a**2)

# case: clamped
Mrr_c = lambda r: p/16*(a**2*(1+nu) - r**2*(3+nu))
Mtt_c = lambda r: -p/16*(a**2*(1+nu) - r**2*(1+3*nu))


r_vec = np.linspace(0, a, 100)

plt.figure()
plt.plot(r_vec, N_zr(r_vec), label="Nzr")
plt.plot(r_vec, Mrr_ss(r_vec), label="Mrr simply supported")
plt.plot(r_vec, Mtt_ss(r_vec), label="Mtt simply supported")

plt.plot(r_vec, Mrr_c(r_vec), label="Mrr clamped")
plt.plot(r_vec, Mtt_c(r_vec), label="Mtt clamped")

plt.grid()
plt.xlabel("r")
plt.legend()
plt.show()
