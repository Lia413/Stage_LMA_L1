import numpy as np
from math import *


# Données initiales
k1 = 2 * (59 / 24)
k2 = 3 * (-7 / 6)
k3 = 4 * (7 / 48)

Delta = k2 ** 2 - 4 * k1 * k3
print("Delta:", Delta)


u_stable = ((-k2) + np.sqrt(Delta)) / (2 * k3)
u_sup_stable = abs(u_stable) + 5
u_inf_stable = abs(u_stable) - 5

u_instable = ((-k2) - np.sqrt(Delta)) / (2 * k3)
u_sup_instable = abs(u_instable) + 5
u_inf_instable = abs(u_instable) - 5


print("u1:", u_stable)
print("u2:", u_instable)

m = 1.0
u_0 = u_instable * 0.9  # perturber aussi avec 1.01
v_0 = 0
omega_0 = np.sqrt(k1 / m)
T = 2 * pi / omega_0
eps = 0.01

t_0 = 0
t_fin = 5 * T
print("t_fin:", t_fin)
nombre_de_points = 1000
delta_temps = (t_fin - t_0) / nombre_de_points
print("delta_temps:", delta_temps)

Mat_Id = np.eye(2)

# Notation pour alléger le calcul et améliorer la compréhension
alpha = 0
beta = alpha / (2 * m * omega_0)

S = np.array(
    [[0, 1], [-(omega_0 ** 2), -2 * beta * omega_0]]
)  # Matrice telle que F(U) = SU
