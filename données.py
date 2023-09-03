import numpy as np
from math import *

# Données initiales
k = 2
m = 4
u_0 = 2
v_0 = 0
omega_0 = np.sqrt(k / m)
T = 2 * pi / omega_0
print("T", T)
t_0 = 0
t_fin = 10 * T
print("t_fin", t_fin)
nombre_de_points = 300
Niter=nombre_de_points
# Notation pour alléger le calcul et améliorer la compréhension
alpha = 0
beta = alpha / (2 * m * omega_0)

# Données nécessaires aux différentes méthodes
Delta = beta ** 2 * 4 * omega_0 ** 2 - 4 * omega_0 ** 2
S = np.array(
    [[0, 1], [-(omega_0 ** 2), -2 * beta * omega_0]]
)  # Matrice telle que F(U) = SU
Id = np.eye(2)
Id_3 = np.eye(3)
delta_temps = (t_fin - t_0) / nombre_de_points
print("delta_temps", delta_temps)


theta = 1
omega_1 = omega_0  # * theta
print(omega_0, "Omega0")
print(omega_1, "omega1")
f_0 = 1  # Paramètre de force

# Grille pour l'évolution de l'erreur
N = [10 * 2 ** i for i in range(1, 10)]


nddl = 20  # Taille de la matrice K
Mat_K = np.zeros((nddl, nddl))  # Initialisation de la matrice K

for i in range(nddl):
    Mat_K[i, i] = 2  # Diagonale 

    if i > 0:
        Mat_K[i, i - 1] = -1  # Termes à gauche de la diagonale 

    if i < nddl - 1:
        Mat_K[i, i + 1] = -1  # Termes à droite de la diagonale 
