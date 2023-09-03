import numpy as np
from math import *


# données initiales


m = 0.1
g = 9.80
l = 0.449
l_0 = 0.400
k = m * g / (l - l_0)
print(k)


omega_0 = np.sqrt(k / m)
theta = 0.5
omega_1 = omega_0 * theta


T = 2 * pi / omega_0
t_0 = 0
t_fin = 10 * T
print("t_fin", t_fin)

nombre_de_points = 200
alpha = 0

beta = alpha / (2 * m * omega_0)
delta_temps = (t_fin - t_0) / nombre_de_points
print("delta_temps", delta_temps)
nombre_ressorts = 4


u_0 = 2
v_0 = -5
f_0 = 1  # Paramètre de force

# création de matrice utile
Mat_I2 = np.eye(2 * nombre_ressorts)
Mat_I = np.eye(nombre_ressorts)
Mat_M = np.dot(m, Mat_I)

Nulle = np.zeros((nombre_ressorts, nombre_ressorts))

periodique = True
k1 = 1
k2 = 1

Mat_A = np.zeros((nombre_ressorts, nombre_ressorts))
if periodique == True:
    Mat_A[0, 0] = 2 * k1
    Mat_A[0, 1] = -1 * k2
    Mat_A[nombre_ressorts - 1, nombre_ressorts - 1] = 2 * k1
    Mat_A[nombre_ressorts - 1, nombre_ressorts - 2] = -1 * k2
elif periodique == False:
    Mat_A[0, 0] = 2
    Mat_A[0, 1] = -1
    Mat_A[nombre_ressorts - 1, nombre_ressorts - 1] = 2
    Mat_A[nombre_ressorts - 1, nombre_ressorts - 2] = -1
for i in range(1, nombre_ressorts - 1):
    Mat_A[i, i - 1] = -1
    Mat_A[i, i] = 2
    Mat_A[i, i + 1] = -1
Mat_K = np.dot(Mat_A, k)


Mat_S_ligne1 = np.concatenate((Nulle, Mat_I), axis=1)
Mat_S_ligne2 = np.concatenate((-1 / m * Mat_K, -2 * beta * omega_0 * Mat_I), axis=1)
Mat_S = np.concatenate((Mat_S_ligne1, Mat_S_ligne2), axis=0)
# print("mat_S", Mat_S)



# on cherche à connaitre les valeurs propres et vecteurs propres de K pour la démonstration résolvant analytiquement le cas NDDL
inv_Mat_M = np.linalg.inv(m * Mat_I)


valp, vecp = np.linalg.eig(Mat_K)

# print("valp", valp)
# print("vecp", vecp)
omega_2 = np.sqrt(valp[0]) # on teste avec une valeur propre de la matrice cas pour voir si on obtient un cas limite
