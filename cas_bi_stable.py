import numpy as np
from données import *

# Fonctions utiles
def Mat_U():
    U = np.array([[u_0], [v_0]])
    return U


def Mat_F(U):
    return np.array(
        [
            [U[1, 0]],
            [-k1 / m * U[0, 0] - k2 / m * U[0, 0] ** 2 - k3 / m * U[0, 0] ** 3],
        ]
    )


def Jacobienne(U):
    return np.array(
        [[0, 1], [-k1 / m - 2 * k2 / m * U[0, 0] - 3 * k3 / m * U[0, 0] ** 2, 0,],]
    )


##############################################################################################

#                                Jacobienne et Newton Raphson                               #

##############################################################################################


def J_trap(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    for i in range(1, N):
        F_n = Mat_F(U)
        J_n = Jacobienne(U)
        U_n = U + delta_temps * np.dot(
            np.linalg.inv(Mat_Id - 0.5 * delta_temps * J_n), F_n
        )
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def J_trap_lim(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    A = Mat_Id - 0.5 * delta_temps * S
    inv_A = np.linalg.inv(A)
    B = Mat_Id + 0.5 * delta_temps * S

    for i in range(1, N):
        U_n = inv_A @ B @ U
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def J_trap_lim_b(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    A = Mat_Id - 0.5 * delta_temps * S
    inv_A = np.linalg.inv(A)

    for i in range(1, N):
        U_n = U + delta_temps * inv_A @ S @ U
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def J_euler_implicite(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    for i in range(1, N):
        F_n = Mat_F(U)
        J_n = Jacobienne(U)
        U_n = U + delta_temps * np.dot(np.linalg.inv(Mat_Id - delta_temps * J_n), F_n)
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def J_euler_implicite_lim(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    A = Mat_Id - delta_temps * S

    for i in range(1, N):
        U_n = np.linalg.inv(A) @ U
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def J_euler_implicite_lim_b(N):
    U = Mat_U()
    liste_positions = [U[0, 0]]
    liste_vitesses = [U[1, 0]]

    A = Mat_Id - delta_temps * S

    for i in range(1, N):
        U_n = U + delta_temps * np.linalg.inv(A) @ S @ U
        U = U_n
        liste_positions.append(U[0, 0])
        liste_vitesses.append(U[1, 0])

    return liste_positions, liste_vitesses


def energie(nombre_de_points, methode):
    E = []
    Ec = []
    Ep = []

    if methode == "euler_implicite":
        positions = J_euler_implicite(nombre_de_points)[0]
        vitesses = J_euler_implicite(nombre_de_points)[1]
    elif methode == "trapèzes":
        positions = J_trap(nombre_de_points)[0]
        vitesses = J_trap(nombre_de_points)[1]

    for i in range(nombre_de_points):
        Ec_i = 0.5 * m * vitesses[i] ** 2
        Ec.append(Ec_i)

        Ep_i = (
            0.5 * k1 * positions[i] ** 2
            + 1 / 3 * k2 * positions[i] ** 3
            + 1 / 4 * k3 * positions[i] ** 4
        )
        Ep.append(Ep_i)

        E_i = Ec_i + Ep_i
        E.append(E_i)

    return E, Ec, Ep


#####################################################################################

#               Calcul des listes de positions, vitesses, et énergies               #

#####################################################################################


positions_trap, vitesses_trap = J_trap(nombre_de_points)
positions_euler, vitesses_euler = J_euler_implicite(nombre_de_points)

positions_trap_lim, vitesses_trap_lim = J_trap_lim(nombre_de_points)
positions_euler_lim, vitesses_euler_lim = J_euler_implicite_lim(nombre_de_points)

positions_trap_lim_b, vitesses_trap_lim_b = J_trap_lim_b(nombre_de_points)
positions_euler_lim_b, vitesses_euler_lim_b = J_euler_implicite_lim_b(nombre_de_points)

# Calcul des énergies
Em_tr, Ec_tr, Ep_tr = energie(nombre_de_points, "trapèzes")
Em_euler, Ec_euler, Ep_euler = energie(nombre_de_points, "euler_implicite")

# Affichage des résultats
print("Positions trap: ", positions_trap)
print("Vitesses trap: ", vitesses_trap)

print("Positions euler: ", positions_euler)
print("Vitesses euler: ", vitesses_euler)
