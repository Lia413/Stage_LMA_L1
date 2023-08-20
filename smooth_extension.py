from données import *


def Mat_F(U):
    F = np.array([[0.0] for i in range(2 * nombre_de_ressorts)])
    for j in range(nombre_de_ressorts):
        F[j] = U[nombre_de_ressorts + j]
    F[nombre_de_ressorts] = -(
        k1 / m * (2 * U[0] - U[1])
        + k2 / m * (U[0] ** 2 - (U[1] - U[0]) ** 2)
        + k3 / m * (U[0] ** 3 + (U[0] - U[1]) ** 3)
        + 1 / (2 * m) * kse * (-2 * U[0] + U[1]) ** 2 / l ** 2
    )
    F[2 * nombre_de_ressorts - 1] = -(
        k1 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2])
        + k2 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2]) ** 2
        + k3 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2]) ** 3
        + 1
        / (2 * m)
        * kse
        * (
            U[nombre_de_ressorts - 1]
            - 2 * U[nombre_de_ressorts]
            + U[nombre_de_ressorts + 1]
        )
        ** 2
        / l ** 2
    )
    for j in range(1, nombre_de_ressorts - 1):
        F[nombre_de_ressorts + j] = -(
            k1 / m * (-U[j - 1] + 2 * U[j] - U[j + 1])
            + k2 / m * ((U[j] - U[j - 1]) ** 2 - (U[j] - U[j + 1]) ** 2)
            + k3 / m * ((U[j] - U[j - 1]) ** 3 + (U[j] - U[j + 1]) ** 3)
            + 1 / (2 * m) * kse * (U[j - 1] - 2 * U[j] + U[j + 1]) ** 2 / l ** 2
        )
    return F


def Jacobienne(U):
    M11 = np.zeros((nombre_de_ressorts, nombre_de_ressorts))
    M12 = np.eye((nombre_de_ressorts))
    M1 = np.concatenate((M11, M12), axis=1)
    M21 = np.zeros((nombre_de_ressorts, nombre_de_ressorts))

    M21[0, 0] = -(
        2 * k1 / m
        + 2 * k2 / m * ((U[0] - U[1]) + U[0])
        + 3 * k3 / m * ((U[0] - U[1]) ** 2 + U[0] ** 2)
        + 1 / (m) * kse * (-2 * U[0] + U[1]) / l ** 2
    )
    M21[0, 1] = (
        k1 / m
        + 2 * k2 / m * (U[1] - U[0])
        + 3 * k3 / m * (U[0] - U[1]) ** 2
        + 1
        / (m)
        * kse
        * (
            U[nombre_de_ressorts - 1]
            - 2 * U[nombre_de_ressorts]
            + U[nombre_de_ressorts + 1]
        )
        / l ** 2
    )

    M21[nombre_de_ressorts - 1, nombre_de_ressorts - 1] = -(
        k1 / m
        + 2 * k2 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2])
        + 3 * k3 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2]) ** 2
        + 1
        / (m)
        * kse
        * (
            U[nombre_de_ressorts - 1]
            - 2 * U[nombre_de_ressorts]
            + U[nombre_de_ressorts + 1]
        )
        ** 2
        / l ** 2
    )

    M21[nombre_de_ressorts - 1, nombre_de_ressorts - 2] = (
        k1 / m
        + 2 * k2 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2])
        + 3 * k3 / m * (U[nombre_de_ressorts - 1] - U[nombre_de_ressorts - 2]) ** 2
        + 1
        / (m)
        * kse
        * (
            U[nombre_de_ressorts - 1]
            - 2 * U[nombre_de_ressorts]
            + U[nombre_de_ressorts + 1]
        )
        / l ** 2
    )

    for i in range(1, nombre_de_ressorts - 1):

        M21[i, i - 1] = (
            k1 / m
            + 2 * k2 / m * (U[i] - U[i - 1])
            + 3 * k3 / m * (U[i] - U[i - 1]) ** 2
            + 1 / (m) * kse * (U[i - 1] - 2 * U[i] + U[i + 1]) / l ** 2
        )

        M21[i, i] = -(
            2 * k1 / m
            + 2 * k2 / m * (U[i + 1] - U[i - 1])
            + 3 * k3 / m * ((U[i] - U[i - 1]) ** 2 + (U[i] - U[i + 1]) ** 2)
            + 1 / (m) * kse * (U[i - 1] - 2 * U[i] + U[i + 1]) / l ** 2
        )

        M21[i, i + 1] = (
            k1 / m
            - 2 * k2 / m * (U[i] - U[i + 1])
            + 3 * k3 / m * (U[i] - U[i + 1]) ** 2
            + 1 / (m) * kse * (U[i - 1] - 2 * U[i] + U[i + 1]) / l ** 2
        )

    M2 = np.concatenate((M21, M11), axis=1)
    M = np.concatenate((M1, M2), axis=0)
    return M


def trap(N):
    U = U_0.copy()
    liste_positions = liste_u
    liste_vitesses = liste_v
    liste_F = []

    for i in range(1, N):
        F_n = Mat_F(U)
        liste_F.append(F_n)
        J_n = Jacobienne(U)
        U_n = U + delta_temps * np.dot(
            np.linalg.inv(Mat_Id - 0.5 * delta_temps * J_n), F_n
        )
        U = U_n
        for j in range(nombre_de_ressorts):
            position = U[j, 0]
            vitesse = U[j + nombre_de_ressorts, 0]
            liste_positions[i][j] = position
            liste_vitesses[i][j] = vitesse
    return liste_positions, liste_vitesses, liste_F


def energie(liste_pos, liste_vit):
    liste_Ep = np.zeros(nombre_de_points)
    liste_Ec = np.zeros(nombre_de_points)
    liste_Em = []

    for i in range(nombre_de_points):
        positions = np.array(liste_pos[i], dtype=object)
        vitesses = np.array(liste_vit[i], dtype=object)  # Convertir a un array numpy
        for j in range(nombre_de_ressorts):
            if j == 0:
                liste_Ep[i]=(                 
                    a1 * positions[j] ** 2
                    + a2 * positions[j] ** 3
                    + a3 * positions[j] ** 4                    
                )
                liste_Ep[i]=liste_Ep[i]+(
                    + kse / 2 * (-2 * positions[j] + positions[j+1]) ** 2 / l ** 2
                )
            else :
                liste_Ep[i]=liste_Ep[i]+(
                     a1 * (positions[j] - positions[j-1]) ** 2
                    + a2 * (positions[j] - positions[j-1]) ** 3
                    + a3 * (positions[j] - positions[j-1]) ** 4     
                )
                if j == nombre_de_ressorts or j+2 >= len(positions):
                    liste_Ep[i] = liste_Ep[i]
                else:
                    liste_Ep[i] = liste_Ep[i] + (
                        kse / 2 * (positions[j-1] - 2 * positions[j+1] + positions[j+2]) ** 2 / l ** 2)

                
            liste_Ec[i]=liste_Ec[i]+(
                1 / 2 * m * vitesses[j] ** 2
            )
        liste_Em.append(liste_Ec[i] + liste_Ep[i])
    return liste_Ep, liste_Ec, liste_Em


pos, vit, liste_F = trap(nombre_de_points)
Ep, Ec, Em = energie(pos, vit)
#print(pos)
print(Ep)
# print(liste_F)
# print(len(pos))
# print(len(Ep))
# print(len(Ec))
# print(len(Em))
