from données import *


def F(t, f, i):
    # Fonction calculant la force F à un instant t donné
    force = np.zeros(2 * nombre_ressorts)
    force[i] = (f * np.cos(omega_2 * t)) / m
    return force


# Affichage de la force pour un instant donné
print(F(delta_temps, f_0, 1))


def ressort_NDDL(methode, N, force, i):
    t1 = t_0
    t2 = t_0 + delta_temps

    U = np.zeros(2 * nombre_ressorts)
    U[0] = u_0
    U[nombre_ressorts] = v_0

    liste_positions = [[u_0] for _ in range(nombre_ressorts)]
    liste_vitesses = [[v_0] for _ in range(nombre_ressorts)]

    A = 2 * Mat_I2 + Mat_S * delta_temps
    B = np.linalg.inv(2 * Mat_I2 - Mat_S * delta_temps)
    C = Mat_I2 - Mat_S * delta_temps
    Q = np.linalg.inv(2 * Mat_I2 - Mat_S * delta_temps)
    inv_C = np.linalg.inv(Mat_I2 - Mat_S * delta_temps)

    if methode == "trapèzes" and force == "force":
        # Méthode des trapèzes avec force externe
        M0 = np.dot(A, B)
        for _ in range(1, N):
            F_t2 = F(t2, f_0, i)
            F_t1 = F(t1, f_0, i)
            U = np.dot(M0, U) + Q @ (F_t2 + F_t1) * delta_temps
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])
            t1 += delta_temps
            t2 += delta_temps
    elif methode == "trapèzes" and not (force == "force"):
        # Méthode des trapèzes sans force externe
        M0 = np.dot(A, B)
        for _ in range(1, N):
            U = np.dot(M0, U)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])

    if methode == "euler_implicite" and force == "force":
        # Méthode d'Euler implicite avec force externe
        for _ in range(1, N):
            U = np.dot(inv_C, U + F(t2, f_0, i) * delta_temps)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])
            t2 += delta_temps
    elif methode == "euler_implicite" and not (force == "force"):
        # Méthode d'Euler implicite sans force externe
        for _ in range(1, N):
            U = np.dot(inv_C, U)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])

    return (liste_positions, liste_vitesses)


# Affichage de la force pour un instant donné
print(F(delta_temps, f_0, 1))

# Affichage des positions pour la méthode des trapèzes avec force externe pour le 5ème ressort
print(ressort_NDDL("trapèzes", nombre_de_points, np.array("force"), 1)[0])
