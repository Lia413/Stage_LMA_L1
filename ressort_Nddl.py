from données import *


def F(t, f, i):
    # Fonction calculant la force F à un instant t donné
    force = np.zeros(2 * nombre_ressorts)
    force[i] = (f * np.cos(omega_2 * t)) / m
    return force


# Affichage de la force pour un instant donné
# print(F(delta_temps, f_0, 1))


def energie(W):
    U = W[0:nombre_ressorts]
    V = W[nombre_ressorts : 2 * nombre_ressorts]
    Ec = 0.5 * V @ Mat_M @ V
    Ep = 0.5 * U @ Mat_K @ U
    Em = Ec + Ep
    return (Ec, Ep, Em)


def ressort_NDDL(methode, N, force, i):
    # Initialisation des conditions initiales                      # N = nb de points
    U = np.zeros(2 * nombre_ressorts)
    U[0] = u_0
    U[nombre_ressorts] = v_0
    t1 = t_0
    t2 = t_0 + delta_temps

    # Initialisation des listes pour stocker les positions et les vitesses
    liste_positions = [[u_0]] + [[0] for k in range(1, nombre_ressorts)]
    liste_vitesses = [[v_0]] + [[0] for k in range(1, nombre_ressorts)]

    liste_Ec = []
    liste_Ep = []
    liste_Em = []
    (Ec, Ep, Em) = energie(U)

    A = 2 * Mat_I2 + Mat_S * delta_temps
    B = np.linalg.inv(2 * Mat_I2 - Mat_S * delta_temps)
    C = Mat_I2 - Mat_S * delta_temps
    Q = np.linalg.inv(2 * Mat_I2 - Mat_S * delta_temps)
    inv_C = np.linalg.inv(Mat_I2 - Mat_S * delta_temps)
    M0 = B @ A

    if methode == "trapèzes" and not (force == "force"):
        for k in range(1, N + 1):
            U = M0 @ U
            (Ec, Ep, Em) = energie(U)
            liste_Ec.append(Ec)
            liste_Ep.append(Ep)
            liste_Em.append(Em)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])
        return (liste_positions, liste_vitesses, liste_Ec, liste_Ep, liste_Em)

    elif methode == "trapèzes" and force == "force":
        # Méthode des trapèzes avec force externe
        for k in range(1, N):
            F_t2 = F(t2, f_0, i)
            F_t1 = F(t1, f_0, i)
            U = np.dot(M0, U) + Q @ (F_t2 + F_t1) * delta_temps
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])
            t1 += delta_temps
            t2 += delta_temps
        return (liste_positions, liste_vitesses)

    elif methode == "euler_implicite" and not (force == "force"):
        for k in range(1, N + 1):
            U = C @ U
            (Ec, Ep, Em) = energie(U)
            liste_Ec.append(Ec)
            liste_Ep.append(Ep)
            liste_Em.append(Em)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])

        return (liste_positions, liste_vitesses, liste_Ec, liste_Ep, liste_Em)

    elif methode == "euler_implicite" and force == "force":
        # Méthode d'Euler implicite avec force externe
        for k in range(1, N):
            U = np.dot(inv_C, U + F(t2, f_0, i) * delta_temps)
            for j in range(nombre_ressorts):
                liste_positions[j].append(U[j])
                liste_vitesses[j].append(U[j + nombre_ressorts])
            t2 += delta_temps
        return (liste_positions, liste_vitesses)


# Résultats
resultat = ressort_NDDL("trapèzes", nombre_de_points, "", 0)
liste_positions = resultat[0]
liste_vitesses = resultat[1]


def solution_analytique_Nddl(t):
    V = []
    for i in range(nombre_ressorts):
        # Calcul de la solution analytique pour chaque ressort
        X = vecp[i] * np.cos(np.sqrt(abs(valp[i])) * t) + vecp[i] * np.sin(
            np.sqrt(abs(valp[i])) * t
        )
        V.append(X)
    return np.array(V)


# Affichage de la solution analytique au temps t_0
# print("V", solution_analytique_Nddl(t_0))


def erreur_quadratique_trapèzes(N):
    erreur_quad = 0.0
    positions_trapèzes = ressort_NDDL("trapèzes", N, "", 0)[0]
    for i in range(N):
        t = i * delta_temps
        positions_analytique = solution_analytique_Nddl(t)
        for j in range(nombre_ressorts):
            # Calcul de l'erreur quadratique entre les positions calculées et les positions analytiques
            diff = positions_trapèzes[j][i] - positions_analytique[j]
            erreur_quad += diff ** 2

    erreur_quad /= N * nombre_ressorts
    return erreur_quad


# Calcul et affichage de l'erreur quadratique pour la méthode des trapèzes
# print("erreur", erreur_quadratique_trapèzes(nombre_de_points))
