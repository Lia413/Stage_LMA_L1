from données import *
import numpy as np

# Méthode analytique
def solution_analytique(t):
    # Initialisation des énergies
    Ep = 0
    Ec = 0
    Em = 0

    if beta ** 2 > 1:
        # Cas 1: oscillations amorties
        r1 = omega_0 * (np.sqrt(beta ** 2 - 1) - beta)
        r2 = omega_0 * (-np.sqrt(beta ** 2 - 1) - beta)
        # Calcul des coefficients A et B
        A = (v_0 - r2 * u_0) / (r1 - r2)
        B = (u_0 * r1 - v_0) / (r1 - r2)
        # Calcul de u et v
        u = A * np.exp(r1 * t) + B * np.exp(r2 * t)
        v = r1 * A * np.exp(r1 * t) + r2 * B * np.exp(r2 * t)
        # Calcul des énergies
        Ec += 1 / 2 * m * v ** 2
        Ep += 1 / 2 * k * u ** 2
        Em += Ep + Ec
        return (u, v, Ec, Ep, Em)

    if beta ** 2 == 1:
        # Cas 2: oscillations critiques
        r = -beta * omega_0
        # Calcul des coefficients A et B
        A = v_0 + beta * omega_0 * u_0
        B = u_0
        # Calcul de u et v
        u = (A * t + B) * np.exp(r * t)
        v = A * np.exp(r * t) + A * t * r * np.exp(r * t) + r * B * np.exp(r * t)
        # Calcul des énergies
        Ec += 1 / 2 * m * v ** 2
        Ep += 1 / 2 * k * u ** 2
        Em += Ep + Ec
        return (u, v, Ec, Ep, Em)

    if beta ** 2 < 1:
        # Cas 3: oscillations non amorties
        mu = omega_0 * np.sqrt(1 - beta ** 2)
        # Calcul des coefficients A et B
        A = (v_0 + u_0 * beta * omega_0) / mu
        B = u_0
        # Calcul de u et v
        u = (A * np.sin(mu * t) + B * np.cos(mu * t)) * np.exp((-beta * omega_0) * t)
        v = (
            A * mu * np.cos(mu * t) * np.exp((-beta * omega_0) * t)
            + (-beta * omega_0) * A * np.sin(mu * t) * np.exp((-beta * omega_0) * t)
            - mu * B * np.sin(mu * t) * np.exp((-beta * omega_0) * t)
            + B * (-beta * omega_0) * np.cos(mu * t) * np.exp((-beta * omega_0) * t)
        )
        # Calcul des énergies
        Ec += 1 / 2 * m * v ** 2
        Ep += 1 / 2 * k * u ** 2
        Em += Ep + Ec
        return (u, v, Ec, Ep, Em)


def solution_analytique_force(t):
    # Solution analytique avec une force externe
    if omega_0 != omega_1:
        u = (1 / (omega_0 ** 2 - omega_1 ** 2)) * (
            np.cos(omega_1 * t) - np.cos(omega_0 * t)
        )
        return u / m
    else:
        u = (1 / 2) * t * np.sin(omega_0 * t) / omega_0
        return u / m


def F(t, f):
    # Fonction calculant la force F à un instant t donné
    F = np.array([[0], [(f * np.cos(omega_0 * t)) / m]])
    return F


def euler_explicite(force,Niter):
    t1 = t_0

    if force == "force":
        # Cas avec une force externe
        U = np.array([[u_0], [v_0]])
        A = Id + S * delta_temps
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = np.dot(A, U + (F(t1, f_0)) * delta_temps)
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])
            t1 += delta_temps
        return (liste_positions, liste_vitesses)

    else:
        # Cas sans force externe
        U = np.array([[u_0], [v_0]])
        A = Id + S * delta_temps
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = np.dot(A, U)
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])
        return (liste_positions, liste_vitesses)


def euler_implicite(force,Niter):
    t2 = t_0 + delta_temps

    if force == "force":
        # Cas avec une force externe
        U = np.array([[u_0], [v_0]])
        A = Id - S * delta_temps
        inv_A = np.linalg.inv(A)
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = np.dot(inv_A, U + (F(t2, f_0)) * delta_temps)
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])
            t2 += delta_temps
        return (liste_positions, liste_vitesses)

    else:
        # Cas sans force externe
        U = np.array([[u_0], [v_0]])
        A = Id - S * delta_temps
        inv_A = np.linalg.inv(A)
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = np.dot(inv_A, U)
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])
        return (liste_positions, liste_vitesses)


def trapèzes(force,Niter):
    t1 = t_0
    t2 = t_0 + delta_temps

    if force == "force":
        # Cas avec une force externe
        U = np.array([[u_0], [v_0]])
        P = np.dot(
            (2 * Id + S * delta_temps), np.linalg.inv((2 * Id - S * delta_temps))
        )
        Q = np.linalg.inv(2 * Id - S * delta_temps)
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = (
                P @ U + Q @ ((F(t2, f_0) + F(t1, f_0))) * delta_temps
            )  # "@"" est un raccourcis de "np.dot" et ici c'est utile.
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])
            t1 += delta_temps
            t2 += delta_temps

        return (liste_positions, liste_vitesses)
    else:
        # Cas sans force externe
        U = np.array([[u_0], [v_0]])
        P = np.dot(
            (2 * Id + S * delta_temps), np.linalg.inv((2 * Id - S * delta_temps))
        )
        liste_positions = [U[0, 0]]
        liste_vitesses = [U[1, 0]]

        for i in range(1, Niter):
            U = np.dot(P, U)
            liste_positions.append(U[0, 0])
            liste_vitesses.append(U[1, 0])

        return (liste_positions, liste_vitesses)


def energie_ressort(nom, bool, force,Niter):
    Ep_i = []
    Ec_i = []
    if nom == "trapèzes" and bool == False:
        liste_positions = trapèzes(force,Niter)[0]
        for i in range(len(liste_positions)):
            Ep_i.append(1 / 2 * k * (liste_positions[i]) ** 2)

    if nom == "trapèzes" and bool == True:
        liste_vitesses = trapèzes(force,Niter)[1]
        for i in range(len(liste_vitesses)):
            Ec_i.append(1 / 2 * m * (liste_vitesses[i]) ** 2)

    if nom == "euler_implicite" and bool == False:
        liste_positions = euler_implicite(force,Niter)[0]
        for i in range(len(liste_positions)):
            Ep_i.append(1 / 2 * k * (liste_positions[i]) ** 2)

    if nom == "euler_implicite" and bool == True:
        liste_vitesses = euler_implicite(force,Niter)[1]
        for i in range(len(liste_vitesses)):
            Ec_i.append(1 / 2 * m * (liste_vitesses[i]) ** 2)

    if nom == "euler_explicite" and bool == False:
        liste_positions = euler_explicite(force,Niter)[0]
        for i in range(len(liste_positions)):
            Ep_i.append(1 / 2 * k * (liste_positions[i]) ** 2)

    if nom == "euler_explicite" and bool == True:
        liste_vitesses = euler_explicite(force,Niter)[1]
        for i in range(len(liste_vitesses)):
            Ec_i.append(1 / 2 * m * (liste_vitesses[i]) ** 2)

    if nom == "analytique" and bool == False:
        liste_positions = solution_analytique(force,Niter)[0]
        for i in range(len(liste_positions)):
            Ep_i.append(1 / 2 * k * (liste_positions[i]) ** 2)

    if nom == "analytique" and bool == True:
        liste_vitesses = solution_analytique(force,Niter)[1]
        for i in range(len(liste_vitesses)):
            Ec_i.append(1 / 2 * m * (liste_vitesses[i]) ** 2)
    return (Ec_i, Ep_i)


def energie_mécanique(méthode, force,Niter):
    E = []
    v = méthode(force,Niter)[1]
    u = méthode(force,Niter)[0]
    for i in range(len(u)):
        E_i = 1 / 2 * (m * v[i] ** 2 + k * u[i] ** 2)
        E.append(E_i)
    return E


def erreur(nbr, methode, force):
    erreur_i_true=[]
    erreur_i=[]
    for i in range(nbr):
        erreur_i.append(
            methode(force,nbr)[0][i])
        erreur_i_true.append(solution_analytique(i * delta_temps)[0]
        )
    return np.log(np.square(np.subtract(erreur_i_true,erreur_i)).mean())


# Mettez un hashtag devant les lignes que vous ne voulez pas afficher

print("logN", np.log(Niter))


print("erreur e imp", erreur(Niter, euler_implicite, ""))
print("erreur trap", erreur(Niter, trapèzes, ""))
