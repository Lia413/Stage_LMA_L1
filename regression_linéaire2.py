from données import *
from ressort_1DDL_lineaire3 import *
import matplotlib.pyplot as plt
import numpy as np

log_Nombre = erreur_e_imp[2]


E_euler_imp = erreur_e_imp[0]
E_trap = erreur_trap[0]


def SS(sigma, n):
    s = 0
    for i in range(len(n)):
        s += 1 / (sigma ** 2)
    return s


def Sxx(n):
    S_x_x = 0
    for elmt in n:
        S_x_x += elmt ** 2
    return S_x_x


def Sx(n):
    S_x = 0
    for elmt in n:
        S_x += elmt
    return S_x


def Sy(E):
    S_y = 0
    for elmt in E:
        S_y += elmt
    return S_y


def Sxy(Sx, Sy):
    S_x_y = 0
    for i in range(len(Sx)):
        S_x_y += Sx[i] * Sy[i]
    return S_x_y


def Delta(S, Sx, S_x_x):
    return S * S_x_x - (Sx) ** 2


def a(Delta, Sx, Sy, Sxy, Sxx):
    return (Sxx * Sy - Sx * Sxy) / Delta


def b(Delta, S, Sxy, Sx, Sy):
    return (S * Sxy - Sx * Sy) / Delta


def coeff(Nombre, E):
    s = SS(1, Nombre)
    S_x = Sx(Nombre)
    S_x_x = Sxx(Nombre)
    S_y = Sy(E)
    S_x_y = Sxy(Nombre, E)
    delta = Delta(s, S_x, S_x_x)
    a_coeff = a(delta, S_x, S_y, S_x_y, S_x_x)
    b_coeff = b(delta, s, S_x_y, S_x, S_y)
    return (a_coeff, b_coeff)


# Les coefficients pour chaque méthode

euler_imp_a_coeff = coeff(log_Nombre, E_euler_imp)[0]
euler_imp_b_coeff = coeff(log_Nombre, E_euler_imp)[1]
trap_a_coeff = coeff(log_Nombre, E_trap)[0]
trap_b_coeff = coeff(log_Nombre, E_trap)[1]

# Test d'efficacité
def Chi_2(a, b, sigma, N, E):
    chi_2 = 0
    for i in range(len(N)):
        chi_2 += ((E[i] - a - b * N[i]) / sigma) ** 2
    return chi_2


print("chi2", Chi_2(euler_imp_a_coeff, euler_imp_b_coeff, 1, log_Nombre, E_euler_imp))

# Régression linéaire
def regression_lineaire(E, Nombre, methode, a_coeff, b_coeff):
    Nombre = np.array(Nombre)
    E = np.array(E)
    plt.plot(Nombre, E, "yo", label="valeurs")
    plt.plot(
        Nombre,
        b_coeff * Nombre + a_coeff,
        "--k",
        label=f"{round(b_coeff, 2)}x + {round(a_coeff, 2)}",
    )  # tronque les coefficient à 2 chiffres après la virgule
    plt.xlabel("log du nombre de points")
    plt.ylabel("log des erreurs")
    if methode == "euler_implicite":
        plt.title(
            "Erreurs en fonction du nombre de points avec la méthode d'Euler implicite"
        )
    if methode == "trapezes":
        plt.title(
            "Erreurs en fonction du nombre de points avec la méthode des trapèzes"
        )
    plt.legend()
    plt.show()
    plt.close()


# Affichage des graphiques de régression linéaire

regression_lineaire(
    E_euler_imp, log_Nombre, "euler_implicite", euler_imp_a_coeff, euler_imp_b_coeff
)
regression_lineaire(E_trap, log_Nombre, "trapezes", trap_a_coeff, trap_b_coeff)
