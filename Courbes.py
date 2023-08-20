from données import *
from ressort_1DDL_lineaire3 import *
import matplotlib.pyplot as plt

temps = np.linspace(t_0, t_fin, nombre_de_points)


def courbes_comparaison(force):
    if force != "force":
        plt.plot(
            temps, solution_analytique(temps)[0], "r", label="Solution analytique",
        )
    elif force == "force":
        plt.plot(
            temps,
            solution_analytique_force(temps),
            "r",
            label="Solution analytique avec force",
        )
    plt.plot(temps, trapèzes(force,Niter)[0], "k--", label="Trapèzes")
    if force != "force":
        plt.plot(temps, euler_explicite(force,Niter)[0], "c*", label="Euler explicite")
    plt.plot(temps, euler_implicite(force,Niter)[0], "g|", label="Euler implicite")
    plt.xlabel("Temps")
    plt.ylabel("Déplacement")
    plt.title("Déplacement d'une masse")
    plt.legend()
    plt.show()
    plt.close()


def courbes_euler_explicite(force):
    plt.plot(temps, euler_explicite(force,Niter)[0], "c*", label="Euler explicite")
    if force != "force":
        plt.plot(
            temps, solution_analytique(temps)[0], "r", label="Solution analytique",
        )
    elif force == "force":
        plt.plot(
            temps,
            solution_analytique_force(temps),
            "r",
            label="Solution analytique avec force",
        )
    plt.xlabel("Temps")
    plt.ylabel("Déplacement")
    plt.title("Déplacement d'une masse")
    plt.legend()
    plt.show()
    plt.close()


def courbes_euler_implicite(force):
    plt.plot(temps, euler_implicite(force,Niter)[0], "g|", label="Euler implicite")
    if force != "force":
        plt.plot(
            temps, solution_analytique(temps)[0], "r", label="Solution analytique",
        )
    elif force == "force":
        plt.plot(
            temps,
            solution_analytique_force(temps),
            "r",
            label="Solution analytique avec force",
        )
    plt.xlabel("Temps")
    plt.ylabel("Déplacement")
    plt.title("Déplacement d'une masse")
    plt.legend()
    plt.show()
    plt.close()


def courbes_trapèzes(force):
    plt.plot(temps, trapèzes(force,Niter)[0], "k--", label="Trapèzes")
    if force != "force":
        plt.plot(
            temps, solution_analytique(temps)[0], "r", label="Solution analytique",
        )
    elif force == "force":
        plt.plot(
            temps,
            solution_analytique_force(temps),
            "r",
            label="Solution analytique avec force",
        )
    plt.xlabel("Temps")
    plt.ylabel("Déplacement")
    plt.title("Déplacement d'une masse")
    plt.legend()
    plt.show()
    plt.close()


def courbes_solu_analytique(force):
    if force != "force":
        plt.plot(
            temps, solution_analytique(temps)[0], "r", label="Solution analytique",
        )
    elif force == "force":
        plt.plot(
            temps,
            solution_analytique_force(temps),
            "r",
            label="Solution analytique avec force",
        )
    plt.xlabel("Temps")
    plt.ylabel("Déplacement")
    plt.title("Déplacement d'une masse")
    plt.legend()
    plt.show()
    plt.close()


# Mettre un hashtag devant celle qu'on ne veut pas
print(courbes_comparaison("force"))
print(courbes_euler_explicite("force"))
print(courbes_euler_implicite("force"))
print(courbes_trapèzes("force"))
# print(courbes_solu_analytique(""))
