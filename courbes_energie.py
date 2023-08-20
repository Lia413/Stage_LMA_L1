import matplotlib.pyplot as plt
import numpy as np
from smooth_extension import *


temps = np.linspace(t_0, t_fin, nombre_de_points)




def diag_énergétique_Ep(Ep):
    plt.plot(temps, Ep, label="Ep")


    plt.xlabel("Temps")
    plt.ylabel("Energie")
    plt.title("Evolution énergétique par la méthode des trapèzes")

    plt.legend()

    plt.show()


diag_énergétique_Ep(Ep)


def diag_énergétique(Ec, Ep, Em):
    plt.plot(temps, Ec, label="Ec")
    plt.plot(temps, Ep, label="Ep")
    plt.plot(temps, Em, "r", label="Em")

    plt.xlabel("Temps")
    plt.ylabel("Energie")
    plt.title("Evolution énergétique par la méthode des trapèzes")

    plt.legend()


    plt.show()


diag_énergétique(Ec, Ep, Em)


def tracer_F():
    t = np.linspace(t_0, t_fin - delta_temps, nombre_de_points - 1).reshape(-1, 1)


    liste_F_x = [F[nombre_de_ressorts] for F in liste_F]


    plt.plot(t, liste_F_x)
    plt.xlabel("Temps")
    plt.ylabel("Force")
    plt.title("Force en fonction du temps")
    plt.show()


tracer_F()
