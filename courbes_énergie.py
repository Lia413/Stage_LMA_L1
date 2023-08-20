from données import *
from ressort_Nddl import *
import matplotlib.pyplot as plt

temps = np.linspace(t_0, t_fin, nombre_de_points)


def energie(methode, N):
    plt.plot(temps, ressort_NDDL(methode, N, "", 0)[2], label="Ep")
    plt.plot(temps, ressort_NDDL(methode, N, "", 0)[3], label="Ec")
    plt.plot(temps, ressort_NDDL(methode, N, "", 0)[4], label="Em")
    plt.xlabel("temps")
    plt.ylabel("energies")
    if methode == "trapèzes":
        plt.title("diagramme énergétique par la méthode des trapèzes")
    if methode == "euler_implicite":
        plt.title("diagramme énergétique par la méthode d'Euler implicite")
    plt.legend()
    plt.show()
    plt.close()


energie("euler_implicite", nombre_de_points)
energie("trapèzes", nombre_de_points)
