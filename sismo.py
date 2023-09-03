from données import *

# from ressort_Nddl import *
from tester_forces import *
import matplotlib.pyplot as plt


def sismo(méthode, force, i):
    temps = np.linspace(t_0, t_fin, 3 * nombre_de_points)
    masses = ressort_NDDL(méthode, 3 * nombre_de_points, force, i)[0]
    for i in range(0, nombre_ressorts):
        for j in range(3 * nombre_de_points):
            masses[i][j] += 2*i  # perturbation d'une masse
        plt.plot(
                temps, masses[i],'c.', linewidth=0.2,
            )  
    plt.xlabel("temps")
    plt.ylabel("numéro de la masse")
    if méthode == "trapèzes":
        plt.title("Sismogramme par la méthode des trapèzes")
    elif méthode == "euler_implicite":
        plt.title("Sismogramme par la méthode d'Euler implicite")
    plt.show()
    plt.close("all")


print(sismo("trapèzes", "", 0))
print(sismo("euler_implicite", "", 0))
