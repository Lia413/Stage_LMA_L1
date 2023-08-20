from données import *
from cas_bi_stable import *
import matplotlib.pyplot as plt


temps = np.linspace(t_0, t_fin, nombre_de_points)


def Trapèzes():
    plt.plot(
        temps, positions_trap, label="Trapèzes ",
    )
    plt.xlabel("Temps (t)")
    plt.ylabel("Position (y)")
    plt.title("Déplacement d'une masse d'un système bistable - Trapèzes")
    plt.legend()
    plt.show()
    plt.close()


Trapèzes()


def Euler_implicite():
    plt.plot(
        temps, positions_euler, label="Euler implicite ",
    )
    plt.xlabel("Temps (t)")
    plt.ylabel("Position (y)")
    plt.title("Déplacement d'une masse d'un système bistable - Euler implicite")
    plt.legend()
    plt.show()
    plt.close()


Euler_implicite()




def diag_énergétique():
    plt.plot(temps, Em_tr, "k", label="Em - Newton Raphson - Trapèzes")
    plt.plot(temps, Ec_tr, "m", label="Ec - Newton Raphson - Trapèzes")
    plt.plot(temps, Ep_tr, "c", label="Ep - Newton Raphson - Trapèzes")

    plt.xlabel("Temps (t)")
    plt.ylabel("Énergie (y)")
    plt.title("Évolution énergétique d'un système bistable 1DDL - Trapèzes")
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.close()


diag_énergétique()
