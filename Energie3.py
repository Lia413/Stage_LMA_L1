from données import *
from ressort_1DDL_lineaire3 import *
import matplotlib.pyplot as plt

# bool=True =>Ec
# bool=False =>Ep
déplacement = np.linspace(t_0, t_fin, nombre_de_points)


def energie_analytique():
    plt.plot(déplacement, solution_analytique(déplacement)[2], label="Ec_analytique")
    plt.plot(déplacement, solution_analytique(déplacement)[3], label="Ep_analytique")
    plt.plot(déplacement, solution_analytique(déplacement)[4], label="Em_analytique")
    plt.xlabel("temps")
    plt.ylabel("energies")
    plt.title("diagramme énergétique de la solution analytique")
    plt.legend()
    plt.show()
    plt.close()


def energie_euler_explicite(force):
    plt.plot(déplacement, energie_ressort("euler_explicite", True,force,nombre_de_points)[0], label="Ec")
    plt.plot(déplacement, energie_ressort("euler_explicite", False,force,nombre_de_points)[1], label="Ep")
    plt.plot(déplacement, energie_mécanique(euler_explicite,force,nombre_de_points), label="Em")
    plt.xlabel("temps")
    plt.ylabel("énergies")
    plt.title("diagramme énergétique avec la méthode euler_explicite")
    plt.legend()
    plt.show()
    plt.close()


def energie_euler_implicite(force):
    plt.plot(
        déplacement, energie_ressort("euler_implicite", True, force,nombre_de_points)[0], label="Ec"
    )
    plt.plot(
        déplacement, energie_ressort("euler_implicite", False, force,nombre_de_points)[1], label="Ep"
    )
    plt.plot(déplacement, energie_mécanique(euler_implicite, force,nombre_de_points), label="Em")
    plt.xlabel("temps")
    plt.ylabel("energies")
    plt.title("diagramme énergétique avec la méthode euler_implicite")
    plt.legend()
    plt.show()
    plt.close()


def energie_trapèzes(force):
    plt.plot(déplacement, energie_ressort("trapèzes", True,nombre_de_points,nombre_de_points)[0], label="Ec")
    plt.plot(déplacement, energie_ressort("trapèzes", False,nombre_de_points,nombre_de_points)[1], label="Ep")
    plt.plot(déplacement, energie_mécanique(trapèzes,force,nombre_de_points), label="Em")
    plt.xlabel("temps")
    plt.ylabel("énergies")
    plt.title("diagramme énergétique avec la méthode des trapèzes")
    plt.legend()
    plt.show()
    plt.close()


def energie_meca_comparaison(force):
    plt.plot(déplacement, energie_mécanique(trapèzes, force,nombre_de_points), "k+", label="Em_trapèzes")
    plt.plot(
        déplacement,
        energie_mécanique(euler_implicite, force,nombre_de_points),
        "*",
        label="Em_implicite",
    )
    plt.plot(
        déplacement,
        energie_mécanique(euler_explicite, force,nombre_de_points),
        "_",
        label="Em_explicite",
    )
    plt.plot(déplacement, solution_analytique(déplacement)[4], label="Em_analytique")

    plt.xlabel("temps")
    plt.ylabel("énergies")
    plt.title("diagramme énergétique avec la méthode des trapèzes")
    plt.legend()
    plt.show()
    plt.close()


# mettre un hasthag devant celle qu'on ne veut pas
print(energie_trapèzes(""))
print(energie_euler_explicite(""))
print(energie_analytique())
print(energie_euler_implicite(""))
print(energie_meca_comparaison(""))
