from données import *
from ressort_1DDL_lineaire3 import *
import matplotlib.pyplot as plt

temps = np.linspace(t_0, t_fin, nombre_de_points)


def euler_explicite_portrait(Niter):
    plt.plot(
        euler_explicite("",Niter)[0], euler_explicite("",Niter)[1], "c", label="euler explicite"
    )
    plt.plot(
        solution_analytique(temps)[0],
        solution_analytique(temps)[1],
        "r",
        linewidth=4,
        label="solution analytique",
    )
    plt.xlabel("deplacement")
    plt.ylabel("vitesse")
    plt.title("Portrait de phase d'un systèmes masses-ressorts à 1DDL")
    plt.legend()
    plt.show()
    plt.close()


def euler_implicite_portrait(Niter):
    plt.plot(
        euler_implicite("",Niter)[0], euler_implicite("",Niter)[1], "g", label="euler implicite"
    )
    plt.plot(
        solution_analytique(temps)[0],
        solution_analytique(temps)[1],
        "r",
        linewidth=4,
        label="solution analytique",
    )
    plt.xlabel("deplacement")
    plt.ylabel("vitesse")
    plt.title("Portrait de phase d'un systèmes masses-ressorts à 1DDL")
    plt.legend()
    plt.show()
    plt.close()


def trapèzes_portrait(Niter):
    plt.plot(trapèzes("",Niter)[0], trapèzes("",Niter)[1], "ko", label="trapèzes")
    plt.plot(
        solution_analytique(temps)[0],
        solution_analytique(temps)[1],
        "r",
        label="solution analytique",
    )
    plt.xlabel("deplacement")
    plt.ylabel("vitesse")
    plt.title("Portrait de phase d'un systèmes masses-ressorts à 1DDL")
    plt.legend()
    plt.show()
    plt.close()


# mettre un hashtag devant les print que l'on ne veut pas
print(trapèzes_portrait(nombre_de_points))
print(euler_explicite_portrait(nombre_de_points))
print(euler_implicite_portrait(nombre_de_points))
