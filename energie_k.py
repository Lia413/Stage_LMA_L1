from données import *
import matplotlib.pyplot as plt

x = np.linspace(-2, 6, nombre_de_points)
y = [[0] for i in range(len(x))]


def energie_pot():
    Ep = []
    a1 = 59 / 24
    a2 = -7 / 6
    a3 = 7 / 48
    for i in range(len(x)):
        Ep_i = a1 * x[i] ** 2 + a2 * x[i] ** 3 + a3 * x[i] ** 4
        Ep.append(Ep_i)

    return Ep


déplacement = np.linspace(t_0, t_fin, nombre_de_points)
plt.plot(x, energie_pot(), label="Ep_trap")
plt.plot(x, y, "r", label="y=0")
plt.xlabel("u")
plt.ylabel("y")
plt.title("Déplacement d'une masse d'un système bistable")
plt.legend()
intersection_points = []
for i in range(len(x)):
    if abs(energie_pot()[i] - y[i][0]) < 0.0001:
        intersection_points.append((x[i], y[i][0]))


if intersection_points:
    intersection_x, intersection_y = zip(*intersection_points)
    plt.scatter(
        intersection_x,
        intersection_y,
        marker="x",
        color="green",
        label="Points of intersection",
    )

    for point in intersection_points:
        plt.axvline(x=point[0], color="green", linestyle="--")

plt.show()
plt.close()
