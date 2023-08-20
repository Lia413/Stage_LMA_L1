from données import *
import matplotlib.pyplot as plt

x = np.linspace(0, 4, nombre_de_points)
y = [[0] for i in range(len(x))]


def d_energie_pot():
    dEp = []

    a1 = 59 / 24
    a2 = -7 / 6
    a3 = 7 / 48

    for i in range(len(x)):
        dEp_i = 2 * a1 + 6 * a2 * x[i] + 12 * a3 * x[i] ** 2
        dEp.append(dEp_i)

    return dEp


déplacement = np.linspace(t_0, t_fin, nombre_de_points)
plt.plot(x, d_energie_pot(), label="d2Ep_trap")
plt.plot(x, y, "r", label="y=0")
plt.xlabel("u")
plt.ylabel("y")
plt.title("Evolution de la dérivée seconde de l'énergie")
plt.legend()


intersection_points = []
for i in range(len(x)):
    if abs(d_energie_pot()[i] - y[i][0]) < 0.001:
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
        plt.text(
            point[0],
            point[1],
            f"({point[0]:.2f})",
            verticalalignment="bottom",
            horizontalalignment="right",
            color="black",
        )

plt.show()
plt.close()
