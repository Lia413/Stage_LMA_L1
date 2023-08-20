from données import *
import matplotlib.pyplot as plt

x = np.linspace(0, 4, nombre_de_points)
y = [[0] for _ in range(len(x))]

a1 = 59 / 24
a2 = -7 / 6
a3 = 7 / 48


def d_energie_pot():
    dEp = []

    for i in range(len(x)):
        dEp_i = 2 * a1 * x[i] + 3 * a2 * x[i] ** 2 + 4 * a3 * x[i] ** 3
        dEp.append(dEp_i)

    return dEp


plt.plot(x, d_energie_pot(), label="dEp_trap")
plt.plot(x, y, "r", label="y=0")
plt.xlabel("u")
plt.ylabel("y")
plt.title("Evolution de la dérivée de l'énergie")
plt.legend()


intersection_points = []
for i in range(len(x)):
    if abs(d_energie_pot()[i] - y[i][0]) < 0.005:
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

dEp_i = 2 * a1 * u_stable + 3 * a2 * u_stable ** 2 + 4 * a3 * u_stable ** 3
print(dEp_i)
plt.show()
plt.close()
