import numpy as np

m=1

k1 = 2 * (59 / 24)
k2 = 3 * (-7 / 6)
k3 = 4 * (7 / 48)
# k1=1
# k2=1
# k3=1
kse=0.0

# a1=1/2
# a2=1/2
# a3=1/2

a1 = 59 / 24
a2 = -7 / 6
a3 = 7 / 48

nombre_de_ressorts=20
nombre_de_points=100

u_0=0.01
v_0=0

l=1

omega_0 = np.sqrt(k1 / m)
T = 2 * np.pi / omega_0
print(T)

t_0 = 0
t_fin = T

delta_temps = (t_fin - t_0) / nombre_de_points
print("deltatemps",delta_temps)

U_0 = np.zeros((2 * nombre_de_ressorts, 1))
U_0[nombre_de_ressorts - 1, 0] = u_0
U_0[2 * nombre_de_ressorts - 1, 0] = v_0


Mat_Id = np.eye(2 * nombre_de_ressorts)


liste_u=[[0]for i in range(nombre_de_points)]
liste_v=[[0]for i in range(nombre_de_points)]

for i in range(nombre_de_points):
    liste_u[i] = [[u_0] for k in range(nombre_de_ressorts)] 
    liste_v[i] = [[v_0] for k in range(nombre_de_ressorts)] 
