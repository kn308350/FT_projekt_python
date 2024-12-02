#ANALIZA PRĘTA W 1D
import numpy as np
import matplotlib.pyplot as plt
import meshio

# Parametry pręta (to zależy od użytkownika)
length = 1.0       # Długość pręta [m]
n_elements = 10    # Liczba elementów
E = 200e9          # Moduł Younga [Pa]
A = 0.01           # Pole przekroju [m^2]
force = 1000      # Siła [N]

# Dyskretyzacja (punkty siatki)
x = np.linspace(0, length, n_elements + 1)

# Analiza MES
stiffness = E * A / (length / n_elements)  # Sztywność elementu
f = np.zeros(len(x))  # Wektor sił - wstępna macierz z samych 0
u = np.zeros(len(x))  # Wektor przemieszczeń - wstępna macierz z samych 0
K = np.zeros((len(x), len(x)))  # Macierz sztywności globalnej - wstępna macierz z samych 0

# Tworzenie macierzy sztywności - na podstawie mniejszych (lokalnych macierzy sztywnosci)
for i in range(n_elements):
    K[i, i] += stiffness
    K[i, i + 1] -= stiffness
    K[i + 1, i] -= stiffness
    K[i + 1, i + 1] += stiffness

# Siła, która jest przyłożona na końcu pręta (w zaleznoscie gdzie będzie przyłożona, będzie różny indeks)
f[-1] = force

# Warunki brzegowe: u[0] = 0 (dla uwierdzonego początku pręta)
K = K[1:, 1:]  # Usuń pierwszy wiersz i kolumnę (wiązanie w 0)
f = f[1:]      # Usuń siłę w zerowym węźle

# Rozwiązanie układu równań
u[1:] = np.linalg.solve(K, f)

# Naprężenia w elementach
strain = np.diff(u) / (length / n_elements)  # Deformacja
stress = E * strain  # Naprężenie

# Tworzymy siatkę 1D
cells = [("line", np.array([[i, i + 1] for i in range(n_elements)]))]
mesh = meshio.Mesh(points=x[:, None], cells=cells)
#cells to lista, która definiuje elementy siatki łączące węzły. "line" -> typ elementu (tutaj 1D)
#np.array([[i, i + 1] for i in range(n_elements)]) -> tablica definiująca połaczenia węzłów w postaci [węzeł i, węzeł i+1]
#mesh to fukncja z biblioteki meshio. x -. współrzedne węzłów (definicja linijka 13)
# Dodanie przemieszczeń i naprężeń jako danych punktowych i elementowych
mesh.point_data["displacement"] = u
#point_data -> słownik w mesh; związany z węzłami siatki - przypisuje im przemieszczenie -> "displacement klucz w słowniku
mesh.cell_data["stress"] = [stress]
#cell_data -. słownik; związany z elementami siatki, "stress" klucz w słowniku, [stress] -> lista naprężeń

# Wizualizacja - wykresy w matplotlib
plt.figure(figsize=(10, 5))

# Przemieszczenia
plt.subplot(1, 2, 1)
plt.plot(x, u, marker='o', label='Przemieszczenia')
plt.xlabel('Długość pręta [m]')
plt.ylabel('Przemieszczenie [m]')
plt.title('Rozkład przemieszczeń')
plt.grid()
plt.legend()

# Naprężenia
plt.subplot(1, 2, 2)
x_mid = (x[:-1] + x[1:]) / 2  # Środek każdego elementu
plt.plot(x_mid, stress, marker='s', color='red', label='Naprężenia')
plt.xlabel('Długość pręta [m]')
plt.ylabel('Naprężenie [Pa]')
plt.title('Rozkład naprężeń')
plt.grid()
plt.legend()

# Wyświetlenie wykresów
plt.tight_layout()
plt.show()
