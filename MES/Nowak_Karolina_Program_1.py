import numpy as np
import matplotlib.pyplot as plt
"""
Uzytkownik podaje wartosci parametrow takie jak dlugosc belki, liczba elementow do analizy, modul Younga, srednice 
przekroju kolowego, dlugosc boku przekroju kwadratowego oraz sile przylozona w srodku belki. Program dokonuje analizy,
ktora belka odksztalci sie mniej przy podanych warunkach i wybiera ja jako optymalna. Dodatkowo wyswietla wartosc 
maksymalnego przemieszczenia pionowego oraz wykres z porownaniem ugiecia belki dla obu przekrojow.
"""
# Pobranie danych od użytkownika z ustawieniem wartosci domyslnych
try:
    BLen = float(input("Podaj długość belki [mm] (domyślnie 500): ") or 500)
    NEl = int(input("Podaj liczbę elementów (musi być parzysta, domyślnie 6): ") or 6)
    if NEl % 2 != 0:
        raise ValueError("Liczba elementów musi być parzysta!")
    E_input = input("Podaj moduł Younga [Pa] (domyślnie 2.1e5): ")
    E = float(E_input) if E_input else 2.1e5

    dimension_circle = float(input("Podaj średnicę dla przekroju kołowego [mm] (domyślnie 10): ") or 10)
    dimension_square = float(input("Podaj długość boku dla przekroju kwadratowego [mm] (domyślnie 20): ") or 20)

    F_mid_input = input("Podaj siłę przyłożoną w środku belki [N] (domyślnie 1000): ")
    F_mid = float(F_mid_input) if F_mid_input else 1000
except ValueError as e:
    print(f"Błąd: {e}")
    exit()

##Czesc obliczeniowa kodu

# Współrzędne węzłów
XNod = np.linspace(0, BLen, NEl + 1)

NNod = len(XNod)  # Liczba węzłów
Dof = 2 * NNod  # Liczba stopni swobody

# Definiowanie elementów
Elem = np.array([[j, j + 1] for j in range(1, NEl + 1)])

# Podpory zgodnie z zalozeniem, ze belka wsparta jest podporami na obu koncach
Sup = np.zeros(Dof, dtype=int)  # Inicjalizacja
Sup[0] = 1  # Pierwszy węzeł
Sup[Dof - 2] = 1  # Ostatni węzeł

# Wektor obciążeń
F = np.zeros(Dof)
F[Dof // 2] = F_mid  # Obciążenie w środku belki

# Siła trakcyjna
q = np.full(NEl, 2)


# Funkcja do obliczania momentu bezwładności
# Typy przekroju: 'circle' (kołowy) lub 'square' (kwadratowy)
def calculate_inertia(cross_section, dimension):
    """
    Funkcja oblicza wartosc momentu bezwladnosci dla obu podanych przekorjow w zaleznosci od podanych wymiarow
    """
    if cross_section == 'circle':
        radius = dimension / 2
        return (np.pi * radius ** 4) / 4
    elif cross_section == 'square':
        return dimension ** 4 / 12
    else:
        raise ValueError("Nieznany typ przekroju")


# Definiowanie właściwości przekroju
cross_sections = ['circle', 'square']
dimensions = [dimension_circle, dimension_square]

# Przechowywanie wyników
results = {}
displacement_data = {}

for cross_section, dimension in zip(cross_sections, dimensions):
    # Moment bezwładności przekroju
    I = np.array([calculate_inertia(cross_section, dimension) for _ in range(NEl)])

    # Obliczenia
    Le = np.zeros(NEl)
    for j in range(NEl):
        Le[j] = abs(XNod[Elem[j, 0] - 1] - XNod[Elem[j, 1] - 1])
        F[Elem[j, 0] * 2 - 2] += 0.5 * q[j] * Le[j]
        F[Elem[j, 1] * 2 - 2] += 0.5 * q[j] * Le[j]

    # Inicjalizacja globalnej macierzy sztywności
    K = np.zeros((Dof, Dof))

    # Agregacja macierzy
    for j in range(NEl):
        L = Le[j]
        k = np.zeros((4, 4))
        k[0, 0] = 12 / L ** 3
        k[0, 1] = 6 / L ** 2
        k[0, 2] = -12 / L ** 3
        k[0, 3] = 6 / L ** 2
        k[1, 1] = 4 / L
        k[1, 2] = -6 / L ** 2
        k[1, 3] = 2 / L
        k[2, 2] = 12 / L ** 3
        k[2, 3] = -6 / L ** 2
        k[3, 3] = 4 / L
        k = k + np.triu(k, 1).T
        k = k * E * I[j]

        index1 = np.arange(Elem[j, 0] * 2 - 2, Elem[j, 0] * 2)
        index2 = np.arange(Elem[j, 1] * 2 - 2, Elem[j, 1] * 2)
        K[np.ix_(index1, index1)] += k[:2, :2]
        K[np.ix_(index2, index2)] += k[2:, 2:]
        K[np.ix_(index1, index2)] += k[:2, 2:]
        K[np.ix_(index2, index1)] += k[2:, :2]

    # Warunki brzegowe
    K_copy = K.copy()
    F_copy = F.copy()
    for j in range(Dof):
        if Sup[j]:
            K_copy[j, :] = 0
            K_copy[:, j] = 0
            K_copy[j, j] = 1
            F_copy[j] = 0

    # Rozwiązanie układu równań
    U = np.linalg.solve(K_copy, F_copy)

    # Maksymalne przemieszczenie
    Umax = np.max(U[::2])

    # Zapis wyników
    results[cross_section] = Umax
    displacement_data[cross_section] = U

# Wybór optymalnego przekroju
optimal_cross_section = min(results, key=results.get)
optimal_displacement = results[optimal_cross_section]

# Wyświetlenie wyników
print(f"Optymalny przekrój: {optimal_cross_section}")
print(f"Maksymalne przemieszczenie pionowe: {optimal_displacement} mm")

# Wizualizacja dla obu przekrojów
plt.figure()
for cross_section in cross_sections:
    plt.plot(XNod, displacement_data[cross_section][::2], 'o-', label=f"Ugięcie belki ({cross_section})")

plt.xlabel("x [mm]")
plt.ylabel("y [mm]")
plt.gca().invert_yaxis()
plt.title("Belka - przemieszczenia dla różnych przekrojów")
plt.legend()
plt.grid()
plt.show()

