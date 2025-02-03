import numpy as np
"""
Uzytkownik moze wprowadzic wartosci dlugosci belki, sily w srodkowym wezle, modulu Younga.
Program oblicza optymalne wymiary dla przekroju kolowego oraz kwadratowego po podaniu maksymalnego dopuszczalnego
przemieszczenia pionowego. Jako optymalny zostaje przyjety najmniejszy przekroj belki ktory spelni zalozenie 
o maksymalnym dopuszczalnym przemieszczeniu. Program zwraca wymiary przekroju z dokladnoscia 2 miejsc po przecinku.

"""
# Pobranie danych od użytkownika
dlugosc_belki = float(input("Podaj długość belki [mm] (domyślnie 500mm): ") or 500)
sila_w_srodku = float(input("Podaj siłę w środkowym węźle [N] (domyślnie 1000N): ") or 1000)
modul_Younga = float(input("Podaj moduł Younga [Pa] (domyślnie 2.1e5): ") or 210000)

# Liczba elementów (musi być parzysta)
NEl = 8

# Współrzędne węzłów
XNod = np.linspace(0, dlugosc_belki, NEl + 1)

NNod = len(XNod)  # Liczba węzłów
Dof = 2 * NNod  # Liczba stopni swobody

# Definicja elementów
Elem = np.array([[j, j + 1] for j in range(1, NEl + 1)])

# Moduł Younga dla każdego elementu
E = np.full(NEl, modul_Younga)

# Podpory
Sup = np.zeros(Dof, dtype=int)  # Inicjalizacja
Sup[0] = 1  # Pierwszy węzeł
Sup[Dof - 2] = 1  # Ostatni węzeł

# Wektor obciążeń
F = np.zeros(Dof)
F[Dof // 2] = sila_w_srodku  # Środkowy węzeł

# Siła rozłożona (obciążenie trakcyjne)
q = np.full(NEl, 2)


# Funkcja do obliczania momentu bezwładności
def calculate_inertia(cross_section, dimension):
    """
    Funkcja oblicza moment bezwladnosci dla obu rozpatrywanych przekrojow - kolowego oraz kwadratowego.

    """
    if cross_section == 'circle':
        return (np.pi * (dimension / 2) ** 4) / 4
    elif cross_section == 'square':
        return dimension ** 4 / 12
    else:
        raise ValueError("Nieznany typ przekroju")


# Funkcja do obliczania przemieszczeń
def calculate_displacement(cross_section, dimension):
    I = np.array([calculate_inertia(cross_section, dimension) for _ in range(NEl)])

    Le = np.zeros(NEl)
    for j in range(NEl):
        Le[j] = abs(XNod[Elem[j, 0] - 1] - XNod[Elem[j, 1] - 1])
        F[Elem[j, 0] * 2 - 2] += 0.5 * q[j] * Le[j]
        F[Elem[j, 1] * 2 - 2] += 0.5 * q[j] * Le[j]

    K = np.zeros((Dof, Dof))
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
        k = k * E[j] * I[j]
        index1 = np.arange(Elem[j, 0] * 2 - 2, Elem[j, 0] * 2)
        index2 = np.arange(Elem[j, 1] * 2 - 2, Elem[j, 1] * 2)
        K[np.ix_(index1, index1)] += k[:2, :2]
        K[np.ix_(index2, index2)] += k[2:, 2:]
        K[np.ix_(index1, index2)] += k[:2, 2:]
        K[np.ix_(index2, index1)] += k[2:, :2]

    K_copy = K.copy()
    F_copy = F.copy()
    for j in range(Dof):
        if Sup[j]:
            K_copy[j, :] = 0
            K_copy[:, j] = 0
            K_copy[j, j] = 1
            F_copy[j] = 0

    U = np.linalg.solve(K_copy, F_copy)
    return U, np.max(U[::2])


# Pobranie maksymalnego dopuszczalnego przemieszczenia
max_displacement = float(input("Podaj maksymalne dopuszczalne przemieszczenie [mm]: "))

cross_sections = ['circle', 'square']
optimal_dimensions = {}
displacements = {}

dimension_range = np.linspace(1, 100, 1000)
for cross_section in cross_sections:
    optimal_dimensions[cross_section] = None
    displacements[cross_section] = None
    for dimension in dimension_range:
        _, Umax = calculate_displacement(cross_section, dimension)
        if Umax <= max_displacement:
            optimal_dimensions[cross_section] = dimension
            displacements[cross_section] = Umax
            break

# Wyświetlenie wyników
for cross_section in cross_sections:
    if optimal_dimensions[cross_section] is not None:
        print(f"Optymalne wymiary dla {cross_section}: {optimal_dimensions[cross_section]:.2f} [mm]")
        print(f"Maksymalne przemieszczenie dla {cross_section}: {displacements[cross_section]:.2f} [mm]")
    else:
        print(f"Brak odpowiednich wymiarów dla {cross_section} z podanymi warunkami.")
