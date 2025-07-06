import numpy as np
import matplotlib.pyplot as plt

# Analytische Lösungen
def analytical_w(xi, n):
    if n == 0:
        return 1 - (1/6) * xi**2
    elif n == 1:
        return np.sin(xi)/xi
    elif n == 5:
        return 1 / np.sqrt(1 + (xi**2) / 3)

# Plot
n_values = [0, 1, 5]

plt.figure(figsize=(10, 6))

for i, n in enumerate(n_values):
    # Lade numerische Daten
    filename = f"lane_emden_n{n}.dat"
    data = np.loadtxt(filename, comments="#", delimiter=";")
    xi_num, w_num = data[:, 0], data[:, 1]

    # Berechne analytische Lösung auf demselben xi-Bereich
    xi_ana = xi_num
    w_ana = analytical_w(xi_ana, n)

    # Plot numerisch und analytisch
    plt.plot(xi_num, w_num, label=f"Numerische Lsg. n={n}", linestyle="-")
    plt.plot(xi_ana, w_ana, label=f"Analytische Lsg. n={n}", color="red", linestyle=":")

plt.xlabel("$\\xi$")
plt.ylabel("$w(\\xi)$")
plt.title("Numerische und analytische Lösungen der Lane-Emden-Gleichung für $n=0,1,5$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()