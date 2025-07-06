import matplotlib.pyplot as plt
import numpy as np

gammas = [3, 5/3, 7/5]
n_values = [0.5, 1.5, 2.5] # selbst berechnet, um Rundugsfehler zu vermeiden
radien = [2.7527, 3.6538, 5.3553] # aus der main.cpp Datei

plt.figure(figsize=(10, 6))

for gamma, n, radius in zip(gammas, n_values, radien):

    filename = f"lane_emden_n{int(n*10)}.dat"
    data = np.loadtxt(filename, comments="#", delimiter=";")

    xi = data[:, 0]
    w = data[:, 1]
    rho = w**n # dimensionslos, da rho_c weggelassen

    # Berechnung Masse (siehe pdf für Berechnung warum das gilt)
    dw_dxi = (w[-1] - w[-2]) / (xi[-1] - xi[-2])
    mass = -radius**2 * dw_dxi

    plt.plot(xi, rho, label=f"$n = {n:.1f}$ $(\\gamma = {gamma:.2f})$, radius $R={radius:.3f}$, mass $M={mass:.3f}$")

plt.xlabel("$\\xi$")
plt.ylabel("$\\rho(\\xi) \propto w^n$")
plt.title("Star structure (density) for different adiabatic exponents")
plt.grid(True)
plt.legend()
plt.show()