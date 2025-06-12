import matplotlib.pyplot as plt
import numpy as np


x0s = [0.21, 0.24, 0.26, 0.27, 0.4, 0.5, 0.6, 0.8]



for x0 in x0s:
    data = np.loadtxt(f"poincare_x0_{x0:.6f}.csv", delimiter=",")
    plt.figure(figsize=(10, 6))
    plt.scatter(data[:, 0], data[:, 1], label=f"x0 = {x0:.2f}")
    plt.title("$x - \\dot{x}$ for different initial conditions")
    plt.xlabel("x")
    plt.ylabel("$\\dot{x}$")
    plt.legend()
    plt.savefig(f"plots_poincare/poincare_x0_{x0:.2f}.png")