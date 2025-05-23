import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Load the data
mercury = np.loadtxt(r"data/mercury_1000_ite.txt", skiprows=1, delimiter=	";")
halley = np.loadtxt(r"data/halleys_comet_1000_ite.txt", skiprows=1, delimiter=";")
# Create a DataFrame for better handling
mercury_df = pd.DataFrame(mercury, columns=["t", "x", "y", "r", "f", "E", "M", "P"])
halley_df = pd.DataFrame(halley, columns=["t", "x", "y", "r", "f", "E", "M", "P"])

plt.figure(figsize=(10, 5))
# Plot Halley's comet trajectory
plt.subplot(1, 2, 1)
plt.scatter(halley[:, 1], halley[:, 2], label="Halley's Comet", color="red", marker="x")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Halley's Comet Trajectory")
plt.legend()
plt.grid()
# Plot Mercury's trajectory
plt.subplot(1, 2, 2)
plt.scatter(mercury[:, 1], mercury[:, 2], label="Mercury", color="blue", marker="x")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Mercury Trajectory")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()




