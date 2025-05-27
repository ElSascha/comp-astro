import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
# Load the data
mercury = np.loadtxt(os.path.join(script_dir,r"..", r"data", r"mercury_fixed_point.txt"), skiprows=1, delimiter=	";")
halley = np.loadtxt(os.path.join(script_dir,r"..", r"data", r"halleys_comet_fixed_point.txt"), skiprows=1, delimiter=";")
# Create a DataFrame for better handling
mercury_df = pd.DataFrame(mercury, columns=["t", "x", "y", "r", "f", "E", "M", "P"])
halley_df = pd.DataFrame(halley, columns=["t", "x", "y", "r", "f", "E", "M", "P"])

plt.figure(figsize=(10, 5))
# Plot Halley's comet trajectory
plt.subplot(1, 2, 1)
plt.scatter(halley[:, 1], halley[:, 2], label="Halley's Comet", color="red", marker="x")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.title("Halley's Comet Trajectory")
plt.legend()
plt.grid()
# Plot Mercury's trajectory
plt.subplot(1, 2, 2)
plt.scatter(mercury[:, 1], mercury[:, 2], label="Mercury", color="blue", marker="x")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.title("Mercury Trajectory")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
# distance between Mercury and sun in days
plt.figure(figsize=(10, 5))
plt.scatter(mercury_df["t"], mercury_df["r"], label="Mercury", color="blue", marker="x")
plt.xlabel("Time (days)")
plt.ylabel("Distance from Sun (AU)")
plt.title("Distance of Mercury from the Sun Over Time")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
# distance between Halley's comet and sun in days
plt.figure(figsize=(10, 5))
plt.scatter(halley_df["t"], halley_df["r"], label="Halley's Comet", color="red", marker="x")
plt.xlabel("Time (days)")
plt.ylabel("Distance from Sun (AU)")
plt.title("Distance of Halley's Comet from the Sun Over Time")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()




