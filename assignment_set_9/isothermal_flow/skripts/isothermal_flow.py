import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Arbeitsverzeichnis setzen
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Daten einlesen
data = np.loadtxt(r'../output/isothermal_flow.txt', delimiter=';', skiprows=1)
# Parameter
num_steps = 1000     # Anzahl der Zeitschritte
num_points = 500     # Anzahl der Gitterzellen

# Daten vorbereiten
x_values = data[0:num_points, 1]  # x bleibt für alle Schritte gleich
density = data[:, 2].reshape((num_steps, num_points))

# Plot vorbereiten
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], lw=2)
ax.set_xlim(x_values[0], x_values[-1])
ax.set_ylim(np.min(density), np.max(density))
ax.set_xlabel("x")
ax.set_ylabel("Dichte")
ax.set_title("Isotherme Gasdichte über Zeit")

# Initialisierungsfunktion
def init():
    line.set_data([], [])
    return line,

# Animationsfunktion
def update(frame):
    line.set_data(x_values, density[frame])
    ax.set_title(f"Isotherme Gasdichte – Zeitschritt {frame}")
    return line,

# Animation erstellen
ani = animation.FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=True, interval=100)
ani.save(r'../plots/isothermal_flow_animation.mp4', writer='ffmpeg', fps=100)
plt.tight_layout()
plt.show()
