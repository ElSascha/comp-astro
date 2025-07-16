import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# Arbeitsverzeichnis setzen
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Daten einlesen
data = np.loadtxt(r'../output/non_isothermal_flow_viscosity.txt', delimiter=';', skiprows=1)
# Parameter
num_steps = 2000    # Anzahl der Zeitschritte
num_points = 500     # Anzahl der Gitterzellen

# Daten vorbereiten
x_values = data[0:num_points, 1]  # x bleibt für alle Schritte gleich
density = data[:, 2].reshape((num_steps, num_points))
momentum = data[:, 3].reshape((num_steps, num_points))
internal_energy = data[:, 4].reshape((num_steps, num_points))

# Plot vorbereiten
fig, ax = plt.subplots(figsize=(10, 6))
line1, = ax.plot([], [], lw=2, label = 'density')
line2, = ax.plot([], [], lw=2, label = 'momentum')
line3, = ax.plot([], [], lw=2, label = 'internal energy')
ax.set_xlim(x_values[0], x_values[-1])
ax.set_ylim(1,2)
ax.set_xlabel("x")
ax.set_ylabel("Value")
ax.set_title("non isothermal flow with viscosity")
ax.legend()
ax.grid()

# Initialisierungsfunktion
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line1, line2, line3

# Animationsfunktion
def update(frame):
    current_time = data[frame * num_points, 0]  # Get the time for this frame
    ax.set_title(f"non isothermal flow with viscosity – t = {current_time:.2f} s")
    line1.set_data(x_values, density[frame])
    #line2.set_data(x_values, momentum[frame])
    #line3.set_data(x_values, internal_energy[frame])
    return line1, line2, line3

# Animation erstellen
ani = animation.FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=True, interval=100)
ani.save(r'../plots/isothermal_flow_viscosity_animation_density.mp4', writer='ffmpeg', fps=100)
#plt.tight_layout()
plt.show()
