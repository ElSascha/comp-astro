import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

def u0(x):
    return np.where(np.abs(x) < 1/3, 1.0, 0.0)


def analytical_solution_periodic(x, t, a=1.0, L=2.0):
    # L ist Länge des Gebiets, z.B. L=2 wenn x ∈ [-1,1]
    return u0(((x - a * t + L/2) % L) - L/2)
t = 4.0 # Time at which to evaluate the analytical solution
a = 1.0

# Set working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
print("Current working directory:", os.getcwd())

# Load data
centered_differencing = np.loadtxt(r'../centered_differencing_sigma_0.8.txt', delimiter=';')
upwind = np.loadtxt(r'../upwind_sigma_0.8.txt', delimiter=';')
lax_wendroff = np.loadtxt(r'../lax_wendroff_sigma_0.8.txt', delimiter=';')

# x-Werte
x = np.linspace(-1, 1, centered_differencing.shape[1])
# Analytische Lösung
u_analytical = analytical_solution_periodic(x, t, a, L=2.0)
# Set up figure and axes
fig, ax = plt.subplots(figsize=(10, 6))
#line1, = ax.plot(x, centered_differencing[999], label='Centered Differencing', color='blue')
#line2, = ax.plot(x, upwind[999], label='Upwind', color='red')
line3, = ax.plot(x, lax_wendroff[999], label='Lax-Wendroff', color='green')
line4, = ax.plot(x, u_analytical, label='Analytical Solution', color='black', linestyle='--')

ax.set_xlim(x.min(), x.max())
#ax.set_ylim(-0.1, 1.1)  # Passe an deine Daten an
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.set_title('Advection Equation Solutions (Animation)')
ax.legend()
ax.grid()
plt.show()

# # Init-Funktion
# def init():
#     line1.set_data([], [])
#     line2.set_data([], [])
#     line3.set_data([], [])
#     return line1, line2, line3

# # Update-Funktion für Frame i
# def update(i):
#     line1.set_data(x, centered_differencing[i])
#     #line2.set_data(x, upwind[i])
#     #line3.set_data(x, lax_wendroff[i])
#     return line1, line2, line3

# # Erzeuge Animation
# ani = animation.FuncAnimation(
#     fig, update, frames=1000, init_func=init, blit=True, interval=4  # 4ms ≈ 250 fps
# )

# # Zeige Animation
# plt.show()
# ani.save("advection_centered_differencing.mp4", fps=250)
