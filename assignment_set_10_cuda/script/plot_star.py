import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

Radius = 0.5 # Radius of the star
N = 1000  # Number of points in the plot
K = 0.1 
n = 1
lamda = 2.01203286081606
total_time = 20
time_step = 0.01
num_time_steps = int(total_time / time_step) 

def rho(r):
    """Density profile of the star."""
    return ((lamda/(2*K(1+n))) * (R**2 - r**2) )**n

#load data from .dat file
data = np.loadtxt(r"../data/particles.dat", skiprows=1, delimiter=';')
# time_step;x;y;z;vx;vy;vz;mass;pressure;rho;cs;linear_acc_fore_x;linear_acc_fore_y;linear_acc_fore_z;damping_force_x;damping_force_y;damping_force_z
positions = data[:, 1:4].reshape(num_time_steps, N, 3)
velocities = data[:, 4:7].reshape(num_time_steps, N, 3)
masses = data[:, 7].reshape(num_time_steps, N)
pressures = data[:, 8].reshape(num_time_steps, N)
densities = data[:, 9].reshape(num_time_steps, N)
cs = data[:, 10].reshape(num_time_steps, N)
linear_acc_fore = data[:, 11:14].reshape(num_time_steps, N, 3)
damping_force = data[:, 14:17].reshape(num_time_steps, N, 3)
# Plotting the star


fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')

scat = ax.scatter(positions[0,:,0], positions[0,:,1], positions[0,:,2], s=1)

ax.set_xlim(-Radius*1.2, Radius*1.2)
ax.set_ylim(-Radius*1.2, Radius*1.2)
ax.set_zlim(-Radius*1.2, Radius*1.2)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

def update(frame):
    scat._offsets3d = (positions[frame,:,0], positions[frame,:,1], positions[frame,:,2])
    current_time = frame * time_step
    ax.set_title(f'SPH Partikel - Stern Animation (3D)\nZeit: {current_time:.2f} s')
    return scat,

# interval so setzen, dass die Animation ca. 20 s läuft (in ms pro Frame)
interval_ms = (total_time / num_time_steps) * 1000

ani = FuncAnimation(fig, update, frames=num_time_steps, interval=interval_ms, blit=False)
plt.show()