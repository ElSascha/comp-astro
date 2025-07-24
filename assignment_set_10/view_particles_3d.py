#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_particle_data(filename):
    """Load particle data from file"""
    data = []
    time = 0.0
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('# Time:'):
                time = float(line.split(':')[1].strip())
            elif not line.startswith('#'):
                values = list(map(float, line.strip().split()))
                if len(values) == 7:  # x, y, z, vx, vy, vz, m
                    data.append(values)
    
    if data:
        data = np.array(data)
        return time, data[:, :3], data[:, 3:6], data[:, 6]  # positions, velocities, masses
    else:
        return time, np.array([]), np.array([]), np.array([])

def interactive_3d_plot(filename="particle_data_final.dat"):
    """Create an interactive 3D plot"""
    time, positions, velocities, masses = load_particle_data(filename)
    
    if len(positions) == 0:
        print(f"No data found in {filename}")
        return
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Scatter plot of particles, color-coded by mass
    scatter = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                        c=masses, cmap='plasma', s=30, alpha=0.8)
    
    # Add velocity vectors (subsample for clarity)
    subsample = slice(None, None, 20)  # Show every 20th velocity vector
    ax.quiver(positions[subsample, 0], positions[subsample, 1], positions[subsample, 2],
              velocities[subsample, 0], velocities[subsample, 1], velocities[subsample, 2],
              length=0.15, alpha=0.7, color='red', arrow_length_ratio=0.1)
    
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title(f'3D SPH Particle System at t = {time:.2f}\n(Red arrows show velocity vectors)')
    
    # Add colorbar for mass
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.1)
    cbar.set_label('Particle Mass')
    
    # Improve the appearance
    ax.grid(True, alpha=0.3)
    
    # Set equal aspect ratio if possible
    max_range = np.array([positions[:,0].max()-positions[:,0].min(),
                         positions[:,1].max()-positions[:,1].min(),
                         positions[:,2].max()-positions[:,2].min()]).max() / 2.0
    mid_x = (positions[:,0].max()+positions[:,0].min()) * 0.5
    mid_y = (positions[:,1].max()+positions[:,1].min()) * 0.5
    mid_z = (positions[:,2].max()+positions[:,2].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Interactive 3D Particle Viewer")
    print("Showing final state of particle system...")
    interactive_3d_plot()
