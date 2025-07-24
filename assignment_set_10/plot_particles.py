#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob
import os
import re
import matplotlib.animation as animation

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

def plot_3d_particles(positions, velocities, masses, time, save_name=None, axis_limits=None):
    """Create a 3D plot of particle positions with velocity vectors"""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Scatter plot of particles, color-coded by mass
    scatter = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                        c=masses, cmap='viridis', s=20, alpha=0.7, 
                        vmin=masses.min(), vmax=masses.max())
    
    # Add velocity vectors (subsample for clarity)
    subsample = slice(None, None, 10)  # Show every 10th velocity vector
    ax.quiver(positions[subsample, 0], positions[subsample, 1], positions[subsample, 2],
              velocities[subsample, 0], velocities[subsample, 1], velocities[subsample, 2],
              length=0.1, alpha=0.6, color='red', arrow_length_ratio=0.1)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'3D Particle Distribution at t = {time:.2f}')
    
    # Set consistent axis limits if provided
    if axis_limits is not None:
        ax.set_xlim(axis_limits[0])
        ax.set_ylim(axis_limits[1])
        ax.set_zlim(axis_limits[2])
    
    # Add colorbar for mass
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('Mass')
    
    # Make the plot look nice
    ax.grid(True, alpha=0.3)
    
    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {save_name}")
    
    return fig, ax

def sort_data_files(filenames):
    """Sort data files by their time step number"""
    def extract_step_number(filename):
        # Extract step number from filename like "particle_data_t123.dat"
        if "t0.dat" in filename:
            return 0
        elif "final.dat" in filename:
            return float('inf')  # Put final at the end
        else:
            match = re.search(r't(\d+)\.dat', filename)
            return int(match.group(1)) if match else 0
    
    return sorted(filenames, key=extract_step_number)

def get_global_axis_limits():
    """Calculate global axis limits from all data files"""
    all_files = glob.glob("particle_data_*.dat")
    
    x_min, x_max = float('inf'), float('-inf')
    y_min, y_max = float('inf'), float('-inf')
    z_min, z_max = float('inf'), float('-inf')
    
    for filename in all_files:
        _, positions, _, _ = load_particle_data(filename)
        if len(positions) > 0:
            x_min = min(x_min, positions[:, 0].min())
            x_max = max(x_max, positions[:, 0].max())
            y_min = min(y_min, positions[:, 1].min())
            y_max = max(y_max, positions[:, 1].max())
            z_min = min(z_min, positions[:, 2].min())
            z_max = max(z_max, positions[:, 2].max())
    
    # Add some padding
    padding = 0.1
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    
    return ([x_min - padding * x_range, x_max + padding * x_range],
            [y_min - padding * y_range, y_max + padding * y_range],
            [z_min - padding * z_range, z_max + padding * z_range])

def create_animation_frames():
    """Create individual frames for all available data files"""
    # Find all particle data files and sort them properly
    data_files = sort_data_files(glob.glob("particle_data_*.dat"))
    
    if not data_files:
        print("No particle data files found. Run the C++ simulation first.")
        return
    
    print(f"Found {len(data_files)} data files")
    print("Processing files in order:")
    for i, f in enumerate(data_files):
        print(f"  {i:2d}: {f}")
    
    # Get global axis limits for consistent scaling
    axis_limits = get_global_axis_limits()
    print(f"Global axis limits: X{axis_limits[0]}, Y{axis_limits[1]}, Z{axis_limits[2]}")
    
    # Create output directory for plots
    os.makedirs("plots", exist_ok=True)
    
    for i, filename in enumerate(data_files):
        print(f"Processing {filename}...")
        time, positions, velocities, masses = load_particle_data(filename)
        
        if len(positions) > 0:
            plot_name = f"plots/particles_frame_{i:03d}.png"
            fig, ax = plot_3d_particles(positions, velocities, masses, time, plot_name, axis_limits)
            plt.close(fig)  # Close to save memory
    
    print(f"Created {len(data_files)} plot frames in the 'plots' directory")
    return data_files

    """Create comparison plot of initial and final states"""
    # Get global axis limits for consistent scaling
    axis_limits = get_global_axis_limits()
    
    # Load initial data
    if os.path.exists("particle_data_t0.dat"):
        time_init, pos_init, vel_init, mass_init = load_particle_data("particle_data_t0.dat")
    else:
        print("Initial data file not found")
        return
    
    # Load final data
    if os.path.exists("particle_data_final.dat"):
        time_final, pos_final, vel_final, mass_final = load_particle_data("particle_data_final.dat")
    else:
        print("Final data file not found")
        return
    
    # Create side-by-side comparison
    fig = plt.figure(figsize=(20, 8))
    
    # Initial state
    ax1 = fig.add_subplot(121, projection='3d')
    scatter1 = ax1.scatter(pos_init[:, 0], pos_init[:, 1], pos_init[:, 2], 
                          c=mass_init, cmap='viridis', s=20, alpha=0.7,
                          vmin=mass_init.min(), vmax=mass_init.max())
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title(f'Initial State (t = {time_init:.2f})')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(axis_limits[0])
    ax1.set_ylim(axis_limits[1])
    ax1.set_zlim(axis_limits[2])
    
    # Final state
    ax2 = fig.add_subplot(122, projection='3d')
    scatter2 = ax2.scatter(pos_final[:, 0], pos_final[:, 1], pos_final[:, 2], 
                          c=mass_final, cmap='viridis', s=20, alpha=0.7,
                          vmin=mass_final.min(), vmax=mass_final.max())
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_title(f'Final State (t = {time_final:.2f})')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(axis_limits[0])
    ax2.set_ylim(axis_limits[1])
    ax2.set_zlim(axis_limits[2])
    
    # Add colorbar
    cbar = plt.colorbar(scatter2, ax=[ax1, ax2], shrink=0.8)
    cbar.set_label('Mass')
    
    plt.tight_layout()
    plt.savefig("particle_evolution_comparison.png", dpi=300, bbox_inches='tight')
    plt.show()
    print("Saved comparison plot to particle_evolution_comparison.png")

def create_mp4_animation():
    """Create an MP4 animation from the particle data"""
    try:
        import matplotlib.animation as animation
        from matplotlib.animation import FFMpegWriter
    except ImportError:
        print("Error: matplotlib.animation not available. Install ffmpeg: sudo apt install ffmpeg")
        return
    
    # Find all particle data files and sort them properly
    data_files = sort_data_files(glob.glob("particle_data_*.dat"))
    
    if not data_files:
        print("No particle data files found. Run the C++ simulation first.")
        return
    
    print(f"Creating animation from {len(data_files)} frames...")
    
    # Get global axis limits for consistent scaling
    axis_limits = get_global_axis_limits()
    
    # Load all data first
    all_data = []
    for filename in data_files:
        time, positions, velocities, masses = load_particle_data(filename)
        if len(positions) > 0:
            all_data.append((time, positions, velocities, masses))
    
    if not all_data:
        print("No valid data found")
        return
    
    # Create figure and axis
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get global mass range for consistent coloring
    all_masses = np.concatenate([data[3] for data in all_data])
    mass_min, mass_max = all_masses.min(), all_masses.max()
    
    def animate(frame_num):
        ax.clear()
        time, positions, velocities, masses = all_data[frame_num]
        
        # Scatter plot of particles
        scatter = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                           c=masses, cmap='viridis', s=20, alpha=0.7,
                           vmin=mass_min, vmax=mass_max)
        
        # Add velocity vectors (subsample for clarity)
        subsample = slice(None, None, 20)  # Show every 20th velocity vector
        ax.quiver(positions[subsample, 0], positions[subsample, 1], positions[subsample, 2],
                  velocities[subsample, 0], velocities[subsample, 1], velocities[subsample, 2],
                  length=0.1, alpha=0.6, color='red', arrow_length_ratio=0.1)
        
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_zlabel('Z Position')
        ax.set_title(f'3D SPH Particle System - t = {time:.2f}s')
        
        # Set consistent axis limits
        ax.set_xlim(axis_limits[0])
        ax.set_ylim(axis_limits[1])
        ax.set_zlim(axis_limits[2])
        
        ax.grid(True, alpha=0.3)
        
        return scatter,
    
    # Create animation
    print("Generating animation...")
    anim = animation.FuncAnimation(fig, animate, frames=len(all_data), 
                                 interval=200, blit=False, repeat=True)
    
    # Save as MP4
    print("Saving animation as particle_evolution.mp4...")
    try:
        writer = FFMpegWriter(fps=30, metadata=dict(artist='SPH Simulation'), bitrate=1800)
        anim.save('particle_evolution.mp4', writer=writer, dpi=150)
        print("Animation saved successfully!")
    except Exception as e:
        print(f"Error saving MP4: {e}")
        print("Trying to save as GIF instead...")
        try:
            anim.save('particle_evolution.gif', writer='pillow', fps=5, dpi=100)
            print("Animation saved as GIF!")
        except Exception as e2:
            print(f"Error saving GIF: {e2}")
    
    plt.close(fig)

def main():
    """Main function to create visualizations"""
    print("3D Particle Visualization Script")
    print("=" * 40)
    
    # Check if data files exist
    data_files = glob.glob("particle_data_*.dat")
    if not data_files:
        print("No particle data files found.")
        print("Please run the C++ simulation first to generate data files.")
        return
    
    print(f"Found {len(data_files)} data files")
    
    # Create comparison plot of initial and final states
    # print("\nCreating initial vs final comparison plot...")
    # plot_initial_and_final()
    
    # Create animation frames with proper sorting and consistent scaling
    print("\nCreating individual frame plots with consistent scaling...")
    sorted_files = create_animation_frames()
    
    # Create MP4 animation
    print("\nCreating MP4 animation...")
    create_mp4_animation()
    
    print("\nVisualization complete!")
    print("- Check 'particle_evolution_comparison.png' for initial vs final comparison")
    print("- Check 'plots/' directory for individual frame images (properly sorted)")
    print("- Check 'particle_evolution.mp4' for animation (or .gif if MP4 failed)")

if __name__ == "__main__":
    main()
