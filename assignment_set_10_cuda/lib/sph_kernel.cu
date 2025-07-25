#include "sph_kernel.cuh"

__device__ double cubic_bspline(double r, double h){
    double cons = NORMALIZATION_CONSTANT/pow(h, 3);
    if (0 <= r/h && r/h < 0.5) {
        return 8/(M_PI * h * h * h) * (6.0 * pow(r/h, 3) - 6.0 * pow(r/h, 2) + 1.0);
    } else if (0.5 <= r/h && r/h <= 1) {
        return 8/(M_PI * h * h * h) * (2.0 * pow(1.0 - r/h, 3));
    } else {
        return 0.0;
    }
}

__device__ double cubic_bspline_derivative(double r, double h){
    if (0 <= r/h && r/h < 0.5) {
        return 6*8/(M_PI * h * h * h * h) * (3 * pow(r/h, 2) - 2 * (r/h));
    } else if (0.5 <= r/h && r/h <= 1) {
        return -6 * 8/(M_PI * h * h * h * h) * pow(1.0 - r/h, 2);
    } else {
        return 0.0;
    }
}

__device__ double3 nabla_cubic_bspline(double3 pos_i, double3 pos_j, double h) {
    double r = length(sub(pos_i, pos_j));
    double factor = cubic_bspline_derivative(r, h) / r; // Normalize by r to get the gradient
    return scal_mul(sub(pos_i, pos_j), factor);
}

// instead of having 2 seperate loops, I will start a thread for each particle and compute the density and pressure in one go
// every thread will compute the density and pressure for one particle
// this will be more efficient than having 2 loops, as we can use shared memory to store the particles and avoid redundant memory accesses


__global__ void compute_density(ParticleData particles, int N, double smoothing_length) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of bounds
    particles.rho[i] = 0.0; // Initialize density to zero
    for (int j = 0; j < N; ++j) {
        if (i != j) {
            double r = length(sub(particles.pos[i], particles.pos[j]));
            particles.rho[i] += particles.mass[j] * cubic_bspline(r, smoothing_length);
        }
    }
    if (particles.rho[i] < 1e-12) particles.rho[i] = 1e-12; // Avoid division by zero
}

__global__ void compute_cs(ParticleData particles, int N, double GAMMA) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.cs[i] = sqrt(GAMMA * particles.pressure[i] / particles.rho[i]); // Compute sound speed using the equation cs = sqrt(gamma * pressure / density)
}

__global__ void compute_acceleration(ParticleData particles, int N, double h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return;

    particles.acc[i] = make_double3(0.0, 0.0, 0.0); // Reset acceleration

    for (int j = 0; j < N; ++j) {
        if (j != i) {
            double3 rij = sub(particles.pos[i], particles.pos[j]);
            double r = length(rij);
            if (r < 1e-6) continue;

            double3 gradW = nabla_cubic_bspline(particles.pos[i], particles.pos[j], h);

            
            double pij = (particles.pressure[i] / (particles.rho[i] * particles.rho[i]) +
                          particles.pressure[j] / (particles.rho[j] * particles.rho[j]));

            double3 acc_contrib = scal_mul(gradW, -particles.mass[j] * pij);

            particles.acc[i] = add(particles.acc[i], acc_contrib);
        }
    }
}

__global__ void compute_pressure(ParticleData particles, int N, double GAMMA, double K) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.pressure[i] = 0.0; // Initialize pressure to zero
    particles.pressure[i] = K * pow(particles.rho[i], GAMMA); // Compute pressure using the equation pressure = pi = K*ρ_i^γ
}

__global__ void compute_linear_acceleration_force(ParticleData particles, int N, double lambda) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.linear_acc_force[i] = make_double3(0.0, 0.0, 0.0); // Initialize linear acceleration force to zero
    particles.linear_acc_force[i] = scal_mul(particles.pos[i],-1*lambda); // Compute linear acceleration force using the equation a_i = -λ * x_i
}

__global__ void compute_damping_force(ParticleData particles, int N, double damping_coefficient) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.damping_force[i] = make_double3(0.0, 0.0, 0.0); // Initialize damping force to zero
    particles.damping_force[i] = scal_mul(particles.vel[i], -damping_coefficient); // Compute damping force using the equation f_i = -nu * v_i
}