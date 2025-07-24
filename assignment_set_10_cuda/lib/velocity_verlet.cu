#include "velocity_verlet.cuh"


__global__ void update_position(ParticleData particles, int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of bounds

    // Update position using the velocity and acceleration
    particles.pos[i] = add(add(particles.pos[i], scal_mul(particles.vel[i], dt)),scal_mul(particles.acc[i], 0.5 * dt * dt));
}

__global__ void update_velocity(ParticleData particles, int N, double dt, double3* acc_old) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of bounds

    // Update velocity using the acceleration
    particles.vel[i] = add(particles.vel[i], scal_mul(add(acc_old[i], particles.acc[i]), 0.5 * dt));
}

void time_step(ParticleData& particles, double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient, double smoothing_length) {
    int blocks = (N + 255) / 256; // Calculate number of blocks needed
    int threads = 256; // Number of threads per block
    // 1. Allocate device memory for acc_old
    double3* acc_old;
    cudaMalloc(&acc_old, N * sizeof(double3));
    // 2. Copy current acceleration to acc_old (device to device copy)
    cudaMemcpy(acc_old, particles.acc, N * sizeof(double3), cudaMemcpyDeviceToDevice);
    // Calculate new positions
    update_position<<<blocks, threads>>>(particles, N, dt);
    cudaDeviceSynchronize();
    // Update acceleration
    compute_acceleration<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    // Calculate new velocities
    update_velocity<<<blocks, threads>>>(particles, N, dt, acc_old);
    cudaDeviceSynchronize();
    // Free old acceleration memory
    cudaFree(acc_old);
    // Compute density, pressure, sound speed, linear acceleration, and damping force
    compute_density<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_pressure<<<blocks, threads>>>(particles, N, GAMMA, K);
    cudaDeviceSynchronize();
    compute_cs<<<blocks, threads>>>(particles, N, GAMMA);
    cudaDeviceSynchronize();
    compute_linear_acceleration_fore<<<blocks, threads>>>(particles, N, lambda);
    cudaDeviceSynchronize();
    compute_damping_force<<<blocks, threads>>>(particles, N, damping_coefficient);
    cudaDeviceSynchronize();
}