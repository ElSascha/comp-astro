#include "velocity_verlet.cuh"


__global__ void update_position(ParticleData particles, int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of boundsd
    // Update position using the Velocity Verlet algorithm
    double3 total_acc = add(add(particles.acc[i], particles.linear_acc_force[i]), particles.damping_force[i]);
    particles.pos[i] = add(particles.pos[i], add(scal_mul(particles.vel[i], dt), scal_mul(total_acc, 0.5 * dt * dt)));
}
__global__ void update_velocity(ParticleData particles, int N, double dt, double3* acc_old, double3* linear_acc_old, double3* damping_force_old) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of bounds
    double3 acc_avg = scal_mul(add(acc_old[i], particles.acc[i]), 0.5);
    double3 lin_acc_avg = scal_mul(add(linear_acc_old[i], particles.linear_acc_force[i]), 0.5);
    double3 damp_avg = scal_mul(add(damping_force_old[i], particles.damping_force[i]), 0.5);
    double3 total_avg = add(add(acc_avg, lin_acc_avg), damp_avg);

    particles.vel[i] = add(particles.vel[i], scal_mul(total_avg, dt));
}

void time_step_verlet(ParticleData& particles, double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient, double smoothing_length) {
    int blocks = (N + 255) / 256; // Calculate number of blocks needed
    int threads = 256; // Number of threads per block
    // 1. Allocate device memory for acc_old
    compute_density<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_pressure<<<blocks, threads>>>(particles, N, GAMMA, K);
    cudaDeviceSynchronize();
    compute_cs<<<blocks, threads>>>(particles, N, GAMMA);
    cudaDeviceSynchronize();
    // Update acceleration
    compute_acceleration<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_linear_acceleration_force<<<blocks, threads>>>(particles, N, lambda);
    cudaDeviceSynchronize();
    compute_damping_force<<<blocks, threads>>>(particles, N, damping_coefficient);
    cudaDeviceSynchronize();
    double3* acc_old;
    double3* linear_acc_old;
    double3* damping_force_old;
    cudaMalloc(&linear_acc_old, N * sizeof(double3));
    cudaMalloc(&damping_force_old, N * sizeof(double3));
    cudaMalloc(&acc_old, N * sizeof(double3));
    // 2. Copy current acceleration to acc_old (device to device copy)
    cudaMemcpy(acc_old, particles.acc, N * sizeof(double3), cudaMemcpyDeviceToDevice);
    cudaMemcpy(linear_acc_old, particles.linear_acc_force, N * sizeof(double3), cudaMemcpyDeviceToDevice);
    cudaMemcpy(damping_force_old, particles.damping_force, N * sizeof(double3), cudaMemcpyDeviceToDevice);
    // Calculate new positions
    update_position<<<blocks, threads>>>(particles, N, dt);
    cudaDeviceSynchronize();
    compute_acceleration<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_linear_acceleration_force<<<blocks, threads>>>(particles, N, lambda);
    cudaDeviceSynchronize();
    compute_damping_force<<<blocks, threads>>>(particles, N, damping_coefficient);
    cudaDeviceSynchronize();
    // Calculate new velocities
    update_velocity<<<blocks, threads>>>(particles, N, dt, acc_old, linear_acc_old, damping_force_old);
    cudaDeviceSynchronize();
    // Free old acceleration memory
    cudaFree(acc_old);
    cudaFree(linear_acc_old);
    cudaFree(damping_force_old);
    // Compute density, pressure, sound speed, linear acceleration, and damping force
  
}