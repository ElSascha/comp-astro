#include "euler.cuh"

__global__ void euler_update_position(ParticleData particles, int N, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double3 acc_total = add(add(particles.acc[i], particles.linear_acc_force[i]), particles.damping_force[i]);
    particles.pos[i] = add(particles.pos[i], scal_mul(particles.vel[i], dt));
    particles.vel[i] = add(particles.vel[i], scal_mul(acc_total, dt));
}

void time_step_euler(ParticleData& particles, double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient, double smoothing_length) {
    int threads = 256;
    int blocks = (N + threads - 1) / threads;

    // Compute hydrodynamic quantities
    compute_density<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_pressure<<<blocks, threads>>>(particles, N, GAMMA, K);
    cudaDeviceSynchronize();
    compute_cs<<<blocks, threads>>>(particles, N, GAMMA);
    cudaDeviceSynchronize();

    // Compute forces
    compute_acceleration<<<blocks, threads>>>(particles, N, smoothing_length);
    cudaDeviceSynchronize();
    compute_linear_acceleration_force<<<blocks, threads>>>(particles, N, lambda);
    cudaDeviceSynchronize();
    compute_damping_force<<<blocks, threads>>>(particles, N, damping_coefficient);
    cudaDeviceSynchronize();

    // Euler update
    euler_update_position<<<blocks, threads>>>(particles, N, dt);
    cudaDeviceSynchronize();
}