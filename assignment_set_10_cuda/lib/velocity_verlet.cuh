#pragma once
#include "sph_kernel.cuh"

__global__ void update_position(ParticleData particles, int N, double dt);
void time_step(ParticleData& particles, double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient, double smoothing_length);
__global__ void update_velocity(ParticleData particles, int N, double dt, double3* acc_old, double3* linear_acc_old, double3* damping_force_old);

