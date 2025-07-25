#pragma once
#include "sph_kernel.cuh"

void time_step_euler(ParticleData& particles, double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient, double smoothing_length);
__global__ void euler_update_position(ParticleData particles, int N, double dt);