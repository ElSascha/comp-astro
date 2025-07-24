#include "sph_kernel.cuh"

__global__ void update_position(ParticleData particles, int N, double dt);
void time_step(ParticleData& particles,double dt, int N, double GAMMA, double K, double lambda, double damping_coefficient);
__global__ void update_velocity(ParticleData particles, int N, double dt, double3* acc_old);

