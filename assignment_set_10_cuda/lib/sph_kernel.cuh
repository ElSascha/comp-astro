#pragma once
#include "sph_functions.cuh"


// contants for SPH simulation
constexpr double RADIUS = 0.5;
constexpr double MASS = 2.0;
constexpr double POLYTROPIC_INDEX = 1.0;
constexpr double POLYTROPIC_EXPONENT = 1.0 + 1.0 / POLYTROPIC_INDEX;
constexpr double POLYTROPIC_CONSTANT = 0.1;
constexpr double DAMPING_COEFFICIENT = 0.1;
constexpr double SMOOTHING_LENGTH = 0.2;
constexpr int NUM_PARTICLES = 1000;
constexpr double TIME = 20.0;
constexpr double TIME_STEP = 0.01;
constexpr double NORMALIZATION_CONSTANT = 8/M_PI;
constexpr double LAMBDA = 2.01203286081606;

// Use SoA (Structure of Arrays) for particle data
struct ParticleData {
    double3* pos;      // position
    double3* vel;      // velocity
    double3* acc;      // acceleration
    double3* linear_acc_force;
    double3* damping_force;      
    double* mass;
    double* rho;
    double* pressure;
    double* cs;       // sound speed
    
};

// CUDA-Kernel-Deklaration

__device__ double cubic_bspline(double r, double h); 
__global__ void compute_density(ParticleData particles, int N, double smoothing_length);
__global__ void compute_pressure(ParticleData particles, int N, double GAMMA, double K);
__global__ void compute_cs(ParticleData particles, int N, double GAMMA);
__global__ void compute_acceleration(ParticleData particles, int N, double smoothing_length);
__device__ double cubic_bspline_derivative(double r, double h);
__device__ double3 nabla_cubic_bspline(double3 pos_i, double3 pos_j, double h);
__global__ void compute_linear_acceleration_force(ParticleData particles, int N, double lambda);
__global__ void compute_damping_force(ParticleData particles, int N, double damping_coefficient);
