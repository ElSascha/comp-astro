#ifndef SPH_KERNEL_CUH
#define SPH_KERNEL_CUH
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
constexpr double NORMALAIZATION_CONSTANT = 8/M_PI;

// Partikel-Datenstruktur
struct Particle {
    double3 pos;      // position
    double3 vel;      // velocity
    double3 acc;      // acceleration
    double mass;
    double rho;
    double pressure;
};

// CUDA-Kernel-Deklaration

__device__ double cubic_bspline(double r, double h); 

__global__ void compute_density(Particle* particles);

__global__ void compute_pressure(Particle* particles);

#endif // SPH_KERNEL_CUH