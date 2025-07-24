#include "sph_kernel.cuh"

__device__ double cubic_bspline(double r, double h){
    double cons = NORMALIZATION_CONSTANT/pow(h, 3);
    if(0 <= r/h && r/h < 1.0/2.0){
        return cons * (6*(pow(r/h,3)) - 6*(pow(r/h,2)) + 1);
    }
    else if(1.0/2.0 <= r/h && r/h <= 1.0){
        return cons * 2*(pow(1-r/h,3));
    }
    else{
        return 0.0;
    }
}

__device__ double3 cubic_bspline_derivative(double3 a, double3 b, double h){
    double r = length(sub(a, b));
    if (r == 0.0) return make_double3(0.0, 0.0, 0.0); // Avoid division by zero
    double3 norm_r = scal_div(sub(a,b),r);
    double cons = 6*NORMALIZATION_CONSTANT/pow(h, 4);
    if(0 <= r/h && r/h < 1.0/2.0){
        return scal_mul(norm_r, cons * (3*(pow(r/h,2)) - 2*(r/h)));
    }
    else if(1.0/2.0 <= r/h && r/h <= 1.0){
        return scal_mul(norm_r, cons * -1*(pow(1-r/h,2)));
    }
    else{
        return make_double3(0.0, 0.0, 0.0); // No contribution outside the smoothing length
    }
}

// instead of having 2 seperate loops, I will start a thread for each particle and compute the density and pressure in one go
// every thread will compute the density and pressure for one particle
// this will be more efficient than having 2 loops, as we can use shared memory to store the particles and avoid redundant memory accesses


__global__ void compute_density(ParticleData particles, int N, double smoothing_length) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Calculate the global thread index
    if (i >= N) return; // Ensure we do not access out of bounds

    particles.rho[i] = 0.0;
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

__global__ void compute_acceleration(ParticleData particles, int N, double smoothing_length) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return;

    particles.acc[i] = make_double3(0.0, 0.0, 0.0); // Initialize acceleration
    for(int j = 0; j < N ; j++){
        if(j != i){
            particles.acc[i] = add(particles.acc[i], 
                scal_mul(cubic_bspline_derivative(particles.pos[i], particles.pos[j], smoothing_length),
                 particles.mass[j] * (particles.pressure[j]/pow(particles.rho[j], 2) + particles.pressure[i]/pow(particles.rho[i], 2))));
        }
    }
    particles.acc[i] = add(add(particles.acc[i], particles.linear_acc_fore[i]), particles.damping_force[i]); // Combine acceleration with linear acceleration fore and damping force
}

__global__ void compute_pressure(ParticleData particles, int N, double GAMMA, double K) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.pressure[i] = K * pow(particles.rho[i], GAMMA); // Compute pressure using the equation pressure = pi = K*ρ_i^γ
}

__global__ void compute_linear_acceleration_fore(ParticleData particles, int N, int lambda) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.linear_acc_fore[i] = scal_mul(particles.pos[i],-1*lambda); // Compute linear acceleration fore using the equation a_i = -λ * x_i
}

__global__ void compute_damping_force(ParticleData particles, int N, double damping_coefficient) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    if (i >= N) return; 
    particles.damping_force[i] = scal_mul(particles.vel[i], -damping_coefficient); // Compute damping force using the equation f_i = -nu * v_i
}