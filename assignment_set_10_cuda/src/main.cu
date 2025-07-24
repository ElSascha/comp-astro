#include "../lib/sph_kernel.cuh"
#include <fstream>

//SPH toy Star with gravity using CUDA
void init_particles(ParticleData& particles, int N) {
    cudaMalloc(&particles.pos, N * sizeof(double3));
    cudaMalloc(&particles.vel, N * sizeof(double3));
    cudaMalloc(&particles.acc, N * sizeof(double3));
    cudaMalloc(&particles.linear_acc_fore, N * sizeof(double3));
    cudaMalloc(&particles.damping_force, N * sizeof(double3));

    cudaMalloc(&particles.mass, N * sizeof(double));
    cudaMalloc(&particles.rho, N * sizeof(double));
    cudaMalloc(&particles.pressure, N * sizeof(double));
    cudaMalloc(&particles.cs, N * sizeof(double));
}

void free_particles(ParticleData& particles) {
    cudaFree(particles.pos);
    cudaFree(particles.vel);
    cudaFree(particles.acc);
    cudaFree(particles.linear_acc_fore);
    cudaFree(particles.damping_force);

    cudaFree(particles.mass);
    cudaFree(particles.rho);
    cudaFree(particles.pressure);
    cudaFree(particles.cs);
}

int main(){
    ParticleData particles;
    cudaMalloc(&particles.pos, NUM_PARTICLES * sizeof(double3));
    cudaMalloc(&particles.vel, NUM_PARTICLES * sizeof(double3));
    cudaMalloc(&particles.acc, NUM_PARTICLES * sizeof(double3));
    cudaMalloc(&particles.pressure, NUM_PARTICLES * sizeof(double));
    cudaMalloc(&particles.rho, NUM_PARTICLES * sizeof(double));
    cudaMalloc(&particles.cs, NUM_PARTICLES * sizeof(double));
    cudaMalloc(&particles.linear_acc_fore, NUM_PARTICLES * sizeof(double3));
    cudaMalloc(&particles.damping_force, NUM_PARTICLES * sizeof(double3));
    // Initialize host data
    std::vector<double3> h_pos(NUM_PARTICLES);
    std::vector<double3> h_vel(NUM_PARTICLES);
    std::vector<double3> h_acc(NUM_PARTICLES);
    std::vector<double> h_pressure(NUM_PARTICLES);
    std::vector<double> h_rho(NUM_PARTICLES);
    std::vector<double> h_cs(NUM_PARTICLES);
    std::vector<double> h_mass(NUM_PARTICLES);
    std::vector<double3> h_linear_acc_fore(NUM_PARTICLES);
    std::vector<double3> h_damping_force(NUM_PARTICLES);
    // Initialize host data with distribution from file
    std::ifstream infile("../initial_dis/random_distribution.dat");
    if(!infile) {
        std::cerr << "Error opening initial data file\n";
        return 1;
    }
    for(int i = 0; i < NUM_PARTICLES; ++i) {
        infile >> h_pos[i].x >> h_pos[i].y >> h_pos[i].z >> h_vel[i].x >> h_vel[i].y >> h_vel[i].z >> h_mass[i];
        h_acc[i] = {0.0, 0.0, 0.0}; // Initial acceleration
        h_pressure[i] = 0.0; // Initial pressure
        h_rho[i] = 1.0; // Initial density
        h_cs[i] = 1.0; // Initial sound speed
        h_linear_acc_fore[i] = {0.0, 0.0, 0.0}; // Initial linear acceleration
        h_damping_force[i] = {0.0, 0.0, 0.0}; // Initial damping force
    }


}