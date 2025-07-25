#include "../lib/sph_kernel.cuh"
#include "../lib/euler.cuh"
#include "../lib/velocity_verlet.cuh"
#include <fstream>


//SPH toy Star with gravity using CUDA
void init_particles(ParticleData& particles, int N) {
    cudaMalloc(&particles.pos, N * sizeof(double3));
    cudaMalloc(&particles.vel, N * sizeof(double3));
    cudaMalloc(&particles.acc, N * sizeof(double3));
    cudaMalloc(&particles.linear_acc_force, N * sizeof(double3));
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
    cudaFree(particles.linear_acc_force);
    cudaFree(particles.damping_force);

    cudaFree(particles.mass);
    cudaFree(particles.rho);
    cudaFree(particles.pressure);
    cudaFree(particles.cs);
}

int main(){
    ParticleData particles;
    // Initialize particles
    init_particles(particles, NUM_PARTICLES);
    // Initialize host data
    std::vector<double3> h_pos(NUM_PARTICLES);
    std::vector<double3> h_vel(NUM_PARTICLES);
    std::vector<double3> h_acc(NUM_PARTICLES);
    std::vector<double> h_pressure(NUM_PARTICLES);
    std::vector<double> h_rho(NUM_PARTICLES);
    std::vector<double> h_cs(NUM_PARTICLES);
    std::vector<double> h_mass(NUM_PARTICLES);
    std::vector<double3> h_linear_acc_force(NUM_PARTICLES);
    std::vector<double3> h_damping_force(NUM_PARTICLES);
    // Initialize host data with distribution from file
    std::ifstream infile("initial_dis/random_distribution.dat");
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
        h_linear_acc_force[i] = {0.0, 0.0, 0.0}; // Initial linear acceleration
        h_damping_force[i] = {0.0, 0.0, 0.0}; // Initial damping force
    }
    // Copy host data to device
    cudaMemcpy(particles.pos, h_pos.data(), NUM_PARTICLES * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.vel, h_vel.data(), NUM_PARTICLES * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.acc, h_acc.data(), NUM_PARTICLES * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.pressure, h_pressure.data(), NUM_PARTICLES * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.rho, h_rho.data(), NUM_PARTICLES * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.cs, h_cs.data(), NUM_PARTICLES * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.linear_acc_force, h_linear_acc_force.data(), NUM_PARTICLES * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.damping_force, h_damping_force.data(), NUM_PARTICLES * sizeof(double3), cudaMemcpyHostToDevice);
    cudaMemcpy(particles.mass, h_mass.data(), NUM_PARTICLES * sizeof(double), cudaMemcpyHostToDevice);
    // start counter to measure time
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    int steps = TIME / TIME_STEP;
    std::ofstream outfile("data/particles.dat");
    if(!outfile) {
            std::cerr << "Error opening output file"<< "\n";
            return 1;
        }
    outfile << "time_step;x;y;z;vx;vy;vz;mass;pressure;rho;cs;linear_acc_fore_x;linear_acc_fore_y;linear_acc_fore_z;damping_force_x;damping_force_y;damping_force_z\n";
    for(int i = 0; i < steps; i++){
         time_step_verlet(particles, TIME_STEP, NUM_PARTICLES, POLYTROPIC_EXPONENT, POLYTROPIC_CONSTANT, LAMBDA, DAMPING_COEFFICIENT, SMOOTHING_LENGTH);

        // Copy results back to host for verification or further processing
        cudaMemcpy(h_pos.data(), particles.pos, NUM_PARTICLES * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_vel.data(), particles.vel, NUM_PARTICLES * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_acc.data(), particles.acc, NUM_PARTICLES * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_pressure.data(), particles.pressure, NUM_PARTICLES * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_rho.data(), particles.rho, NUM_PARTICLES * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cs.data(), particles.cs, NUM_PARTICLES * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_linear_acc_force.data(), particles.linear_acc_force, NUM_PARTICLES * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_damping_force.data(), particles.damping_force, NUM_PARTICLES * sizeof(double3), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_mass.data(), particles.mass, NUM_PARTICLES * sizeof(double), cudaMemcpyDeviceToHost);
        // write into file

        for(int j = 0; j < NUM_PARTICLES; ++j) {
            outfile << i * TIME_STEP << ";"
                    << h_pos[j].x << ";" << h_pos[j].y << ";" << h_pos[j].z << ";"
                    << h_vel[j].x << ";" << h_vel[j].y << ";" << h_vel[j].z << ";"
                    << h_mass[j] << ";"
                    << h_pressure[j] << ";"
                    << h_rho[j] << ";"
                    << h_cs[j] << ";"
                    << h_linear_acc_force[j].x << ";" << h_linear_acc_force[j].y << ";" << h_linear_acc_force[j].z << ";"
                    << h_damping_force[j].x << ";" << h_damping_force[j].y << ";" << h_damping_force[j].z << "\n";
        }
        outfile.flush();
        if (i % 100 == 0) {
            std::cout << "Step " << i << " done\n";
        }
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "Total time: " << milliseconds << " ms\n";

    // Free device memory
    free_particles(particles);
    return 0;

}