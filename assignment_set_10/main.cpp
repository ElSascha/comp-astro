
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

// Parameters
double n = 1.0;                 // polytropic index
double Gamma = 1.0 + 1.0 / n;   // polytropic exponent
double K = 0.1;                 // polytropic constant
double nu = 0.1;                // damping coefficient
double h = 0.2;                 // smoothing length
double N = 1000;                // number of particles
double t_end = 20.0;            // simulation time
double dt = 0.01;               // fixed time step size

const double lambda = 2.01203286081606;
double const pi = 3.141592653589793;

std::vector<double> x, y, z, vx, vy, vz, m;

void loadData() {
    // Load particle data from file
    std::ifstream file("random_distribution.dat");
    if (!file.is_open()) {
        std::cerr << "Error opening file: random_distribution.dat" << std::endl;
        return;
    }

    double temp;

    while (file >> temp) {
        x.push_back(temp);
        file >> temp;
        y.push_back(temp);
        file >> temp;
        z.push_back(temp);
        file >> temp;
        vx.push_back(temp);
        file >> temp;
        vy.push_back(temp);
        file >> temp;
        vz.push_back(temp);
        file >> temp;
        m.push_back(temp);
    }

    file.close();
}

double kernelFunction(double r_tilde) {
    if (0 <= r_tilde/h && r_tilde/h < 0.5) {
        return 8/(pi * h * h * h) * (6.0 * pow(r_tilde/h, 3) - 6.0 * pow(r_tilde/h, 2) + 1.0);
    } else if (0.5 <= r_tilde/h && r_tilde/h <= 1) {
        return 8/(pi * h * h * h) * (2.0 * pow(1.0 - r_tilde/h, 3));
    } else {
        return 0.0;
    }
}

double kernelFunctionDerivative(double r_tilde) {
    if (0 <= r_tilde/h && r_tilde/h < 0.5) {
        return 6*8/(pi * h * h * h * h) * (3 * pow(r_tilde/h, 2) - 2 * (r_tilde/h));
    } else if (0.5 <= r_tilde/h && r_tilde/h <= 1) {
        return -6 * 8/(pi * h * h * h * h) * pow(1.0 - r_tilde/h, 2);
    } else {
        return 0.0;
    }
}

std::vector<double> nablaKernelFunction(double x1, double y1, double z1, double x2, double y2, double z2) {
    std::vector<double> result(3, 0.0);
    double r_tilde = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
    double coeff = kernelFunctionDerivative(r_tilde);
    result[0] = coeff * (x1 - x2) / r_tilde;
    result[1] = coeff * (y1 - y2) / r_tilde;
    result[2] = coeff * (z1 - z2) / r_tilde;
    return result;
}

std::vector<double> computeDensity() {
    std::vector<double> density(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double r_tilde = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2) + pow(z[i] - z[j], 2));
            density[i] += m[j] * kernelFunction(r_tilde);
        }
    }
    return density;
}

std::vector<double> computeSoundSpeed(const std::vector<double>& density) {
    std::vector<double> soundSpeed(N, 0.0);
    for (int i = 0; i < N; ++i) {
        soundSpeed[i] = sqrt(K * Gamma * pow(density[i], Gamma - 1.0));
    }
    return soundSpeed;
}

std::vector<double> computePressure(const std::vector<double>& density) {
    std::vector<double> pressure(N, 0.0);
    for (int i = 0; i < N; ++i) {
        pressure[i] = K * pow(density[i], Gamma);
    }
    return pressure;
}

// this is the acceleration due to pressure
std::vector<std::vector<double>> computeAcceleration(const std::vector<double>& density, const std::vector<double>& pressure) {
    std::vector<std::vector<double>> accelerations(N, std::vector<double>(3, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::vector<double> nablaK = nablaKernelFunction(x[i], y[i], z[i], x[j], y[j], z[j]);
            accelerations[i][0] -= m[j] * (pressure[i]/pow(density[i], 2) + pressure[j]/pow(density[j], 2)) * nablaK[0];
            accelerations[i][1] -= m[j] * (pressure[i]/pow(density[i], 2) + pressure[j]/pow(density[j], 2)) * nablaK[1];
            accelerations[i][2] -= m[j] * (pressure[i]/pow(density[i], 2) + pressure[j]/pow(density[j], 2)) * nablaK[2];
        }
    }
    return accelerations;
}

std::vector<std::vector<double>> computeLinearAcceleration(){
    std::vector<std::vector<double>> linearAcceleration(N, std::vector<double>(3, 0.0));
    for (int i = 0; i < N; ++i) {
        linearAcceleration[i][0] = -lambda * x[i];
        linearAcceleration[i][1] = -lambda * y[i];
        linearAcceleration[i][2] = -lambda * z[i];
    }
    return linearAcceleration;
}

std::vector<std::vector<double>> computeDampingForce() {
    std::vector<std::vector<double>> dampingForce(N, std::vector<double>(3, 0.0));
    for (int i = 0; i < N; ++i) {
        dampingForce[i][0] = -nu * vx[i];
        dampingForce[i][1] = -nu * vy[i];
        dampingForce[i][2] = -nu * vz[i];
    }
    return dampingForce;
}

int main() {
    loadData();
    for(double t = 0.0; t < t_end; t += dt) {
        std::vector<double> density = computeDensity();
        std::vector<double> soundSpeed = computeSoundSpeed(density);
        std::vector<double> pressure = computePressure(density);
        std::vector<std::vector<double>> accelerations = computeAcceleration(density, pressure);
        std::vector<std::vector<double>> linearAcceleration = computeLinearAcceleration();
        std::vector<std::vector<double>> dampingForce = computeDampingForce();

        // Euler integration step
        for (int i = 0; i < N; ++i) {
            vx[i] += (accelerations[i][0] + linearAcceleration[i][0] + dampingForce[i][0]) * dt;
            vy[i] += (accelerations[i][1] + linearAcceleration[i][1] + dampingForce[i][1]) * dt;
            vz[i] += (accelerations[i][2] + linearAcceleration[i][2] + dampingForce[i][2]) * dt;

            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;
            z[i] += vz[i] * dt;
        }
    }
    return 0;
}