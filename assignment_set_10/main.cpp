
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



int main() {
    loadData();
    std::cout << "Length of vectors " << x.size() << std::endl;
    return 0;
}