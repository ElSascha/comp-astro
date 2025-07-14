#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>


using namespace std;

const double c_s = 1.0;
const double eps = 1e-8; // for numerical stability
int grid_size = 500;
double sigma = 0.5;
double domain[2] = {0, 100};
int steps = 100;
int ghost_cells = 2; // Number of ghost cells for boundary conditions

double initial_density(double x) {
    return 1.0 + 0.3 * exp(-pow(x - 50, 2) / 10);
}

void upwind_advection_step(vector<double> &rho, vector<double> &momentum, double dt, double dx) {

    static vector<double> rho_half(rho.size()); // Temporary storage for half-step density
    static vector<double> momentum_half(momentum.size()); // Temporary storage for half-step momentum density

    for (int i = ghost_cells; i < grid_size + ghost_cells; ++i) {
        double u_half_pos = 0.5 * (momentum[i] / max(rho[i], eps) + momentum[i + 1] / max(rho[i + 1], eps)); // Average velocity at "half-step" in positive direction
        double u_half_neg = 0.5 * (momentum[i - 1] / max(rho[i - 1], eps) + momentum[i] / max(rho[i], eps));// Average velocity at "half-step" in negative direction

        double F_rho_pos, F_rho_neg, F_mom_pos, F_mom_neg; // Fluxes in positive and negative directions

        // Upwind: positive direction
        if (u_half_pos > 0) {
            F_rho_pos = rho[i] * u_half_pos;
            F_mom_pos = momentum[i] * u_half_pos;
        } else {
            F_rho_pos = rho[i + 1] * u_half_pos;
            F_mom_pos = momentum[i + 1] * u_half_pos;
        }

        // Upwind: negative direction
        if (u_half_neg > 0) {
            F_rho_neg = rho[i - 1] * u_half_neg;
            F_mom_neg = momentum[i - 1] * u_half_neg;
        } else {
            F_rho_neg = rho[i] * u_half_neg;
            F_mom_neg = momentum[i] * u_half_neg;
        }

        // Update half-step
        rho_half[i] = rho[i] - dt / dx * (F_rho_pos - F_rho_neg);
        momentum_half[i] = momentum[i] - dt / dx * (F_mom_pos - F_mom_neg);
    }

    // Final update with pressure force
    for (int i = ghost_cells; i < grid_size +ghost_cells; ++i) {
        rho[i] = rho_half[i];
        momentum[i] = momentum_half[i] - dt * pow(c_s, 2) * (rho_half[i + 1] - rho_half[i - 1]) / (2 * dx);
    }
}

void reflective_boundary_conditions(vector<double> &rho, vector<double> &momentum) {
   
    // left
    for (int i = 0; i < ghost_cells; ++i) {
        rho[i] = rho[2 * ghost_cells - 1 - i];
        momentum[i] = -momentum[2 * ghost_cells - 1 - i];
    }

    // right
    for (int i = 0; i < ghost_cells; ++i) {
        rho[ghost_cells+ grid_size + i] = rho[ghost_cells + grid_size - 1 - i];
        momentum[ghost_cells + grid_size + i] = -momentum[ghost_cells + grid_size - 1 - i];
    }
}

int main() {
    double dx = (domain[1] - domain[0]) / grid_size;
    double dt = sigma * dx / c_s;
    vector<double> rho(grid_size + 2*ghost_cells);       // density with ghost cells
    vector<double> momentum(grid_size + 2*ghost_cells);  // momentum density with ghost cells

    // initial conditions
    for (int i = 0; i < grid_size ; ++i) {
        double x = domain[0] + (i - 0.5) * dx;
        rho[i + ghost_cells] = initial_density(x);
        momentum[i + ghost_cells] = 0.0;
    }

    // boundary conditions
    reflective_boundary_conditions(rho, momentum);

    // save output to file
    std::ostringstream filename;
    filename << "../isothermal_flow/output/isothermal_flow.txt";
    
    std::ofstream file(filename.str());
    if (!file.is_open()) {
        cerr << "error opening file: " << filename.str() << endl;
        return 1;
    }

    file << "Step;Position;Density;Momentum\n";

    // Time integration loop
    for (int step = 0; step < steps; ++step) {
        reflective_boundary_conditions(rho, momentum);
        upwind_advection_step(rho, momentum, dt, dx);
        // output results
        for (int i = 0; i < grid_size; ++i) {
            double x = domain[0] + (i - 0.5) * dx;
            file << step << ";"
                 << fixed << setprecision(6) << x << ";"
                 << rho[i + ghost_cells] << ";"
                 << momentum[i + ghost_cells] << "\n";
        }
        file.flush();  // safety flush to ensure data is written to file
    }

    file.close();
    cout << "simulation complete" << filename.str() << endl;
    return 0;
}
