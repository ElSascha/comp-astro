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
const int grid_size = 500;
const double sigma = 0.5;
const double domain[2] = {0, 100};
const int total_time = 100;
int ghost_cells = 6; // Number of ghost cells for boundary conditions on each side of the grid

double initial_density(double x) {
    return 1.0 + 0.3 * exp(-pow(x - 50, 2) / 10);
}

void reflective_boundary_conditions(vector<double> &density, vector<double> &momentum) {
 // left
    for (int i = 0; i < ghost_cells; ++i) {
        density[i] = density[2 * ghost_cells - 1 - i];
        momentum[i] = -momentum[2 * ghost_cells - 1 - i];
    }

    // right
    for (int i = 0; i < ghost_cells; ++i) {
        density[ghost_cells+ grid_size + i] = density[ghost_cells + grid_size - 1 - i];
        momentum[ghost_cells + grid_size + i] = -momentum[ghost_cells + grid_size - 1 - i];
    } 
}

void upwind_advection_step(vector<double> &density, vector<double> &momentum, double dt, double dx) {

    static vector<double> density_half(density.size()); // Temporary storage for half-step density
    static vector<double> momentum_half(momentum.size()); // Temporary storage for half-step momentum density

    for (int i = ghost_cells; i < grid_size + ghost_cells; ++i) {
        double u_half_pos = 0.5 * (momentum[i] / density[i] + momentum[i + 1] /density[i + 1]); // Average velocity at "half-step" in positive direction
        double u_half_neg = 0.5 * (momentum[i - 1] / density[i-1] + momentum[i] / density[i]);// Average velocity at "half-step" in negative direction

        double F_density_pos, F_density_neg, F_mom_pos, F_mom_neg; // Fluxes in positive and negative directions

        // Upwind: positive direction
        if (u_half_pos > 0) {
            F_density_pos = density[i] * u_half_pos;
            F_mom_pos = momentum[i] * u_half_pos;
        } else {
            F_density_pos = density[i + 1] * u_half_pos;
            F_mom_pos = momentum[i + 1] * u_half_pos;
        }

        // Upwind: negative direction
        if (u_half_neg > 0) {
            F_density_neg = density[i - 1] * u_half_neg;
            F_mom_neg = momentum[i - 1] * u_half_neg;
        } else {
            F_density_neg = density[i] * u_half_neg;
            F_mom_neg = momentum[i] * u_half_neg;
        }

        // Update half-step
        density_half[i] = density[i] - dt / dx * (F_density_pos - F_density_neg);
        momentum_half[i] = momentum[i] - dt / dx * (F_mom_pos - F_mom_neg);
    }
    reflective_boundary_conditions(density_half, momentum_half);
    // Final update with pressure force
    for (int i = ghost_cells; i < grid_size +ghost_cells; ++i) {
        density[i] = density_half[i];
        momentum[i] = momentum_half[i] - dt * pow(c_s, 2) * (density_half[i + 1] - density_half[i - 1]) / (2 * dx);
    }
}



int main() {
    double dx = (domain[1] - domain[0]) / grid_size;
    double dt = sigma * dx / c_s;
    vector<double> density(grid_size + 2*ghost_cells);       // density with ghost cells
    vector<double> momentum(grid_size + 2*ghost_cells);  // momentum density with ghost cells

    // initial conditions
    for (int i = 0; i < grid_size + ghost_cells ; ++i) {
        double x = domain[0] + (i + 0.5) * dx;
        density[i + ghost_cells] = initial_density(x);
        momentum[i + ghost_cells] = 0.0;
    }

    // boundary conditions
    reflective_boundary_conditions(density, momentum);

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
    double current_time = 0.0;
    int steps = 0.0;
    for (int step = 0; current_time < total_time; ++step) {
        reflective_boundary_conditions(density, momentum);
        upwind_advection_step(density, momentum, dt, dx);
        // output results
        for (int i = 0; i < grid_size; ++i) {
            double x = domain[0] + (i + 0.5) * dx;
            file << current_time << ";"
                 << fixed << setprecision(6) << x << ";"
                 << density[i + ghost_cells] << ";"
                 << momentum[i + ghost_cells] << "\n";
        }
        file.flush();  // safety flush to ensure data is written to file
        current_time += dt;
        steps++;
    }

    file.close();
    cout << "simulation complete" << filename.str() << endl;
    cout << "Total steps: " << steps << endl;   
    return 0;
}
