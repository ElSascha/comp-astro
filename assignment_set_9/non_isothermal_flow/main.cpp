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


const double eps = 1e-8; // for numerical stability
const int grid_size = 500;
double domain[2] = {0, 100};
const int steps = 40;
const double g = 1.4; // adiabatic exponent
const int ghost_cells = 2; // Number of ghost cells for boundary conditions
const double CFL = 0.5; // Courant-Friedrichs-Lewy condition
const double dx = (domain[1] - domain[0]) / grid_size;


double initial_density(double x) {
    return (x <= 50)? 2.0 : 1.0;
}

double calc_cs(double epsilon){
    return sqrt(g*(g-1)*epsilon);
}

double calc_dt(vector<double> &density, vector<double> &momentum, vector<double> &total_energy) {
    double v_max = 0.0;
    for (int i = ghost_cells; i < grid_size + ghost_cells; ++i) {
        double u = momentum[i] / max(density[i], eps);
        double e_kin = 0.5 * u * u;
        double epsilon = (total_energy[i] / max(density[i], eps)) - e_kin;
        double cs = calc_cs(epsilon);
        v_max = max(v_max, fabs(u) + cs);
    }
    return CFL * dx / v_max;
}
double p(double density, double epsilon){
    return (g -1)*density*epsilon;
}

double eps_total(double density, double momentum, double epsilon) {
    double v = momentum / max(density, eps);
    return v * v / 2 + epsilon;
}
void upwind_advection_step(vector<double> &rho, vector<double> &momentum,vector<double> &total_energy) {
    double dt = calc_dt(rho, momentum, total_energy);
    static vector<double> rho_half(rho.size()); // Temporary storage for half-step density
    static vector<double> momentum_half(momentum.size()); // Temporary storage for half-step momentum density
    static vector<double> total_energy_half(total_energy.size()); // Temporary storage for half-step total energy

    for (int i = ghost_cells; i < grid_size + ghost_cells; ++i) {
        double u_half_pos = 0.5 * (momentum[i] / max(rho[i], eps) + momentum[i + 1] / max(rho[i + 1], eps)); // Average velocity at "half-step" in positive direction
        double u_half_neg = 0.5 * (momentum[i - 1] / max(rho[i - 1], eps) + momentum[i] / max(rho[i], eps));
        double p_half_pos = 0.5 * (p(rho[i], total_energy[i] / rho[i]) + p(rho[i + 1], total_energy[i + 1] / rho[i + 1]));
        // Average velocity at "half-step" in negative direction

        double F_rho_pos, F_rho_neg, F_mom_pos, F_mom_neg, F_energy_pos, F_energy_neg; // Fluxes in positive and negative directions

        // Upwind: positive direction
        if (u_half_pos > 0) {
            F_rho_pos = rho[i] * u_half_pos;
            F_mom_pos = momentum[i] * u_half_pos;
            F_energy_pos = total_energy[i] * u_half_pos + p_half_pos * u_half_pos;
        } else {
            F_rho_pos = rho[i + 1] * u_half_pos;
            F_mom_pos = momentum[i + 1] * u_half_pos;
            F_energy_pos = total_energy[i + 1] * u_half_pos + p_half_pos * u_half_pos;
        }

        // Upwind: negative direction
        if (u_half_neg > 0) {
            F_rho_neg = rho[i - 1] * u_half_neg;
            F_mom_neg = momentum[i - 1] * u_half_neg;
            F_energy_neg = total_energy[i - 1] * u_half_neg + p(rho[i - 1], total_energy[i - 1] / rho[i - 1]) * u_half_neg;
        } else {
            F_rho_neg = rho[i] * u_half_neg;
            F_mom_neg = momentum[i] * u_half_neg;
            F_energy_neg = total_energy[i] * u_half_neg + p(rho[i], total_energy[i] / rho[i]) * u_half_neg;
        }

        // Update half-step
        rho_half[i] = rho[i] - dt / dx * (F_rho_pos - F_rho_neg);
        momentum_half[i] = momentum[i] - dt / dx * (F_mom_pos - F_mom_neg);
        total_energy_half[i] = total_energy[i] - dt / dx * (F_energy_pos - F_energy_neg);
    }

    // Final update with pressure force
    for (int i = ghost_cells; i < grid_size +ghost_cells; ++i) {
        double u_p = momentum_half[i+1] / max(rho_half[i+1], eps);
        double eps_tot_p = total_energy_half[i+1] / max(rho_half[i+1], eps);
        double eps_kin_p = 0.5 * u_p * u_p;
        double epsilon_p = eps_tot_p - eps_kin_p;
        double pressure_p = p(rho_half[i+1], epsilon_p);
        double u_m = momentum_half[i-1] / max(rho_half[i-1], eps);
        double eps_tot_m = total_energy_half[i-1] / max(rho_half[i-1], eps);
        double eps_kin_m = 0.5 * u_m * u_m;
        double epsilon_m = eps_tot_m - eps_kin_m;
        double pressure_m = p(rho_half[i-1], epsilon_m);

        rho[i] = rho_half[i];
        momentum[i] = momentum_half[i] - (dt / (2*dx))*(pressure_p - pressure_m);
        total_energy[i] = total_energy_half[i] - (dt/(2*dx)) * (pressure_p * u_p - pressure_m * u_m);
    }
}

void reflective_boundary_conditions(vector<double> &rho, vector<double> &momentum, vector<double> &total_energy) {
   
    // left
    for (int i = 0; i < ghost_cells; ++i) {
        rho[i] = rho[2 * ghost_cells - 1 - i];
        momentum[i] = -momentum[2 * ghost_cells - 1 - i];
        total_energy[i] = total_energy[2 * ghost_cells - 1 - i];
    }

    // right
    for (int i = 0; i < ghost_cells; ++i) {
        rho[ghost_cells+ grid_size + i] = rho[ghost_cells + grid_size - 1 - i];
        momentum[ghost_cells + grid_size + i] = -momentum[ghost_cells + grid_size - 1 - i];
        total_energy[ghost_cells + grid_size + i] = total_energy[ghost_cells + grid_size - 1 - i];
    }
}

int main() {
    double initial_eps = 1.0; // Initial internal energy
    vector<double> rho(grid_size + 2*ghost_cells);       // density with ghost cells
    vector<double> momentum(grid_size + 2*ghost_cells);  // momentum density with ghost cells
    vector<double> total_energy(grid_size + 2*ghost_cells); // total energy with ghost cells

    // initial conditions
    for (int i = 0; i < grid_size ; ++i) {
        double x = domain[0] + (i + 0.5) * dx;
        rho[i + ghost_cells] = initial_density(x);
        momentum[i + ghost_cells] = 0.0;
        total_energy[i + ghost_cells] = eps_total(initial_density(x), 0.0 ,initial_eps) * initial_density(x);
    }

    // boundary conditions
    reflective_boundary_conditions(rho, momentum, total_energy);

    // save output to file
    std::ostringstream filename;
    filename << "../non_isothermal_flow/output/non_isothermal_flow.txt";
    std::ofstream file(filename.str());
    if (!file.is_open()) {
        cerr << "error opening file: " << filename.str() << endl;
        return 1;
    }
    file << "Step;Position;Density;Momentum;TotalEnergy\n";

    // Time integration loop
    for (int step = 0; step < steps; ++step) {
        reflective_boundary_conditions(rho, momentum, total_energy);
        upwind_advection_step(rho, momentum, total_energy);
        // output results
        for (int i = 0; i < grid_size; ++i) {
            double x = domain[0] + (i + 0.5) * dx;
            file << step << ";"
                 << fixed << setprecision(6) << x << ";"
                 << rho[i + ghost_cells] << ";"
                 << momentum[i + ghost_cells] << ";"
                 << total_energy[i + ghost_cells] << "\n";

        }
        file.flush();  // safety flush to ensure data is written to file
    }

    file.close();
    cout << "simulation complete" << filename.str() << endl;
    return 0;
}
