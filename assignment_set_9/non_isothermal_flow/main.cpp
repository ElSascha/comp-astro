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


const int grid_size = 500;
double domain[2] = {0, 100};
const int total_time = 40;
const double g = 1.4; // adiabatic exponent
const int ghost_cells = 10; // Number of ghost cells for boundary conditions
const double CFL = 0.1; // Courant-Friedrichs-Lewy condition
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
        double u = momentum[i] / density[i];
        double e_kin = 0.5 * u * u;
        double epsilon = (total_energy[i] / density[i]) - e_kin;
        double cs = calc_cs(epsilon);
        v_max = max(v_max, fabs(u) + cs);
    }
    return CFL * dx / v_max;
}
double p(double density, double epsilon){
    return (g -1)*density*epsilon;
}

double eps_total(double density, double momentum, double epsilon) {
    double v = momentum / density;
    return v * v / 2 + epsilon;
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

void upwind_advection_step(vector<double> &density, vector<double> &momentum,vector<double> &total_energy) {
    double dt = calc_dt(density, momentum, total_energy);
    static vector<double> density_half(density.size()); // Temporary storage for half-step density
    static vector<double> momentum_half(momentum.size()); // Temporary storage for half-step momentum density
    static vector<double> total_energy_half(total_energy.size()); // Temporary storage for half-step total energy

    for (int i = ghost_cells; i < grid_size + ghost_cells; ++i) {
        double u_half_pos = 0.5 * (momentum[i] / density[i] + momentum[i + 1] / density[i+1]); // Average velocity at "half-step" in positive direction
        double u_half_neg = 0.5 * (momentum[i - 1] / density[i-1] + momentum[i] / density[i]);
        double u = momentum[i] / density[i]; // Current velocity
        double eps_tot = total_energy[i] / density[i]; // Current total energy
        double eps_kin = 0.5 * u * u; // Kinetic energy
        double epsilon = eps_tot - eps_kin; // Internal energy
        double p_half_pos = p(density[i + 1], epsilon); // Pressure at "half-step" in positive direction
        double p_half_neg = p(density[i - 1], epsilon); // Pressure at "half-step" in negative direction

        double F_rho_pos, F_rho_neg, F_mom_pos, F_mom_neg, F_energy_pos, F_energy_neg; // Fluxes in positive and negative directions

        // Upwind: positive direction
        if (u_half_pos > 0) {
            F_rho_pos = density[i] * u_half_pos;
            F_mom_pos = momentum[i] * u_half_pos;
            F_energy_pos = (total_energy[i]+p_half_pos)* u_half_pos;
        } else {
            F_rho_pos = density[i + 1] * u_half_pos;
            F_mom_pos = momentum[i + 1] * u_half_pos;
            F_energy_pos = (total_energy[i+1]+p_half_pos) * u_half_pos;
        }

        // Upwind: negative direction
        if (u_half_neg > 0) {
            F_rho_neg = density[i - 1] * u_half_neg;
            F_mom_neg = momentum[i - 1] * u_half_neg;
            F_energy_neg = (total_energy[i-1]+p_half_neg)* u_half_neg;
        } else {
            F_rho_neg = density[i] * u_half_neg;
            F_mom_neg = momentum[i] * u_half_neg;
            F_energy_neg = (total_energy[i]+p_half_neg) * u_half_neg;
        }

        // Update half-step
        density_half[i] = density[i] - dt / dx * (F_rho_pos - F_rho_neg);
        momentum_half[i] = momentum[i] - dt / dx * (F_mom_pos - F_mom_neg);
        total_energy_half[i] = total_energy[i] - dt / dx * (F_energy_pos - F_energy_neg);
    }
    
    reflective_boundary_conditions(density_half, momentum_half, total_energy_half);
    // Final update with pressure force
    for (int i = ghost_cells; i < grid_size +ghost_cells; ++i) {
        double u_p = momentum_half[i+1] / density_half[i+1];
        double u_m = momentum_half[i-1] / density_half[i-1];

        double eps_tot_p = total_energy_half[i+1] /density_half[i+1];
        double eps_kin_p = 0.5 * u_p * u_p;
        double epsilon_p = eps_tot_p - eps_kin_p;

        double eps_tot_m = total_energy_half[i-1] /density_half[i-1];
        double eps_kin_m = 0.5 * u_m * u_m;
        double epsilon_m = eps_tot_m - eps_kin_m;

        double pressure_p = p(density_half[i+1], epsilon_p);
        double pressure_m = p(density_half[i-1], epsilon_m);
        // Update final step
        density[i] = density_half[i];
        momentum[i] = momentum_half[i] - (dt / (2*dx))*(pressure_p - pressure_m);
        total_energy[i] = total_energy_half[i] - (dt/(2*dx)) * (pressure_p * u_p - pressure_m * u_m);
    }
}


int main() {
    double initial_eps = 1.0; // Initial internal energy
    vector<double> density(grid_size + 2*ghost_cells);       // density with ghost cells
    vector<double> momentum(grid_size + 2*ghost_cells);  // momentum density with ghost cells
    vector<double> total_energy(grid_size + 2*ghost_cells); // total energy with ghost cells

    // initial conditions
    for (int i = 0; i < grid_size ; ++i) {
        double x = domain[0] + (i + 0.5) * dx;
        density[i + ghost_cells] = initial_density(x);
        momentum[i + ghost_cells] = 0.0;
        total_energy[i + ghost_cells] = eps_total(initial_density(x), 0.0, initial_eps) * initial_density(x);
    }

    // boundary conditions
    reflective_boundary_conditions(density, momentum, total_energy);

    // save output to file
    std::ostringstream filename;
    filename << "../non_isothermal_flow/output/non_isothermal_flow.txt";
    std::ofstream file(filename.str());
    if (!file.is_open()) {
        cerr << "error opening file: " << filename.str() << endl;
        return 1;
    }
    file << "Step;Position;Density;Momentum;InternalEnergy\n";
    double steps_dt = 0;
    int total = 0;
    while(steps_dt < total_time) {
        double dt = calc_dt(density, momentum, total_energy);
        steps_dt += dt;
        total++;
        reflective_boundary_conditions(density, momentum, total_energy);
        upwind_advection_step(density, momentum, total_energy);

        // output results
        for (int i = 0; i < grid_size; ++i) {
            double x = domain[0] + (i + 0.5) * dx;
            file << fixed << setprecision(6) << steps_dt << ";"
                 << x << ";"
                 << density[i + ghost_cells] << ";"
                 << momentum[i + ghost_cells] << ";"
                 << total_energy[i + ghost_cells] / density[i + ghost_cells] - 0.5 * pow(momentum[i + ghost_cells] / density[i + ghost_cells], 2) << "\n";
        }
        file.flush();  // safety flush to ensure data is written to file
    }
    
    file.close();
    std::cout << total;
    return 0;
}
