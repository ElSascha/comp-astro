#include "lib/Methods.hpp"
#include <iomanip>
#include <sstream>
using namespace Methods;



double u_start(double x){
    if(abs(x) < 1.0/3.0){
        return 1.0;
    } else {
        return 0.0;
    }
}

void periodic_boundary_conditions(vector<double> &u){
    int N = u.size() - 2; // assuming u has ghost cells
    u[0] = u[N]; // left boundary condition (ghost cell)
    u[N + 1] = u[1]; // right boundary condition (ghost cell
}

void save_state(const vector<double> &u, const std::string &filename){
    std::ofstream file(filename);
    if(file.is_open()){
        for(int i = 1; i < (int)u.size() - 1; i++){
            file << u[i] << ";"; // write to file, excluding ghost cells
        }
        file.close();
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
}

void time_step(vector<double> &u, double sigma, std::vector<double> (*method)(const std::vector<double>&, double)){
    u = method(u, sigma);
    periodic_boundary_conditions(u);
}

int main(){
    double boundries[2] = {-1.0, 1.0}; // x_min, x_max
    double a = 1.0; // advection speed
    double sigma = 0.8; // Courant number
    int grid_size = 400; // number of grid points
    double dx = (boundries[1] - boundries[0]) / grid_size; // grid spacing
    double dt = sigma * dx / a; // time step size
    double t_max = 4.0; // maximum time
    int steps = static_cast<int>(t_max / dt);// number of time steps

    vector<double> u(grid_size +2); // grid with ghost cells
    // Initialize grid
    for(int i = 0; i<grid_size; i++){
        u[i] = u_start(boundries[0] + i*dx);
    }
    // boundary conditions
    periodic_boundary_conditions(u);
    
    std::ostringstream filename;
    filename << "lax_wendroff" << "_sigma_" << std::fixed << std::setprecision(1) << sigma << ".txt";
    std::ofstream file(filename.str());
    if(!file.is_open()) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
        return 1;
    }

    // Time loop for advection equation
    for(int step = 0; step < steps; step++){
        std::cout << "Time step: " << step << std::endl;
        // Save state at each time step
        for(int i = 1; i < (int)u.size() - 1; i++){
            file << u[i];
            if(i < (int)u.size() - 2) file << ";";
        }
        file << "\n"; // new line for each time step
        time_step(u, sigma, lax_wendroff); // Change to upwind or lax_wendroff as needed
    }
    return 0;

}