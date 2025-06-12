#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using std::cout;
using std::endl;
using std::vector;
using std::sqrt;

const double mu2 = 1e-3;
const double mu1 = 1 - mu2; // given normalization mu1 + mu2 = 1

const double Omega = 1.0; // angular velocity in dimensionless units is equal to 1
const double dt = 1e-3; // time step 
const int Nsteps = 5e6; // number of steps


// this outputs its results into r1 and r2, which need to be passed as parameters
void compute_r(double x, double y, double& r1, double& r2) {
    r2 = sqrt(std::pow(x - mu1, 2) + y*y);
    r1 = sqrt(std::pow(x + mu2, 2) + y*y);
}

double compute_initial_vy(double x0, double C_J){
    double y0 = 0.0; // initial y position is assumed to be 0.0
    double r1, r2;
    compute_r(x0, y0, r1, r2);
    double vy2 = Omega * Omega * (x0*x0) + 2 * ( (mu1 / r1) + (mu2 / r2) ) - C_J; // (11) rearranged for vy². vx0 = y0 = 0.0 assumed.
    return (vy2 > 0) ? sqrt(vy2) : 0.0; // return the initial y velocity (protected against negative discriminant)
}

// computes the derivatives of the state vector, returns them into existing vec[4] dstate
void derivatives(const vector<double>& state, vector<double>& dstate) {
    double x = state[0];
    double y = state[1];
    double vx = state[2];
    double vy = state[3];

    double r1, r2;
    compute_r(x, y, r1, r2);
    double r1_3 = r1 * r1 * r1;
    double r2_3 = r2 * r2 * r2;

    dstate[0] = vx; // dx/dt = vx
    dstate[1] = vy; // dy/dt = vy
    dstate[2] = 2 * Omega * vy + Omega * Omega * x - mu1 * (x + mu2) / r1_3 - mu2 * (x - mu1) / r2_3; // dvx/dt
    dstate[3] = -2 * Omega * vx + Omega * Omega * y - y * (mu1 / r1_3 + mu2 / r2_3); // dvy/dt
}

void rk4_step(vector<double>& state){
    vector<double> k1(4), k2(4), k3(4), k4(4), temp(4);

    derivatives(state, k1);
    for (int i = 0; i < 4; ++i) 
        temp[i] = state[i] + 0.5 * dt * k1[i];

    derivatives(temp, k2);
    for (int i = 0; i < 4; ++i) 
        temp[i] = state[i] + 0.5 * dt * k2[i];

    derivatives(temp, k3);
    for (int i = 0; i < 4; ++i) 
        temp[i] = state[i] + dt * k3[i];

    derivatives(temp, k4);
    for (int i = 0; i < 4; ++i) 
        state[i] += (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}




int main(){
    double CJ = 3.03; // Jacobi constant value
    vector<double> x0s = {0.21, 0.24, 0.26, 0.27, 0.4, 0.5, 0.6, 0.8};
    cout.precision(15);

    for (double x0 : x0s) {
        vector<double> state = {x0, 0.0, 0.0, compute_initial_vy(x0, CJ)}; // initial state vector [x, y, vx, vy]
        std::ofstream out("orbit_x0_" + std::to_string(x0) + ".csv");

        for (int step = 0; step < Nsteps; ++step) {
            if (step % 500 == 0) { // output every n steps
                cout << "Calc for x0 = " << x0 << ", step = " << step << endl;
                out << state[0] << "," << state[1] << endl; // output current position
            }
            rk4_step(state); // perform RK4 step
        }
        out.close();
        cout << "Orbit for x0 = " << x0 << " computed and saved to orbit_x0_" << x0 << ".dat" << endl;
    }
    return 0;
}