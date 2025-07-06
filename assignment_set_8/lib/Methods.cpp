#include "Methods.hpp"

namespace Methods {
    //  1D linear advection equation
    // using vectors instead of arrays for easier memory management

    vector<double> centered_differencing(const vector<double> &u, double sigma){
        int N = u.size() -2; // assuming u has ghost cells
        vector<double> result(u.size());
        for(int i = 1; i <= N; i++){
            result[i] = u[i] - (0.5 * sigma * (u[i+1] - u[i-1]));
        }
        return result; // return the updated array
    }

    vector<double> upwind(const vector<double> &u, double sigma){
        int N = u.size() -2; // assuming u has ghost cells
        vector<double> result(u.size());
        for(int i = 1; i <= N; i++){
            result[i] = u[i] - (sigma * (u[i] - u[i-1]));
        }
        return result; // return the updated array
    }

    vector<double> lax_wendroff(const vector<double> &u, double sigma){
        int N = u.size() -2; // assuming u has ghost cells
        vector<double> result(u.size());
        for(int i = 1; i <= N; i++){
            result[i] = u[i] - (0.5 * sigma * (u[i+1] - u[i-1])) + (sigma * sigma * 0.5 * (u[i+1] - 2*u[i] + u[i-1]));
        }
        return result; // return the updated array
    }
}