#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <tuple>

double M = 5.974e24; // mass of the Earth in kg
double m = 7.348e22; // mass of the Moon
double omega = 2.662e-6; // orbital frequency of moon and satellite
double R = 3.844e8; // distance earth to moon in m
double G = 6.67430e-11; // gravitational constant

double f_lagrange(double r, void *params) {
    return (G*M/(pow(r,2)) - G*m/(pow(R-r,2)) - pow(omega,2)*r);
}

std::pair<double, int> root_finder(){
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver * s;
    gsl_function F;

    F.function = &f_lagrange;
    F.params = nullptr;

    double r_lower = 6e6; // lower bound: circa der Erdradius
    double r_upper = R-1e6; // upper bound: kurz vor Mond
    double r;

    T = gsl_root_fsolver_brent; // solver type: Brent-Dekker method
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set(s, &F, r_lower, r_upper);

    // Iteration
    int status;
    int iter = 0;
    int max_iter = 100;

    do{
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        r_lower = gsl_root_fsolver_x_lower(s);
        r_upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(r_lower, r_upper, 0, 0.00001);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    if (status != GSL_SUCCESS){
        return std::make_pair(0.0, iter);
    }

    return std::make_pair(r, iter);
}

int main() {
    double r;
    int iterations;
    std::tie(r, iterations) = root_finder();
    if (r==0.0){
        std::cout << "No solution found after " << iterations << " iterations." << std::endl;
    } else {
        std::cout << "Found solution at r = " << r/1000 << " km after " << iterations << " iterations." << std::endl;
    }

    return 0;
}