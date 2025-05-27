#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

const double M = 5.974e24; // kg
const double G = 6.67430e-11; // m^3 kg^-1 s^-2
const double R = 3.844e8; // m
const double m = 7.348e22; // kg
const double w = 2.662e-6; // 1/s

double f(double r, void *params) {
    
    return (G*M/(r*r)) - (G*m / pow(R-r,2)) - (w*w*r);
}

int main(){
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r;
    double r_lo = 6e6; // lower bound: circa der Erdradius
    double r_hi = R - 1e6; // upper bound: kurz vor Mond
    gsl_function F;
    F.function = &f;
    F.params = nullptr;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, r_lo, r_hi);

    printf ("using %s method\n",
          gsl_root_fsolver_name (s));

    do
        {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        r_lo = gsl_root_fsolver_x_lower (s);
        r_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (r_lo, r_hi,
                                       0, 0.001);
        }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    std::cout << "Root found at: " << r << " m" << std::endl; 
    std::cout << "Number of iterations: " << iter << std::endl;
    return 0;
}