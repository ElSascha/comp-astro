#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double Si(double t){
    if(t == 0) {
        return 1.0; // Si(0) is defined as 1
    }
    return sin(t)/t;
}

double C(double t){
    return cos((M_PI / 2 ) * t*t);
}

double trapezoidal_rule_worker(double (*f)(double), double a, double b, int n, double T_old) {
    double h = (b-a)/n;
    double T = 0;
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += (f(a + i*h) + f(a + (i+1)*h));
    }
    T = (h/2) * sum;
    if(n>1&&abs(T - T_old) < 1e-8) { // If the difference is small enough, return the result
        return T;
    }
    return trapezoidal_rule_worker(f,a,b,n*2,T); // Otherwise, double the number of intervals and call the function recursively
}

double trapezoidal_rule(double (*f)(double), double a, double b, int n) {
    if(n < 1) {
        n = 1; // Ensure at least one interval
    }
    double T_old = 0;
    return trapezoidal_rule_worker(f, a, b, n, T_old); // Start the recursive process
}

double simpsons_rule(double (*f)(double), double a, double b, int n) {
    if(n % 2 != 0) {
        n++; // n must be even for Simpson's rule
    }
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    double S = 0;
    for(int i = 1; i<n ; i++){
        double x = a + i * h;
        if(i % 2 == 0) {
            sum += 2 * f(x); // even index
        } else {
            sum += 4 * f(x); // odd index
        }
    }
    S = (h / 3) * sum;
    return S;
}

double adaptive_simpsons_rule_worker(double (*f)(double), double a, double b, double epsilon, 
                                    double S, double fa, double fb, double fc, int steps){
    double m = (a + b) / 2 , h = (b-a)/2; // middle point and half interval width
    double lm = (a + m)/2, rm = (m+b)/2; //left middle and right middle    
    double flm = f(lm), frm = f(rm); // function values at left and right middle points     
    double S_left = (h/6) * (fa + 4 * flm + fc), S_right = (h/6) * (fc + 4 * frm + fb); // left and right Simpson's rule
    double delta = S_left + S_right - S; // difference between the two approximations
    if(fabs(delta) < epsilon *15 || steps > 1000) { // if the difference is small enough or too many steps
        return S_left + S_right + (delta) / 15; // return the sum of left and right approximations with a correction
    } else {
        // Recursive call for left and right intervals
        return adaptive_simpsons_rule_worker(f, a, m, epsilon / 2, S_left, fa, fc,flm, steps + 1) +
               adaptive_simpsons_rule_worker(f, m, b, epsilon / 2, S_right, fc, fb,frm, steps + 1); // return the sum of the two recursive calls
    }
                                    }

double adaptive_simpsons_rule(double (*f)(double), double a, double b, double epsilon) {
    if(a == b) return 0; // If the interval is zero, return 0
    double h = (b - a) / 2;
    double fa = f(a), fb = f(b), fc = f((a + b) / 2); // function values at the endpoints and middle point
    double S = (h / 3) * (fa + 4 * fc + fb); // initial Simpson's rule approximation
    return adaptive_simpsons_rule_worker(f, a, b, epsilon, S, fa, fb, fc, 0); // start the recursive process

}

int main(){
    double a = 0.0;
    double b = 1.0;
    double c = 5.0;
    int n = 1000000; // Number of subintervals
    double epsilon = 1e-8;
    std::cout.precision(15); 
    double Si_result_trapezoidal = trapezoidal_rule(Si, a, b, n); 
    double Si_result_simpsons = simpsons_rule(Si, a, b, n);
    double Si_result_adaptive = adaptive_simpsons_rule(Si, a, b, epsilon);
    std::cout << "Si(1) using Trapezoidal Rule: " << Si_result_trapezoidal << std::endl;
    std::cout << "Si(1) using Simpson's Rule: " << Si_result_simpsons << std::endl;
    std::cout << "Si(1) using Adaptive Simpson's Rule: " << Si_result_adaptive << std::endl;
    double C_result_trapezoidal = trapezoidal_rule(C, a, c, n);
    double C_result_simpsons = simpsons_rule(C, a, c, n);
    double C_result_adaptive = adaptive_simpsons_rule(C, a, c, epsilon);
    std::cout << "C(5) using Trapezoidal Rule: " << C_result_trapezoidal << std::endl;
    std::cout << "C(5) using Simpson's Rule: " << C_result_simpsons << std::endl;
    std::cout << "C(5) using Adaptive Simpson's Rule: " << C_result_adaptive << std::endl;
    return 0;
}