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

double trapezoidal_rule(double (*f)(double), double a, double b, int n) {
    double h = (b-a)/n;
    double T = 0;
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += (f(a + i*h) + f(a * (i+1)*h));
    }
    T = (h/2) * sum;
    return T;
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

int main(){
    double a = 0.0;
    double b = 1.0;
    int n = 1000000; // Number of subintervals
    //double result = trapezoidal_rule(Si, a, b, n);
    double result = simpsons_rule(Si, a, b, n); // Uncomment to use Simpson's rule instead
    cout.precision(10); // Set precision for better output
    cout << fixed; // Use fixed-point notation
    cout << "The integral of Si(t) from " << a << " to " << b << " is approximately: " << result << endl;
    return 0;
}