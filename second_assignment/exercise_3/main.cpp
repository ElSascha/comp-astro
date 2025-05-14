#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <vector>



std::vector<double> quadratic_solution(double a, double b, double c){
    if(a == 0){ // linear case
        std::vector<double> result = {-c/b};
        return result;
    }

    double discriminant = b*b - 4*a*c;
    if(discriminant < 0){
        std::cout<<"No real solution"<<std::endl;
        return std::vector<double>();
    }
    if(discriminant == 0){
        std::vector<double> result = {-b/(2*a)};
        return result;
    }

    std::vector<double> result = {(-b + sqrt(discriminant))/(2*a), (-b - sqrt(discriminant))/(2*a)};
    return result;
}

int main(){
    double a = 1.0;
    double b = -3.0;
    double c = 2.0;

    // Test the function with both solutions
    std::vector<double> solutions = quadratic_solution(a, b, c);

    std::cout << "Quadratic equation: " << a << "x² + " << b << "x + " << c << " = 0" << std::endl;
    std::cout << "Solution 1: " << solutions[0] << std::endl;
    std::cout << "Solution 2: " << solutions[1] << std::endl;

    return 0;
}
