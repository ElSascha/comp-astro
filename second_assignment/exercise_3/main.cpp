#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <cmath>



std::vector<double> quadratic_solution(double a, double b, double c){
    if(a == 0){ // linear case
        std::vector<double> result = {-c/b};
        return result;
    }

    double discriminant = b*b - 4.0*a*c;
    if(discriminant < 0){
        std::cout<<"No real solution"<<std::endl;
        return std::vector<double>();
    }
    if(discriminant == 0){
        std::vector<double> result = {-b/(2.0*a)};
        return result;
    }

    std::vector<double> result = {(-b + sqrt(discriminant))/(2.0*a), (-b - sqrt(discriminant))/(2.0*a)};
    return result;
}

std::vector<double> quadratic_solution_2(double a, double b, double c){
    if(a == 0){ // linear case
        std::vector<double> result = {-c/b};
        return result;
    }

    double discriminant = b*b - 4.0*a*c;
    if(discriminant < 0){
        std::cout<<"No real solution"<<std::endl;
        return std::vector<double>();
    }
    if(discriminant == 0){
        std::vector<double> result = {-b/(2.0*a)};
        return result;
    }

    // here is where the magic happens
    std::vector<double> result = {-(2.0*c)/(b+sqrt(discriminant)), (-b - sqrt(discriminant))/(2.0*a)};
    return result;
}

int main(){
    double a = 1.0;
    double b = 1.0;
    
    std::cout.precision(15);

    for(int exp = 1; exp < 20; exp++){
        std::vector<double> solutions = quadratic_solution(a, b, pow(10, -exp));
        std::cout << "For exp = " << exp << ": " << solutions[0] << ", " << solutions[1] << std::endl;
    }
    std::cout << std::endl;
    for(int exp = 1; exp < 20; exp++){
        std::vector<double> solutions = quadratic_solution_2(a, b, pow(10, -exp));
        std::cout << "For exp = " << exp << ": " << solutions[0] << ", " << solutions[1] << std::endl;
    }
    std::cout << "Machine epsilon for double: " << std::numeric_limits<double>::epsilon() << std::endl;
    // // Test the function with both solutions
    // std::vector<double> solutions = quadratic_solution(a, b, c);

    // std::cout << "Quadratic equation: " << a << "x² + " << b << "x + " << c << " = 0" << std::endl;
    // std::cout << "Solution 1: " << solutions[0] << std::endl;
    // std::cout << "Solution 2: " << solutions[1] << std::endl;

    return 0;
}
