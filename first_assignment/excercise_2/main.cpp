#include <iostream>
#include <string>
#include <cmath>
#include <limits.h>

float find_epsilon_single(){
    float epsilon = 1.0f;
    while (1.0f + epsilon/2.0f > 1.0f){
        epsilon = epsilon/2.0f;
    }
    return epsilon;
}

double find_epsilon_double(){
    double epsilon = 1.0;
    while (1.0 + epsilon/2.0 > 1.0){
        epsilon = epsilon/2.0;
    }
    return epsilon;
}

int main(){
    float epsilon = find_epsilon_single();
    std::cout << "Epsilon for single precision: " << epsilon << std::endl;
    double epsilon_double = find_epsilon_double();
    std::cout << "Epsilon for double precision: " << epsilon_double << std::endl;
    std::cout << "Epsilon Float by lib: " << std::numeric_limits<float>::epsilon() << std::endl;
    std::cout << "Epsilon double by lib: " << std::numeric_limits<double>::epsilon() << std::endl;
    
    return 0;
}