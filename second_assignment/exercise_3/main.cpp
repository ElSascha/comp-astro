#include <iostream>
#include <vector>
#include <cmath>

double* solve_quadratic(double* param){
    if (param[1] * param[1] - 4 * param[0] * param[2] < param[1] * param[1]){
        double* res = new double[2];
        res[1] = (-param[1] - sqrt(param[1]*param[1] - 4*param[0]*param[2])) / (2*param[0]);
        res[0] = (-param[1] + sqrt(param[1]*param[1] - 4*param[0]*param[2])) / (2*param[0]);
        return res;
    }
    else{
        return nullptr;
    }
}


int main(){

    double param[] = {1, 1, 1e-6}; // Coefficients of the quadratic equation ax^2 + bx + c = 0

    double* result = solve_quadratic(param);
    std::cout << "x_1: " << result[0] << ", x_2: " << result[1] << std::endl;
    // Free the allocated memory
    delete[] result;
    return 0;

}