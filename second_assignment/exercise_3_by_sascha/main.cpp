#include <iostream>
#include <vector>
#include <cmath>
#include <limits.h>
/*

*/
/**
 * @brief Solves a quadratic equation of the form ax^2 + bx + c = 0.
 *
 * @param param Pointer to an array of three doubles representing the coefficients [a, b, c].
 * @return Pointer to an array of doubles containing the real roots of the equation:
 *         - If the equation is linear (a == 0, b != 0), returns a pointer to an array of size 1 with the single root.
 *         - If the equation has two real roots, returns a pointer to an array of size 2 with the roots.
 *         - If there are no real roots, returns nullptr.
 * 
 * @note The returned array is dynamically allocated and must be deleted by the caller to avoid memory leaks.
 */

double* solve_quadratic(double* param){
    //Linear case
    if (param[0] == 0 && param[1] != 0){
        double* res = new double[1];
        res[0] = -param[2] / param[1];
        return res;
    }
    double a = param[0];
    double b = param[1];
    double c = param[2];
    double* res = new double[2];
    if(b*b - 4.0*a*c < b*b){
        res[1] = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
        res[0] = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
        return res;
    }
    else{
        res[0] = -(2.0*c/(b + sqrt(b*b - 4.0*a*c)));
        res[1] = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
        return res;
    }
}


int main(){

    double* param = new double[3];
    double* result;
    param[0] = 1; 
    param[1]= 1;
    std::cout.precision(20);
    std::cout<<"a = 1, b = 1, c = 1e-x"<<std::endl;
    for (int i = 1 ; i< 30; i++){
        std::cout << "c = " << std::pow(10.0, -i) << std::endl;
        param[2] = std::pow(10.0, -i);
        result = solve_quadratic(param);
        std::cout << "x_1: " << result[0] << ", x_2: " << result[1] << std::endl;
        // Free the allocated memory
        delete[] result;

    }
    return 0;

}