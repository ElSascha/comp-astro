#include <iostream>
#include <vector>
#include <cmath>
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