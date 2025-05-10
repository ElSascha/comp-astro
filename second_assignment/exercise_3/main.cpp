#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>



/**
 * @brief Solves a quadratic equation of the form ax^2 + bx + c = 0.
 * 
 * This function computes one of the two possible solutions for a quadratic equation
 * based on the value of the solution_selector parameter. If the coefficient 'a' is zero,
 * the equation is treated as linear (bx + c = 0). If the discriminant is negative, 
 * indicating no real solutions, the function returns NaN and prints a message.
 * 
 * @param a Coefficient of the quadratic term (x^2). If a == 0, the equation is treated as linear.
 * @param b Coefficient of the linear term (x).
 * @param c Constant term.
 * @param solution_selector Determines which solution to return:
 *        - 0: Returns the solution using the positive square root of the discriminant.
 *        - 1: Returns the solution using the negative square root of the discriminant.
 * 
 * @return The selected solution to the quadratic equation. If no real solution exists,
 *         returns std::numeric_limits<double>::quiet_NaN().
 * 
 * @note If a == 0 and b == 0, the behavior is undefined.
 * @note If solution_selector is not 0 or 1, the behavior is undefined.
 */
double quadratic_solution(double a, double b, double c, int solution_selector){
    if(a == 0){ // linear case
        return -c/b;
    }

    double discriminant = b*b - 4*a*c;
    if(discriminant < 0){
        std::cout<<"No real solution"<<std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if(solution_selector == 0){
        return (-b + sqrt(discriminant))/(2*a);
    }else if (solution_selector == 1){
        return (-b - sqrt(discriminant))/(2*a);
    }
    else{
        throw std::invalid_argument("solution_selector must be 0 or 1");
    }
}

int main(){
    double a = 1.0;
    double b = -3.0;
    double c = 2.0;

    // Test the function with both solutions
    double solution1 = quadratic_solution(a, b, c, 0);
    double solution2 = quadratic_solution(a, b, c, 1);

    std::cout << "Quadratic equation: " << a << "x^2 + " << b << "x + " << c << " = 0" << std::endl;
    std::cout << "Solution 1: " << solution1 << std::endl;
    std::cout << "Solution 2: " << solution2 << std::endl;

    return 0;
}
