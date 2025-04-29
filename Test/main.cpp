#include <iostream>
#include <string>

double dfunc(double x, double y){
    return x*x + 2*y;
}

int main(){
    double a = 0.2001;
    double b = 0.1001;
    std::cout<<a+b<<std::endl;
    return 0;   
}
