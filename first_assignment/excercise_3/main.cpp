#include <iostream>
#include <string>
#include <cmath>
#include <limits.h>


int main(){
    float a = 1e30;
    float b = 1;
    float c = sqrt(a) * sqrt(a + b*b/a);// thats the trick i guess
    std::cout<<std::numeric_limits<float>::max()<<std::endl;
    std::cout<<c<<std::endl;
    return 0;

}