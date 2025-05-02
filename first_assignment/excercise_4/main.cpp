#include <iostream>
#include <string>
#include <cmath>
#include <iomanip> 
#include <fstream>

double find_pi_recursion_worker(int n,int k, double A_n){
    
    if (n == k){
        return A_n; // retrun A_n+1
    }
    else{
        double A_npo = pow(2.0, k) * sqrt(2.0 * (1.0 - sqrt(1.0 - pow(A_n / (pow(2.0,k)), 2.0))));
        return find_pi_recursion_worker(n,k+1,A_npo); // recurision step
    }
    
}

long double find_pi_Kahan_worker(int n, int k, long double A_n){
    long double Z_n = 2.0*pow(A_n/pow(2.0,k+1),2.0)/(1.0 + sqrt(1.0 - pow(A_n/pow(2.0,k),2.0)));
    long double A_npo = pow(2.0,k)*sqrt(4.0 * Z_n);
    if (n == k){
        return A_n;
    }
    else{
        return find_pi_Kahan_worker(n,k+1,A_npo); // recurision step
    }
}

long double find_pi_Kahan(int n){
    return find_pi_Kahan_worker(n,1,2.0); 
    
}

double find_pi_recursion(int n){
    return find_pi_recursion_worker(n,1,2); // A_0 = 2
}

int main(){
    std::ofstream pi_recursion;
    pi_recursion.open ("pi_recursion_data.txt");
    for (int i = 1; i <=30 ; i++){
        double A_n = find_pi_recursion(i);
        pi_recursion << i << ";" << std::setprecision(30)<< A_n << std::endl;
    }
    pi_recursion.close();

    std::ofstream pi_Kahan;
    pi_Kahan.open ("pi_Kahan_data.txt");
    for (int i = 1; i <=500 ; i++){
        long double A_n = find_pi_Kahan(i);
        pi_Kahan << i << ";" << std::setprecision(3000)<< A_n << std::endl;
    }
    pi_Kahan.close();
}