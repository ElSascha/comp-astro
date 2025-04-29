#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <limits>
using namespace std::chrono;


/*
    @brief function to calculate the harmonic sum for alpha = 2
    @param k number of summations 
    @return the harmonic sum
*/  

float general_harmonic_sum_single(int k){
    float sum = 0.0f;
    int i = 1;
    while(i <= k){
        sum += 1.0f / (static_cast<float>(i) * static_cast<float>(i));
        i++;
    }
    return sum;
}

double general_harmonic_sum_double(int k){ 
    double sum = 0.0;
    for(double i = 1.0; i <= k; i++){
        sum += 1.0/(i*i);
    }
    return sum;
}

double backwards_general_harmonic_sum_double(int k){ 
    double sum = 0.0;
    for(double i = k; i > 0; i--){
        sum += 1.0/(i*i);
    }
    return sum;
}

double backwards_general_harmonic_sum_single(int k){ 
    float sum = 0.0f;
    for(int i = k; i > 0; i--){
        float j = i;
        sum += 1.0f/(j*j);
    }
    return sum;
}

int main(){
    int k1 = 1e3;
    int k2 = 1e6;
    int k3 = 1e7;
    int k4 = 1e8;

    auto start1 = high_resolution_clock::now();
    float sum1_single = general_harmonic_sum_single(k1);
    auto stop1 = high_resolution_clock::now();

    auto start2 = high_resolution_clock::now();
    float sum2_single = general_harmonic_sum_single(k2);
    auto stop2 = high_resolution_clock::now();

    auto start3 = high_resolution_clock::now();
    float sum3_single = general_harmonic_sum_single(k3);
    auto stop3 = high_resolution_clock::now();

    auto start4 = high_resolution_clock::now();
    float sum4_single = general_harmonic_sum_single(k4);
    auto stop4 = high_resolution_clock::now();

    auto start1_backwards = high_resolution_clock::now();
    float sum1_backwards_single = backwards_general_harmonic_sum_single(k1);
    auto stop1_backwards = high_resolution_clock::now();

    auto start2_backwards = high_resolution_clock::now();
    float sum2_backwards_single = backwards_general_harmonic_sum_single(k2);
    auto stop2_backwards = high_resolution_clock::now();

    auto start3_backwards = high_resolution_clock::now();
    float sum3_backwards_single = backwards_general_harmonic_sum_single(k3);
    auto stop3_backwards = high_resolution_clock::now();

    auto start4_backwards = high_resolution_clock::now();
    float sum4_backwards_single = backwards_general_harmonic_sum_single(k4);
    auto stop4_backwards = high_resolution_clock::now();

    auto duration1_backwards = duration_cast<microseconds>(stop1_backwards - start1_backwards);
    auto duration2_backwards = duration_cast<microseconds>(stop2_backwards - start2_backwards);
    auto duration3_backwards = duration_cast<microseconds>(stop3_backwards - start3_backwards);
    auto duration4_backwards = duration_cast<microseconds>(stop4_backwards - start4_backwards);

    auto duration1 = duration_cast<microseconds>(stop1 - start1);
    auto duration2 = duration_cast<microseconds>(stop2 - start2);
    auto duration3 = duration_cast<microseconds>(stop3 - start3);
    auto duration4 = duration_cast<microseconds>(stop4 - start4);

    std::cout << "Single precision:"<<std::endl;
    std::cout << "Harmonic forward sum for k = " << k1 << ": " <<std::setprecision(15)<< sum1_single << std::endl;
    std::cout << "Time needed: " << duration1.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k1 << ": " <<std::setprecision(15)<< sum1_backwards_single << std::endl;
    std::cout << "Time needed: " << duration1_backwards.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k2 << ": " << sum2_single << std::endl;
    std::cout << "Time needed: " << duration2.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k2 << ": " << sum2_backwards_single << std::endl;
    std::cout << "Time needed: " << duration2_backwards.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k3 << ": " << sum3_single << std::endl;
    std::cout << "Time needed: " << duration3.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k3 << ": " << sum3_backwards_single << std::endl;
    std::cout << "Time needed: " << duration3_backwards.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k4 << ": " << sum4_single << std::endl;
    std::cout << "Time needed: " << duration4.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k4 << ": " << sum4_backwards_single << std::endl;
    std::cout << "Time needed: " << duration4_backwards.count() << " microseconds" << std::endl;

    // Test with double precision

    auto start1_double = high_resolution_clock::now();
    double sum1_double = general_harmonic_sum_double(k1);
    auto stop1_double = high_resolution_clock::now();

    auto start2_double = high_resolution_clock::now();
    double sum2_double = general_harmonic_sum_double(k2);
    auto stop2_double = high_resolution_clock::now();

    auto start3_double = high_resolution_clock::now();
    double sum3_double = general_harmonic_sum_double(k3);
    auto stop3_double = high_resolution_clock::now();

    auto start4_double = high_resolution_clock::now();
    double sum4_double = general_harmonic_sum_double(k4);
    auto stop4_double = high_resolution_clock::now();

    auto start1_backwards_double = high_resolution_clock::now();
    double sum1_backwards_double = backwards_general_harmonic_sum_double(k1);
    auto stop1_backwards_double = high_resolution_clock::now();

    auto start2_backwards_double = high_resolution_clock::now();
    double sum2_backwards_double = backwards_general_harmonic_sum_double(k2);
    auto stop2_backwards_double = high_resolution_clock::now();

    auto start3_backwards_double = high_resolution_clock::now();
    double sum3_backwards_double = backwards_general_harmonic_sum_double(k3);
    auto stop3_backwards_double = high_resolution_clock::now();

    auto start4_backwards_double = high_resolution_clock::now();
    double sum4_backwards_double = backwards_general_harmonic_sum_double(k4);
    auto stop4_backwards_double = high_resolution_clock::now();

    auto duration1_backwards_double = duration_cast<microseconds>(stop1_backwards_double - start1_backwards_double);
    auto duration2_backwards_double = duration_cast<microseconds>(stop2_backwards_double - start2_backwards_double);
    auto duration3_backwards_double = duration_cast<microseconds>(stop3_backwards_double - start3_backwards_double);
    auto duration4_backwards_double = duration_cast<microseconds>(stop4_backwards_double - start4_backwards_double);

    auto duration1_double = duration_cast<microseconds>(stop1_double - start1_double);
    auto duration2_double = duration_cast<microseconds>(stop2_double - start2_double);
    auto duration3_double = duration_cast<microseconds>(stop3_double - start3_double);
    auto duration4_double = duration_cast<microseconds>(stop4_double - start4_double);


    std::cout << "Double precision:"<<std::endl;
    std::cout << "Harmonic forward sum for k = " << k1 << ": " << sum1_double << std::endl;
    std::cout << "Time needed: " << duration1_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k1 << ": " << sum1_backwards_double << std::endl;
    std::cout << "Time needed: " << duration1_backwards_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k2 << ": " << sum2_double << std::endl;
    std::cout << "Time needed: " << duration2_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k2 << ": " << sum2_backwards_double << std::endl;
    std::cout << "Time needed: " << duration2_backwards_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k3 << ": " << sum3_double << std::endl;
    std::cout << "Time needed: " << duration3_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k3 << ": " << sum3_backwards_double << std::endl;
    std::cout << "Time needed: " << duration3_backwards_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic forward sum for k = " << k4 << ": " << sum4_double << std::endl;
    std::cout << "Time needed: " << duration4_double.count() << " microseconds" << std::endl;
    std::cout << "Harmonic backwards sum for k = " << k4 << ": " << sum4_backwards_double << std::endl;
    std::cout << "Time needed: " << duration4_backwards_double.count() << " microseconds" << std::endl;
    

}