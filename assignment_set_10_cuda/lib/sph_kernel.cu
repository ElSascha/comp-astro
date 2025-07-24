#pragma once // Ensure this header is included only once
#include "sph_kernel.cuh"

__device__ double cubic_bspline(double r, double h){
    double cons = NORMALAIZATION_CONSTANT/pow(h, 3);
    if(0 <= r/h && r/h < 1/2){
        return cons * (6*(pow(r/h,3)) - 6*(pow(r/h,2)) + 1);
    }

    else if(1/2 <= r/h && r/h <= 1){
        return cons * 2*(pow(1-r/h,3));
    }
    
    else if(r/h > 1){
        return 0.0;
    }

    else{
        //throw exception("Invalid value for r/h in cubic_bspline");
        try {
            throw std::runtime_error("Invalid value for r/h in cubic_bspline");
        } catch (const std::runtime_error& e) {
            std::cout<<"Error: "<< e.what() << std::endl;
            return 0.0;
        } 
    }
    
}

__global__ void compute_density(Particle* particles){
    
}

__global__ void compute_pressure(Particle* particles){

}