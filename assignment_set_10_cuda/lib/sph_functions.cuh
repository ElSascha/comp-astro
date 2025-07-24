#ifndef SPH_FUNCTIONS_CUH
#define SPH_FUNCTIONS_CUH

#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>

__host__ __device__ double3 add(double3 a, double3 b); // host device so it can be called from both CPU and GPU 
__host__ __device__ double3 sub(double3 a, double3 b);
__host__ __device__ double3 scal_mul(double3 a, double b);
__host__ __device__ double3 scal_div(double3 a, double b);
__host__ __device__ double3 cross(double3 a, double3 b);
__host__ __device__ double dot(double3 a, double3 b);
__host__ __device__ double length(double3 a);
__host__ __device__ double3 normalize(double3 a);
__host__ __device__ double dot(double3 a, double3 b);
__host__ __device__ double length(double3 a);

#endif // SPH_FUNCTIONS_CUH