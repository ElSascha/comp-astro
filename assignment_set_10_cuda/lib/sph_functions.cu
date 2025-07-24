#include "sph_kernel.cuh"

__host__ __device__ double3 add(double3 a, double3 b){

    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ double3 sub(double3 a, double3 b){
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ double3 scal_mul(double3 a, double b){
    return make_double3(a.x * b, a.y * b, a.z * b);
}

__host__ __device__ double3 scal_div(double3 a, double b)
{
    return make_double3(a.x / b, a.y / b, a.z / b);
}

__host__ __device__ double3 cross(double3 a, double3 b)
{
    return make_double3(a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x);
}

__host__ __device__ double dot(double3 a, double3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__host__ __device__ double length(double3 a)
{
    return sqrt(dot(a, a));
}

__host__ __device__ double3 normalize(double3 a)
{
    double len = length(a);
    if (len > 0)
    {
        return scal_div(a, len);
    }
    return a;
}
