#include <iostream>
#include <cuda_runtime.h>

__global__ void add(double *a, double *b, double *c) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    c[index] = a[index] + b[index];
}

int main(){
    const int N = 10;
    double h_a[N], h_b[N], h_c[N];

    // Init Host Arrays
    for (int i = 0; i < N; ++i) {
        h_a[i] = i;
        h_b[i] = 100 + i;
    }

    // Device Pointer
    double *d_a, *d_b, *d_c;
    cudaMalloc(&d_a, N * sizeof(double));
    cudaMalloc(&d_b, N * sizeof(double));
    cudaMalloc(&d_c, N * sizeof(double));

    // Kopieren: Host → Device
    cudaMemcpy(d_a, h_a, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, N * sizeof(double), cudaMemcpyHostToDevice);

    // CUDA-Kernel starten
    int threadsPerBlock = 256;
    int blocks = (N + threadsPerBlock - 1) / threadsPerBlock;
    add<<<blocks, threadsPerBlock>>>(d_a, d_b, d_c);

    // Ergebnis zurück: Device → Host
    cudaMemcpy(h_c, d_c, N * sizeof(double), cudaMemcpyDeviceToHost);

    // Ausgabe
    for (int i = 0; i < N; ++i)
        std::cout << h_a[i] << " + " << h_b[i] << " = " << h_c[i] << "\n";

    // Aufräumen
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    return 0;
    
}