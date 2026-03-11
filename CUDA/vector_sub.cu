#include <iostream>
#include <cuda_runtime.h>

using namespace std;

#define N 10

// Error checking macro
#define CUDA_CHECK(call)                                      \
{                                                             \
    cudaError_t err = call;                                   \
    if (err != cudaSuccess)                                   \
    {                                                         \
        cerr << "CUDA Error: " << cudaGetErrorString(err)     \
             << " at line " << __LINE__ << endl;              \
        exit(EXIT_FAILURE);                                   \
    }                                                         \
}

// Kernel for vector subtraction
__global__ void vectorSub(int *a, int *b, int *c)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < N)
    {
        c[i] = a[i] - b[i];
    }
}

int main()
{
    int h_a[N], h_b[N], h_c[N];

    // Initialize host arrays
    for (int i = 0; i < N; i++)
    {
        h_a[i] = i + 1;     // 1 2 3 4 5 6 7 8 9 10
        h_b[i] = i;         // 0 1 2 3 4 5 6 7 8 9
    }

    int *d_a, *d_b, *d_c;

    // Allocate device memory
    CUDA_CHECK(cudaMalloc((void**)&d_a, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_c, N * sizeof(int)));

    // Copy data from host to device
    CUDA_CHECK(cudaMemcpy(d_a, h_a, N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, N * sizeof(int), cudaMemcpyHostToDevice));

    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    vectorSub<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c);

    // Check for kernel launch error
    CUDA_CHECK(cudaGetLastError());

    // Synchronize
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(h_c, d_c, N * sizeof(int), cudaMemcpyDeviceToHost));

    // Print result
    cout << "Result vector: ";
    for (int i = 0; i < N; i++)
    {
        cout << h_c[i] << " ";
    }
    cout << endl;

    // Free device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    return 0;
}