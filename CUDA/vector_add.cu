#include <iostream>
#include <cuda_runtime.h>

#define N 10

// Error checking macro
#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess)
    {
        std::cerr << "CUDA Error: " << cudaGetErrorString(code)
                  << " " << file << " " << line << std::endl;
        exit(code);
    }
}

// Kernel
__global__ void vectorAdd(int *a, int *b, int *c, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < n)
    {
        c[i] = a[i] + b[i];
    }
}

int main()
{
    int h_a[N], h_b[N], h_c[N];

    // Initialize input vectors
    for (int i = 0; i < N; i++)
    {
        h_a[i] = i;
        h_b[i] = i * 2;
    }

    int *d_a, *d_b, *d_c;

    // Allocate GPU memory
    CUDA_CHECK(cudaMalloc((void **)&d_a, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void **)&d_b, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc((void **)&d_c, N * sizeof(int)));

    // Copy data to GPU
    CUDA_CHECK(cudaMemcpy(d_a, h_a, N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, N * sizeof(int), cudaMemcpyHostToDevice));

    // Launch kernel
    vectorAdd<<<1, N>>>(d_a, d_b, d_c, N);

    // Check kernel launch error
    CUDA_CHECK(cudaGetLastError());

    // Wait for GPU to finish
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to CPU
    CUDA_CHECK(cudaMemcpy(h_c, d_c, N * sizeof(int), cudaMemcpyDeviceToHost));

    // ✅ PRINT RESULT
    std::cout << "Result vector: ";
    for (int i = 0; i < N; i++)
    {
        std::cout << h_c[i] << " ";
    }
    std::cout << std::endl;

    // Free GPU memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    cudaDeviceReset();

    return 0;
}
//nvcc -arch=sm_86 -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64" vector_add.cu -o vector_add.exe
