#include <iostream>
#include <cuda_runtime.h>

#define N 3   // Matrix size (N x N)

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

// Kernel for matrix addition
__global__ void matrixAdd(int *A, int *B, int *C, int width)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < width && col < width)
    {
        int index = row * width + col;
        C[index] = A[index] + B[index];
    }
}

int main()
{
    int size = N * N * sizeof(int);

    int h_A[N][N], h_B[N][N], h_C[N][N];

    // Initialize matrices
    std::cout << "Matrix A:\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            h_A[i][j] = i + j;
            std::cout << h_A[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nMatrix B:\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            h_B[i][j] = i * j;
            std::cout << h_B[i][j] << " ";
        }
        std::cout << std::endl;
    }

    int *d_A, *d_B, *d_C;

    // Allocate GPU memory
    CUDA_CHECK(cudaMalloc((void **)&d_A, size));
    CUDA_CHECK(cudaMalloc((void **)&d_B, size));
    CUDA_CHECK(cudaMalloc((void **)&d_C, size));

    // Copy data to GPU
    CUDA_CHECK(cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice));

    // Define block and grid size (2D)
    dim3 block(16, 16);
    dim3 grid((N + block.x - 1) / block.x, (N + block.y - 1) / block.y);

    // Launch kernel
    matrixAdd<<<grid, block>>>(d_A, d_B, d_C, N);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to CPU
    CUDA_CHECK(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

    // Print result
    std::cout << "\nResult Matrix C (A + B):\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::cout << h_C[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Free memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cudaDeviceReset();

    return 0;
}