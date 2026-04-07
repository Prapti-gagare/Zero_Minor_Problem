#include <iostream>
#include <cuda_runtime.h>

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

// CUDA kernel for matrix multiplication
__global__ void matrixMul(int *A, int *B, int *C, int N)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < N && col < N)
    {
        int sum = 0;
        for (int k = 0; k < N; k++)
            sum += A[row * N + k] * B[k * N + col];
        C[row * N + col] = sum;
    }
}

int main()
{
    int N;
    std::cout << "Enter matrix size N: ";
    scanf("%d", &N);

    int size = N * N * sizeof(int);
    int *h_A = new int[N * N];
    int *h_B = new int[N * N];
    int *h_C = new int[N * N];

    // Input matrix A
    std::cout << "Enter Matrix A (" << N << "x" << N << "):\n";
    for (int i = 0; i < N * N; i++)
        scanf("%d", &h_A[i]);

    // Input matrix B
    std::cout << "Enter Matrix B (" << N << "x" << N << "):\n";
    for (int i = 0; i < N * N; i++)
        scanf("%d", &h_B[i]);

    int *d_A, *d_B, *d_C;
    CUDA_CHECK(cudaMalloc((void **)&d_A, size));
    CUDA_CHECK(cudaMalloc((void **)&d_B, size));
    CUDA_CHECK(cudaMalloc((void **)&d_C, size));

    CUDA_CHECK(cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice));

    dim3 block(16, 16);
    dim3 grid((N + block.x - 1) / block.x, (N + block.y - 1) / block.y);

    matrixMul<<<grid, block>>>(d_A, d_B, d_C, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

    std::cout << "\nResult Matrix C (A x B):\n";
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            std::cout << h_C[i * N + j] << " ";
        std::cout << std::endl;
    }

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
    delete[] h_A;
    delete[] h_B;
    delete[] h_C;

    cudaDeviceReset();
    return 0;
}