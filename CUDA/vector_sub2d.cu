#include <iostream>
#include <cuda_runtime.h>

using namespace std;

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

// Kernel for 2D matrix subtraction
__global__ void matrixSub(int *A, int *B, int *C, int rows, int cols)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < rows && col < cols)
    {
        int index = row * cols + col;
        C[index] = A[index] - B[index];
    }
}

int main()
{
    int rows = 3, cols = 3;
    int size = rows * cols * sizeof(int);

    int h_A[9], h_B[9], h_C[9];

    // Initialize Matrix A
    cout << "Matrix A:\n";
    for (int i = 0; i < rows * cols; i++)
    {
        h_A[i] = i + 10;
        cout << h_A[i] << " ";
        if ((i + 1) % cols == 0) cout << endl;
    }

    // Initialize Matrix B
    cout << "\nMatrix B:\n";
    for (int i = 0; i < rows * cols; i++)
    {
        h_B[i] = i + 1;
        cout << h_B[i] << " ";
        if ((i + 1) % cols == 0) cout << endl;
    }

    int *d_A, *d_B, *d_C;

    // Allocate device memory
    CUDA_CHECK(cudaMalloc((void**)&d_A, size));
    CUDA_CHECK(cudaMalloc((void**)&d_B, size));
    CUDA_CHECK(cudaMalloc((void**)&d_C, size));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice));

    // Define block and grid size
    dim3 threadsPerBlock(16, 16);
    dim3 blocksPerGrid((cols + 15) / 16, (rows + 15) / 16);

    // Launch kernel
    matrixSub<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, rows, cols);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back
    CUDA_CHECK(cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost));

    // Print result
    cout << "\nMatrix C (A - B):\n";
    for (int i = 0; i < rows * cols; i++)
    {
        cout << h_C[i] << " ";
        if ((i + 1) % cols == 0) cout << endl;
    }

    // Free memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}