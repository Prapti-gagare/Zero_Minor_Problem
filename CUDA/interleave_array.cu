#include <cuda_runtime.h>
#include <iostream>
#include <vector>

__global__ void interleave_kernel(const float* A, const float* B, float* output, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        output[2 * idx]     = A[idx]; // even indices from A
        output[2 * idx + 1] = B[idx]; // odd indices from B
    }
}

extern "C" void solve(const float* A, const float* B, float* output, int N) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    interleave_kernel<<<blocksPerGrid, threadsPerBlock>>>(A, B, output, N);
    cudaDeviceSynchronize();
}

int main() {
    int N;
    std::cout << "Enter size of arrays: ";
    std::cin >> N;

    std::vector<float> h_A(N), h_B(N), h_output(2 * N);

    std::cout << "Enter elements of array A:\n";
    for (int i = 0; i < N; i++) std::cin >> h_A[i];

    std::cout << "Enter elements of array B:\n";
    for (int i = 0; i < N; i++) std::cin >> h_B[i];

    float *d_A, *d_B, *d_output;
    cudaMalloc(&d_A, N * sizeof(float));
    cudaMalloc(&d_B, N * sizeof(float));
    cudaMalloc(&d_output, 2 * N * sizeof(float));

    cudaMemcpy(d_A, h_A.data(), N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B.data(), N * sizeof(float), cudaMemcpyHostToDevice);

    solve(d_A, d_B, d_output, N);

    cudaMemcpy(h_output.data(), d_output, 2 * N * sizeof(float), cudaMemcpyDeviceToHost);

    std::cout << "Interleaved output:\n";
    for (int i = 0; i < 2 * N; i++)
        std::cout << h_output[i] << " ";
    std::cout << std::endl;

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_output);

    return 0;
}