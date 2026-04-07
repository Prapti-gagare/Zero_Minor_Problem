#include <iostream>
#include <cuda_runtime.h>
#include <math.h>

// GPU kernel
__global__ void swiglu_kernel(const float* input, float* output, int halfN) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < halfN) {
        float x = input[i];              // first half
        float y = input[i + halfN];      // second half

        // sigmoid(x)
        float sigmoid_x = 1.0f / (1.0f + expf(-x));

        // SiLU(x) = x * sigmoid(x)
        float silu = x * sigmoid_x;

        // SWiGLU = SiLU(x) * y
        output[i] = silu * y;
    }
}

// Host function (given)
extern "C" void solve(const float* input, float* output, int N) {
    int halfN = N / 2;

    int threadsPerBlock = 256;
    int blocksPerGrid = (halfN + threadsPerBlock - 1) / threadsPerBlock;

    swiglu_kernel<<<blocksPerGrid, threadsPerBlock>>>(input, output, halfN);

    cudaDeviceSynchronize();
}

int main() {
    int N = 8; // must be even
    int halfN = N / 2;

    // Input: first half = x, second half = y
    float h_input[] = {1, 2, 3, 4,   5, 6, 7, 8};

    // Output size = N/2
    float* h_output = new float[halfN];

    float *d_input, *d_output;

    // Allocate GPU memory
    cudaMalloc((void**)&d_input, N * sizeof(float));
    cudaMalloc((void**)&d_output, halfN * sizeof(float));

    // Copy input to GPU
    cudaMemcpy(d_input, h_input, N * sizeof(float), cudaMemcpyHostToDevice);

    // Call solve (which launches kernel)
    solve(d_input, d_output, N);

    // Copy result back
    cudaMemcpy(h_output, d_output, halfN * sizeof(float), cudaMemcpyDeviceToHost);

    // Print output
    std::cout << "Output: ";
    for (int i = 0; i < halfN; i++) {
        std::cout << h_output[i] << " ";
    }
    std::cout << std::endl;

    // Free memory
    cudaFree(d_input);
    cudaFree(d_output);
    delete[] h_output;

    return 0;
}