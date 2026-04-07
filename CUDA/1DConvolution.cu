#include <iostream>
#include <cuda_runtime.h>

// Kernel function
__global__ void convolution_1d_kernel(const float* input,
                                      const float* kernel,
                                      float* output,
                                      int input_size,
                                      int kernel_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int output_size = input_size - kernel_size + 1;

    if (idx < output_size) {
        float sum = 0.0f;

        for (int j = 0; j < kernel_size; j++) {
            sum += input[idx + j] * kernel[j];
        }

        output[idx] = sum;
    }
}

int main() {
    int input_size = 8;
    int kernel_size = 3;

    float h_input[]  = {1, 2, 3, 4, 5, 6, 7, 8};
    float h_kernel[] = {1, 0, -1};

    int output_size = input_size - kernel_size + 1;

    // ✅ Dynamic allocation (FIXED)
    float* h_output = new float[output_size];

    float *d_input, *d_kernel, *d_output;

    // Allocate memory on GPU
    cudaMalloc((void**)&d_input, input_size * sizeof(float));
    cudaMalloc((void**)&d_kernel, kernel_size * sizeof(float));
    cudaMalloc((void**)&d_output, output_size * sizeof(float));

    // Copy data to GPU
    cudaMemcpy(d_input, h_input, input_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_kernel, h_kernel, kernel_size * sizeof(float), cudaMemcpyHostToDevice);

    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (output_size + threadsPerBlock - 1) / threadsPerBlock;

    convolution_1d_kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_input, d_kernel, d_output, input_size, kernel_size
    );

    // Wait for GPU to finish
    cudaDeviceSynchronize();

    // Copy result back to CPU
    cudaMemcpy(h_output, d_output, output_size * sizeof(float), cudaMemcpyDeviceToHost);

    // Print output
    std::cout << "Output: ";
    for (int i = 0; i < output_size; i++) {
        std::cout << h_output[i] << " ";
    }
    std::cout << std::endl;

    // Free GPU memory
    cudaFree(d_input);
    cudaFree(d_kernel);
    cudaFree(d_output);

    // Free CPU memory
    delete[] h_output;

    return 0;
}