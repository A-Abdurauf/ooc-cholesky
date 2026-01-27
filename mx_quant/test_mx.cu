#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm> // for std::max
#include <cuda_runtime.h>
#include "mx_wrapper.h"


// --- Helper to print IEEE 754 bits ---
void print_float_bits(float val) {
    uint32_t bits = *reinterpret_cast<uint32_t*>(&val);
    
    // Sign bit (1 bit)
    std::cout << ((bits >> 31) & 1) << " | ";
    
    // Exponent (8 bits)
    for (int k = 7; k >= 0; k--) {
        std::cout << ((bits >> (23 + k)) & 1);
    }
    std::cout << " | ";
    
    // Mantissa (23 bits)
    for (int k = 22; k >= 0; k--) {
        std::cout << ((bits >> k) & 1);
    }
}
// -------------------------------------

int main() {
    const int N = 32;
    std::vector<float> h_input(N);
    std::vector<float> h_output(N);
    // We still allocate N for safety, but for a single block (1x32), 
    // the kernel only reads the first element.
    std::vector<float> h_max(N); 

    // 1. Fill input
    float global_max = 0.0f;
    for (int i = 0; i < N; i++) {
        h_input[i] = (i) * 0.01f;
        // Track the absolute max of the *entire* group
        float abs_val = std::fabs(h_input[i]);
        if (abs_val > global_max) global_max = abs_val;
    }

    // 2. --- CRITICAL FIX ---
    // The kernel (with axis_size=32) treats this as ONE block.
    // It expects the max value for the WHOLE block to be at h_max[0].
    // (It does NOT read h_max[1]...h_max[31]).
    for(int i=0; i<N; ++i) {
        h_max[i] = global_max; // Fill all just to be safe, but index 0 is what matters
    }
    
    std::cout << "Global Max for block: " << global_max << std::endl;

    float *d_in, *d_out, *d_max;
    cudaMalloc(&d_in,  N * sizeof(float));
    cudaMalloc(&d_out, N * sizeof(float));
    cudaMalloc(&d_max, N * sizeof(float));

    cudaMemcpy(d_in,  h_input.data(), N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_max, h_max.data(),   N * sizeof(float), cudaMemcpyHostToDevice);

    // 3. Launch MX quantization
    // axis_size = N means "Normalize over the entire vector"
    mx_quantize_fp32_to_mx(
        d_in, d_out, d_max,
        N,          // total elements
        N,          // axis_size (Block size)
        1,          // post_axis
        8,          // scale bits
        4,          // elem_ebits
        3,          // elem_mbits (E4M3)
        448.0f,     // max norm
        true,       // flush subnorms
        0           // rounding mode
    );

    cudaMemcpy(h_output.data(), d_out, N*sizeof(float), cudaMemcpyDeviceToHost);

    // 4. Print Output
    std::cout << "--- Comparison (Input vs Quantized) ---" << std::endl;
    std::cout << std::left << std::setw(10) << "Index" 
              << std::setw(15) << "Input" 
              << std::setw(15) << "Output" << std::endl;

    for (int i = 0; i < N; i++) {
        std::cout << std::left << std::setw(10) << i 
                  << std::setw(15) << h_input[i] 
                  << std::setw(15) << h_output[i] << "\n";
    }

    // 5. Binary Print (Bitwise)
    std::cout << "--- Binary Comparison (Sign | Exponent | Mantissa) ---" << std::endl;
    for (int i = 0; i < 5; i++) {
        std::cout << "Index " << i << ":\n";
        
        std::cout << "  IN : ";
        print_float_bits(h_input[i]);
        std::cout << " (" << h_input[i] << ")\n";

        std::cout << "  OUT: ";
        print_float_bits(h_output[i]);
        std::cout << " (" << h_output[i] << ")\n";
        
        std::cout << "----------------------------------------------------" << std::endl;
    }

    cudaFree(d_in);
    cudaFree(d_out);
    cudaFree(d_max);
    return 0;
}