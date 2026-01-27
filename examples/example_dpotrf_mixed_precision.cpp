#include <iostream>
#include <Eigen/Dense>

// // Must come BEFORE including plasma headers
// #define PLASMA_WITH_HPCASIA24
// #define PLASMA_WITH_CUDA

#include "plasma.h"
#include "plasma_d_mixed.h"
#include <Eigen/Eigen>
#include <cblas.h>
#include <common.h>
#include <context.h>
#include <math.h>
#include <plasma.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#undef max
#undef min
#include <chrono>

#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <core_blas.h>
int main() {
    int cores = 2;
    int N = 4096;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(N, N);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L(N, N);

    // Initialize PLASMA
    PLASMA_Init(cores);

    plasma_context_t *plasma = plasma_context_self();
    plasma->autotuning_enabled = 0;
    plasma->nb = 32;

    // Create a symmetric positive definite matrix
    A.setRandom();
    A = (A * A.transpose()).eval();
    A.diagonal().array() += N;

    // Copy A to L for factorization
    L = A;

    auto start = std::chrono::high_resolution_clock::now();

    // Call the mixed-precision routine (provided by the HPCASIA24 branch)
    PLASMA_dpotrf_gpu_reuse_data_mixed_precision(PlasmaLower, N, L.data(), N);

    cudaDeviceSynchronize();  // <-- important for correct timing
    auto duration = std::chrono::high_resolution_clock::now() - start;
    //print duration 
    std::cout << "Runtime: " << duration.count() << " microseconds" << std::endl;
    // Force lower-triangular view (zero upper part)
    L.triangularView<Eigen::StrictlyUpper>().setZero();

    // Compute reconstruction error: ‖A - LLᵀ‖
    A = A.selfadjointView<Eigen::Lower>();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> LLT = L * L.transpose();
    double error = (A - LLT).norm() / A.norm();  // relative error

    std::cout << "Relative factorization error: " << error << std::endl;

    // Optional: absolute max error
    // double abs_err = (A - LLT).array().abs().maxCoeff();
    double abs_err = (A - L * L.transpose()).array().abs().rowwise().sum().maxCoeff();
    std::cout << "Max absolute error: " << abs_err << std::endl;

    PLASMA_Finalize();
    return EXIT_SUCCESS;
}

// #include "plasma.h"
// #include "plasma_async.h"
// #include "plasma_desc.h"
// #include "plasma_types.h"
// int main() {

//     int N  = 4096;
//     int nb = 256;
//     plasma_desc_t descL;
//     plasma_sequence_t *sequence;
//     plasma_request_t request;

//     // Create descriptor
//     double *L = (double*)malloc(N * N * sizeof(double));
//     plasma_desc_create(&descL, L, PlasmaRealDouble,
//                        nb, nb, nb*nb,
//                        N, N, 0, 0, N, N);

//     // Create sequence & request
//     plasma_context_t *plasma = plasma_context_self();
//     plasma_sequence_create(plasma, &sequence);
//     plasma_request_init(&request);

//     // Call tiled dpotrf
//     PLASMA_dpotrf_gpu_reuse_data_table_no_free_Tile_Async(
//         PlasmaLower, &descL, sequence, &request);

//     plasma_sequence_wait(plasma, sequence);

//     plasma_desc_destroy(&descL);
//     free(L);
// }
