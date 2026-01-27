#include <iostream>
#include <Eigen/Dense>
#include "plasma.h"
#include "plasma_d_mixed.h"
#include <chrono>
#include <lapacke.h>
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

    // --- Create clean SPD matrix ---
    A.setRandom();
    A = (A * A.transpose()).eval();
    A.diagonal().array() += N;
    L = A;

    // --- Allocate contiguous buffer for tiled storage ---
    double *A_data = (double*)malloc(sizeof(double) * N * N);
    if(!A_data) { std::cerr << "Failed to allocate buffer" << std::endl; return EXIT_FAILURE; }

    // --- Allocate tiled descriptor (already-tiled matrix storage) ---
    plasma_desc_t *descA = nullptr;
    int info = PLASMA_Desc_Create(
        &descA,
        A_data,            // contiguous backing buffer
        PlasmaRealDouble,  // tile datatype (keep double for now)
        plasma->nb,        // tile rows
        plasma->nb,        // tile cols
        plasma->nb*plasma->nb,
        N, N,              // full matrix size
        0, 0, N, N         // tile offsets
    );(&descA, A_data, PlasmaRealDouble,
                                  plasma->nb, plasma->nb, plasma->nb*plasma->nb,
                                  N, N, 0, 0, N, N);
    if(info != 0) { std::cerr << "PLASMA_Desc_Create failed: " << info << std::endl; return EXIT_FAILURE; }

    // --- Pack Eigen data into tiles ---
    for(int m=0; m<descA->mt; m++) {
        for(int n=0; n<descA->nt; n++) {
            double* tile = (double*)plasma_getaddr(*descA, m, n);
            for(int i=0; i<descA->mb; i++) {
                for(int j=0; j<descA->nb; j++) {
                    int global_row = m*descA->mb + i;
                    int global_col = n*descA->nb + j;
                    if(global_row<N && global_col<N)
                        tile[i*descA->nb + j] = A(global_row, global_col);
                }
            }
        }
    }

    // --- Launch tiled mixed‑precision GPU dpotrf ---
    PLASMA_sequence *sequence = nullptr;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_Sequence_Create(&sequence);

    // Start timing only for GPU factorization
    auto start_gpu = std::chrono::high_resolution_clock::now();

    PLASMA_dpotrf_gpu_reuse_data_table_no_free_Tile_Async(
        PlasmaLower, descA, sequence, &request);

    PLASMA_Sequence_Wait(sequence);
    cudaDeviceSynchronize();

    auto end_gpu = std::chrono::high_resolution_clock::now();
    double runtime_gpu = std::chrono::duration<double>(end_gpu - start_gpu).count();
    std::cout << "GPU factorization runtime: " << runtime_gpu << " s" << std::endl;

    // --- Compute approximate GFLOPS for Cholesky (1/3*N^3)  --- (1/3*N^3) ---
    double flops = (1.0/3.0) * N * N * N;
    double gflops = flops / runtime_gpu / 1e9;
    std::cout << "Approx. GFLOPS: " << gflops << std::endl;

    // --- Unpack tiled descriptor back to Eigen ---
    for(int m=0; m<descA->mt; m++) {
        for(int n=0; n<descA->nt; n++) {
            double* tile = (double*)plasma_getaddr(*descA, m, n);
            for(int i=0; i<descA->mb; i++) {
                for(int j=0; j<descA->nb; j++) {
                    int global_row = m*descA->mb + i;
                    int global_col = n*descA->nb + j;
                    if(global_row<N && global_col<N)
                        L(global_row, global_col) = tile[i*descA->nb + j];
                }
            }
        }
    }

    L.triangularView<Eigen::StrictlyUpper>().setZero();

    // --- Compute errors ---
    A = A.selfadjointView<Eigen::Lower>();
    Eigen::MatrixXd LLT = L * L.transpose();
    double rel_err = (A - LLT).norm() / A.norm();
    double abs_err = (A - LLT).array().abs().maxCoeff();

    std::cout << "Relative factorization error: " << rel_err << std::endl;
    std::cout << "Max absolute error: " << abs_err << std::endl;

    PLASMA_Sequence_Destroy(sequence);
    PLASMA_Desc_Destroy(&descA);
    free(A_data);
    PLASMA_Finalize();
    return EXIT_SUCCESS;
}
