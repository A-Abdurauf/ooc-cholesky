/**
 *
 * @file example_dpotrf.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example of Cholesky factorization
 *
 * @version 2.8.0
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:03:07 2016
 *
 **/
#include <Eigen/Eigen>

#include <cblas.h>
#include <common.h>
#include <context.h>
#include <math.h>
#include <plasma.h>
#include <plasma_d_mixed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#undef max
#undef min
#include <chrono>
#include <iostream>
#include <fstream>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <core_blas.h>

int check_factorization(int, double *, double *, int, int);

int IONE = 1;
int ISEED[4] = {0, 0, 0, 1}; /* initial seed for dlarnv() */

int main(int argc, char **argv) {
  int cores = 2;
  int N = 1024; // match your covariance size
  int nb = 128;
  std::string bin_path = "/home/abduraa/MX_project/logs/my_cov_weak_1024.bin";
  std::string format;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--n" && i + 1 < argc) {
      N = std::stoi(argv[++i]);
    } else if (arg == "--nb" && i + 1 < argc) {
      nb = std::stoi(argv[++i]);
    } else if (arg == "--bin" && i + 1 < argc) {
      bin_path = argv[++i];
    } else if (arg == "--cores" && i + 1 < argc) {
      cores = std::stoi(argv[++i]);
    } else if (arg == "--format" && i + 1 < argc) {
      format = argv[++i];
    }
  }

  if (!format.empty()) {
    setenv("MX_FORCE_FORMAT", format.c_str(), 1);
  }

  std::vector<double> data(static_cast<size_t>(N) * N);
  std::ifstream in(bin_path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot open input bin");
  }
  in.read(reinterpret_cast<char *>(data.data()),
          static_cast<std::streamsize>(data.size() * sizeof(double)));
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>>
      A(data.data(), N, N);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L = A;

  PLASMA_Init(cores);

  plasma_context_t *plasma = plasma_context_self();
  plasma->autotuning_enabled = 0;
  plasma->nb = nb;  // tile size

  auto start = std::chrono::high_resolution_clock::now();

  // Call the tile-first path to avoid the column-major → tile transform in
  // PLASMA_dpotrf_gpu.
  const auto status =
      PLASMA_dpotrf_gpu_reuse_data_table_mixed_precision(PlasmaLower, N,
                                                         L.data(), N);
  if (status != 0) {
    printf("factorization failed in %s:%d\n", __FILE__, __LINE__);
  }
  auto duration = std::chrono::high_resolution_clock::now() - start;
  // print duration
  std::cout << "Runtime: " << duration.count() << " microseconds" << std::endl;

  L.triangularView<Eigen::StrictlyUpper>().setZero();
  double error = 0;
  A = A.selfadjointView<Eigen::Lower>();
  error = (A - L * L.transpose()).array().abs().rowwise().sum().maxCoeff();

  std::cout << "error: " << error << std::endl;

  PLASMA_Finalize();

  return EXIT_SUCCESS;
}
