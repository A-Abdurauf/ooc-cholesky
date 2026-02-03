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

static bool parse_int_arg(const std::string &name, const std::string &value,
                          int &out) {
  try {
    size_t pos = 0;
    int v = std::stoi(value, &pos);
    if (pos != value.size()) {
      std::cerr << "Invalid value for " << name << ": '" << value << "'\n";
      return false;
    }
    out = v;
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Invalid value for " << name << ": '" << value << "'";
    std::cerr << " (" << e.what() << ")\n";
    return false;
  }
}

int main(int argc, char **argv) {
  int cores = 2;
  int N = 1024; // default; can be inferred from bin if omitted/mismatched
  int nb = 128;
  std::string bin_path = "/home/abduraa/MX_project/logs/my_cov_weak_1024.bin";
  std::string format;
  std::string mx_mode;
  std::string fp32_bucket;
  std::string fp16_bucket;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--n" && i + 1 < argc) {
      if (!parse_int_arg("--n", argv[++i], N)) return 2;
    } else if (arg == "--nb" && i + 1 < argc) {
      if (!parse_int_arg("--nb", argv[++i], nb)) return 2;
    } else if (arg == "--bin" && i + 1 < argc) {
      bin_path = argv[++i];
    } else if (arg == "--cores" && i + 1 < argc) {
      if (!parse_int_arg("--cores", argv[++i], cores)) return 2;
    } else if (arg == "--format" && i + 1 < argc) {
      format = argv[++i];
    } else if (arg == "--mx-mode" && i + 1 < argc) {
      mx_mode = argv[++i];
    } else if (arg == "--fp32-bucket" && i + 1 < argc) {
      fp32_bucket = argv[++i];
    } else if (arg == "--fp16-bucket" && i + 1 < argc) {
      fp16_bucket = argv[++i];
    }
  }

  if (!format.empty()) {
    setenv("MX_FORCE_FORMAT", format.c_str(), 1);
  }
  if (!mx_mode.empty()) {
    setenv("MX_MX_MODE", mx_mode.c_str(), 1);
  }
  if (!fp32_bucket.empty()) {
    setenv("MX_BUCKET_FP32", fp32_bucket.c_str(), 1);
  }
  if (!fp16_bucket.empty()) {
    setenv("MX_BUCKET_FP16", fp16_bucket.c_str(), 1);
  }

  std::ifstream in(bin_path, std::ios::binary | std::ios::ate);
  if (!in) {
    throw std::runtime_error("cannot open input bin");
  }
  const std::streamsize bytes = in.tellg();
  in.seekg(0, std::ios::beg);

  const size_t count = static_cast<size_t>(bytes) / sizeof(double);
  const size_t dim = static_cast<size_t>(std::llround(std::sqrt(count)));
  if (dim * dim != count) {
    throw std::runtime_error("bin file size is not a square matrix");
  }
  if (N <= 0 || static_cast<size_t>(N) * N != count) {
    N = static_cast<int>(dim);
  }

  std::vector<double> data(static_cast<size_t>(N) * N);
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
