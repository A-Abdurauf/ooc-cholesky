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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <core_blas.h>

int check_factorization(int, double *, double *, int, int);

int IONE = 1;
int ISEED[4] = {0, 0, 0, 1}; /* initial seed for dlarnv() */

int main() {
  int cores = 2;
  int N = 256;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A{N, N};
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L{N, N};

  PLASMA_Init(cores);

  plasma_context_t *plasma = plasma_context_self();
  plasma->autotuning_enabled = 0;
  plasma->nb = 128;

  A.setRandom();
  A = (A * A.transpose()).eval();
  A.diagonal().array() += N;
  L = A;

  const auto status = PLASMA_dpotrf_gpu(PlasmaLower, N, L.data(), N);
  if (status != 0) {
    printf("factorization failed in %s:%d\n", __FILE__, __LINE__);
  }

  L.triangularView<Eigen::StrictlyUpper>().setZero();
  double error = 0;
  A = A.selfadjointView<Eigen::Lower>();
  error = (A - L * L.transpose()).array().abs().rowwise().sum().maxCoeff();

  std::cout << "error: " << error << std::endl;

  PLASMA_Finalize();

  return EXIT_SUCCESS;
}
