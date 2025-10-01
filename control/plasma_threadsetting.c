/**
 *
 * @file plasma_threadsetting.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 *
 **/
#include "common.h"
/***************************************************************************//**
 * switch lapack thread_num initialization
 **/
#if defined(PLASMA_WITH_ACML) || defined(PLASMA_WITH_MKL)
#include <omp.h>
#endif

#if defined(PLASMA_WITH_MKL)
#include <mkl_service.h>
#endif

#if defined(PLASMA_WITH_OPENBLAS)
    void openblas_set_num_threads(int threads);
#endif
/////////////////////////////////////////////////////////////
static inline void plasma_setlapack_numthreads(int num_threads)
{
#if defined(PLASMA_WITH_ACML)
    omp_set_num_threads( num_threads );
#elif defined(PLASMA_WITH_MKL)
    mkl_set_num_threads( num_threads );
#elif defined(PLASMA_WITH_OPENBLAS)
    openblas_set_num_threads( num_threads );
#else
    if (num_threads > 1) {
        plasma_warning( "plasma_setlapack_numthreads",
                        "\n"
                        "==========================================================================================\n"
                        "  WARNING you are calling a parallel section without linking with a multithread library   \n"
                        "  please compile with multi thread library and add -fopenmp(gcc)/-openmp(icc) to both     \n"
                        "  compilation and linkage flags                                                           \n"
                        "==========================================================================================\n");
    }
#endif
}
/***************************************************************************//**
 * switch lapack to multi threading and unset plasma affinity
 **/
void plasma_setlapack_multithreads(int num_threads)
{
    plasma_unsetaffinity();
    plasma_setlapack_numthreads(num_threads);
}

/***************************************************************************//**
 * switch lapack to sequential and setting plasma affinity
 **/
void plasma_setlapack_sequential(plasma_context_t *plasma)
{
    /*
     * All threads need to call the parallel section that
     * set omp/mkl numthreads to 1, then only the master (0) will
     * set plasma affinity
     */
    plasma_static_call(plasma_psetlapack_numthreads);
    plasma_setaffinity(plasma->thread_bind[0]);
}

void plasma_psetlapack_numthreads(plasma_context_t *plasma)
{
    plasma_setlapack_numthreads(1);
}
