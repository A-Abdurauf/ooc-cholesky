/**
 *
 * @file control/descriptor.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifdef PLASMA_WITH_CUDA
#define _GNU_SOURCE
#include <cuda_runtime_api.h>
#include <sched.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#endif
#include <stdlib.h>
#include <assert.h>
#include "common.h"

/***************************************************************************//**
 *  Check for descriptor correctness
 **/
int plasma_desc_check(PLASMA_desc *desc)
{
    if (desc == NULL) {
        plasma_error("plasma_desc_check", "NULL descriptor");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (desc->mat == NULL) {
        plasma_error("plasma_desc_check", "NULL matrix pointer");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (desc->dtyp != PlasmaRealFloat &&
        desc->dtyp != PlasmaRealDouble &&
        desc->dtyp != PlasmaComplexFloat &&
        desc->dtyp != PlasmaComplexDouble  ) {
        plasma_error("plasma_desc_check", "invalid matrix type");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->mb <= 0 || desc->nb <= 0) {
        plasma_error("plasma_desc_check", "negative tile dimension");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->bsiz < desc->mb*desc->nb) {
        plasma_error("plasma_desc_check", "tile memory size smaller than the product of dimensions");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->m < 0) || (desc->n < 0)) {
        plasma_error("plasma_desc_check", "negative matrix dimension");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->lm < desc->m) || (desc->ln < desc->n)) {
        plasma_error("plasma_desc_check", "matrix dimensions larger than leading dimensions");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->i > 0 && desc->i >= desc->lm) || (desc->j > 0 && desc->j >= desc->ln)) {
        plasma_error("plasma_desc_check", "beginning of the matrix out of scope");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->i+desc->m > desc->lm || desc->j+desc->n > desc->ln) {
        plasma_error("plasma_desc_check", "submatrix out of scope");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->i % desc->mb != 0) || (desc->j % desc->nb != 0)) {
        plasma_error("plasma_desc_check", "Submatrix has to start at beginning of tiles");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
#ifdef PLASMA_WITH_CUDA
int plasma_desc_mat_alloc(PLASMA_desc *desc)
{
    size_t size;

    size = (size_t)desc->lm * (size_t)desc->ln * (size_t)plasma_element_size(desc->dtyp);
    CHECK_CUDA(cudaMallocHost(&desc->mat, size));
    if (desc->mat == NULL) {
        plasma_error("plasma_desc_mat_alloc", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    return PLASMA_SUCCESS;
}

static void set_current_thread_affinity(int core_id) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t current_thread = pthread_self();
    int result = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
    if (result != 0) {
        perror("pthread_setaffinity_np");
    }
}

static void reset_current_thread_affinity() {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
    for (int i = 0; i < num_cores; i++) {
        CPU_SET(i, &cpuset);
    }

    pthread_t current_thread = pthread_self();
    int result = pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
    if (result != 0) {
        perror("pthread_setaffinity_np");
    }
}

int plasma_desc_mat_alloc_affinity(PLASMA_desc *desc, const plasma_context_t *plasma)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    int affinity[deviceCount];
    int uniqueDevice = 0;
    for (int deviceIdx = 0; deviceIdx < deviceCount; ++deviceIdx) {
        if (plasma->rootRankOfDevice[deviceIdx] == -1) {
          continue;
        }
        affinity[uniqueDevice] = deviceIdx;
        ++uniqueDevice;
    }
    desc->uniqueDevice = uniqueDevice;
    // TODO(Jie): this needs to run exclusively, we can see all devices!
    int processor_gap = sysconf(_SC_NPROCESSORS_ONLN) / deviceCount;
    void **mat = malloc(uniqueDevice * sizeof(void *));
    desc->mat = mat;
    const int nt_each_each_processors = desc->mt / uniqueDevice;
    for (int i = 0; i < uniqueDevice - 1; ++i) {
      size_t size;
      size = (size_t)desc->mb * nt_each_each_processors * (size_t)desc->ln * (size_t)plasma_element_size(desc->dtyp);
      set_current_thread_affinity(i * processor_gap);
      CHECK_CUDA(cudaMallocHost(&mat[i], size));
    }
    size_t size;

    size = (desc->lm - (size_t)desc->mb * nt_each_each_processors * (uniqueDevice - 1)) * (size_t)desc->ln * (size_t)plasma_element_size(desc->dtyp);
    set_current_thread_affinity((uniqueDevice - 1) * processor_gap);
    CHECK_CUDA(cudaMallocHost(&mat[uniqueDevice - 1], size));
    if (desc->mat == NULL) {
        plasma_error("plasma_desc_mat_alloc", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    reset_current_thread_affinity();

    return PLASMA_SUCCESS;
}
#else
int plasma_desc_mat_alloc(PLASMA_desc *desc)
{
    size_t size;

    size = (size_t)desc->lm * (size_t)desc->ln * (size_t)plasma_element_size(desc->dtyp);
    if ((desc->mat = malloc(size)) == NULL) {
        plasma_error("plasma_desc_mat_alloc", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    return PLASMA_SUCCESS;
}
#endif

/***************************************************************************//**
 *
 **/
#ifdef PLASMA_WITH_CUDA
int plasma_desc_mat_free(PLASMA_desc *desc)
{
    if (desc->mat != NULL) {
        CHECK_CUDA(cudaFreeHost(desc->mat));
        desc->mat = NULL;
    }
    return PLASMA_SUCCESS;
}

int plasma_desc_mat_free_affinity(PLASMA_desc *desc)
{
    if (desc->mat != NULL) {
        for (int i = 0; i < desc->uniqueDevice; ++i) {
          CHECK_CUDA(cudaFreeHost(((void **)desc->mat)[i]));
        }
        free(desc->mat);
        desc->mat = NULL;
    }
    return PLASMA_SUCCESS;
}
#else
int plasma_desc_mat_free(PLASMA_desc *desc)
{
    if (desc->mat != NULL) {
        free(desc->mat);
        desc->mat = NULL;
    }
    return PLASMA_SUCCESS;
}
#endif

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Desc_Create - Create matrix descriptor.
 *
 *******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg PlasmaRealFloat:     single precision real (S),
 *          @arg PlasmaRealDouble:    double precision real (D),
 *          @arg PlasmaComplexFloat:  single precision complex (C),
 *          @arg PlasmaComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in elements (mb*nb).
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column index to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Desc_Create(PLASMA_desc **desc, void *mat, PLASMA_enum dtyp, int mb, int nb, int bsiz,
                       int lm, int ln, int i, int j, int m, int n)
{
    plasma_context_t *plasma;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_Desc_Create", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Allocate memory and initialize the descriptor */
    *desc = (PLASMA_desc*)malloc(sizeof(PLASMA_desc));
    if (*desc == NULL) {
        plasma_error("PLASMA_Desc_Create", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    **desc = plasma_desc_init(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n);
    (**desc).mat = mat;
    status = plasma_desc_check(*desc);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_Desc_Create", "invalid descriptor");
        return status;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Desc_Destroy - Destroys matrix descriptor.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Desc_Destroy(PLASMA_desc **desc)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_Desc_Destroy", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        plasma_error("PLASMA_Desc_Destroy", "attempting to destroy a NULL descriptor");
        return PLASMA_ERR_UNALLOCATED;
    }
    free(*desc);
    *desc = NULL;
    return PLASMA_SUCCESS;
}
