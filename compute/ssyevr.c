/**
 * @file ssyevr.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:56 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_ssyevr - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex Hermitian matrix A. The matrix A is
 *  preliminary reduced to tridiagonal form using a two-stage
 *  approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] range
 *          = PlasmaAllVec: all eigenvalues will be found.
 *          = PlasmaVec: all eigenvalues in the half-open interval (VL,VU]
 *                       will be found.
 *          = PlasmaIvec: the IL-th through IU-th eigenvalues will be found.
 *          For range = PlasmaVec or PlasmaIvec and IU - IL < N - 1,
 *          DSTEBZ and ZSTEIN are called.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = PlasmaLower) or the
 *          upper triangle (if uplo = PlasmaUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] vl
 *          see vu.
 *
 * @param[in] vu
 *          If RANGE=PlasmaVec, the lower and upper bounds of the
 *          interval to be searched for eigenvalues. VL < VU.
 *          Not referenced if RANGE = PlasmaAllVec or PlasmaIvec.
 *
 * @param[in] il
 *          see iu.
 *
 * @param[in] iu
 *          If RANGE=PlasmaIvec, the indices (in ascending order) of the
 *          smallest and largest eigenvalues to be returned.
 *          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
 *          Not referenced if RANGE = PlasmaAllVec or PlasmaVec.
 *
 * @param[in] abstol
 *          The absolute error tolerance for the eigenvalues.
 *          An approximate eigenvalue is accepted as converged
 *          when it is determined to lie in an interval [a,b]
 *          of width less than or equal to
 *
 *                  ABSTOL + EPS *   max( |a|,|b| ) ,
 *
 *          where EPS is the machine precision.  If ABSTOL is less than
 *          or equal to zero, then  EPS*|T|  will be used in its place,
 *          where |T| is the 1-norm of the tridiagonal matrix obtained
 *          by reducing A to tridiagonal form.
 *
 *          See "Computing Small Singular Values of Bidiagonal Matrices
 *          with Guaranteed High Relative Accuracy," by Demmel and
 *          Kahan, LAPACK Working Note #3.
 *
 *          If high relative accuracy is important, set ABSTOL to
 *          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
 *          eigenvalues are computed to high relative accuracy when
 *          possible in future releases.  The current code does not
 *          make any guarantees about high relative accuracy, but
 *          furutre releases will. See J. Barlow and J. Demmel,
 *          "Computing Accurate Eigensystems of Scaled Diagonally
 *          Dominant Matrices", LAPACK Working Note #7, for a discussion
 *          of which matrices define their eigenvalues to high relative
 *          accuracy. (hard set to use tryrac=1).
 *
 * @param[in] nbcomputedeig
 *          The total number of eigenvalues found.  0 <= M <= N.
 *          If RANGE = PlasmaAllVec, M = N, and if RANGE = PlasmaIvec, M = IU-IL+1.
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_ssyevr
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssyevr_Tile
 * @sa PLASMA_ssyevr_Tile_Async
 * @sa PLASMA_cheevr
 * @sa PLASMA_dsyevr
 * @sa PLASMA_ssyevr
 *
 ******************************************************************************/
int PLASMA_ssyevr(PLASMA_enum jobz, PLASMA_enum range, PLASMA_enum uplo,
                 int N,
                 float *A, int LDA,
                 float vl, float vu, int il, int iu,
                 float abstol, int *nbcomputedeig,
                 float *W,
                 PLASMA_desc *descT,
                 float *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_ssyevr", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_ssyevr", "illegal value of jobz");
        return -1;
    }
    if (range != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_ssyevr", "illegal value of jobz");
        return -1;
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_ssyevr", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_ssyevr", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_ssyevr", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_ssyevr", "illegal value of LDQ");
        return -9;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SSYEVD, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssyevr", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_ssyevr_Tile_Async(jobz, range, uplo, &descA, vl, vu, il, iu, abstol, nbcomputedeig, W, descT, Q, LDQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_ssyevr_Tile - Computes all eigenvalues and, optionally, eigenvectors of a
 *  complex Hermitian matrix A using a two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] range
 *          = PlasmaAllVec: all eigenvalues will be found.
 *          = PlasmaVec: all eigenvalues in the half-open interval (VL,VU]
 *                       will be found.
 *          = PlasmaIvec: the IL-th through IU-th eigenvalues will be found.
 *          For range = PlasmaVec or PlasmaIvec and IU - IL < N - 1,
 *          DSTEBZ and ZSTEIN are called.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[in] vl
 *          see vu.
 *
 * @param[in] vu
 *          If RANGE=PlasmaVec, the lower and upper bounds of the
 *          interval to be searched for eigenvalues. VL < VU.
 *          Not referenced if RANGE = PlasmaAllVec or PlasmaIvec.
 *
 * @param[in] il
 *          see iu.
 *
 * @param[in] iu
 *          If RANGE=PlasmaIvec, the indices (in ascending order) of the
 *          smallest and largest eigenvalues to be returned.
 *          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
 *          Not referenced if RANGE = PlasmaAllVec or PlasmaVec.
 *
 * @param[in] abstol
 *          The absolute error tolerance for the eigenvalues.
 *          An approximate eigenvalue is accepted as converged
 *          when it is determined to lie in an interval [a,b]
 *          of width less than or equal to
 *
 *                  ABSTOL + EPS *   max( |a|,|b| ) ,
 *
 *          where EPS is the machine precision.  If ABSTOL is less than
 *          or equal to zero, then  EPS*|T|  will be used in its place,
 *          where |T| is the 1-norm of the tridiagonal matrix obtained
 *          by reducing A to tridiagonal form.
 *
 *          See "Computing Small Singular Values of Bidiagonal Matrices
 *          with Guaranteed High Relative Accuracy," by Demmel and
 *          Kahan, LAPACK Working Note #3.
 *
 *          If high relative accuracy is important, set ABSTOL to
 *          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
 *          eigenvalues are computed to high relative accuracy when
 *          possible in future releases.  The current code does not
 *          make any guarantees about high relative accuracy, but
 *          furutre releases will. See J. Barlow and J. Demmel,
 *          "Computing Accurate Eigensystems of Scaled Diagonally
 *          Dominant Matrices", LAPACK Working Note #7, for a discussion
 *          of which matrices define their eigenvalues to high relative
 *          accuracy. (hard set to use tryrac=1).
 *
 * @param[in] nbcomputedeig
 *          The total number of eigenvalues found.  0 <= M <= N.
 *          If RANGE = PlasmaAllVec, M = N, and if RANGE = PlasmaIvec, M = IU-IL+1.
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_ssyevr
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimention of the eigenvectors matrix Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssyevr_Tile
 * @sa PLASMA_ssyevr_Tile_Async
 * @sa PLASMA_cheevr
 * @sa PLASMA_dsyevr
 * @sa PLASMA_ssyevr
 *
 ******************************************************************************/
int PLASMA_ssyevr_Tile(PLASMA_enum jobz, PLASMA_enum range, PLASMA_enum uplo,
                       PLASMA_desc *A,
                       float vl, float vu, int il, int iu,
                       float abstol, int *nbcomputedeig,
                       float *W, PLASMA_desc *T,
                       float *Q, int LDQ)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssyevr_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_ssyevr_Tile_Async(jobz, range, uplo, A, vl, vu, il, iu, abstol, nbcomputedeig, W, T, Q, LDQ, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_ssyevr_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssyevr
 * @sa PLASMA_ssyevr_Tile
 * @sa PLASMA_cheevr_Tile_Async
 * @sa PLASMA_dsyevr_Tile_Async
 * @sa PLASMA_ssyevr_Tile_Async
 *
 ******************************************************************************/
int PLASMA_ssyevr_Tile_Async(PLASMA_enum jobz, PLASMA_enum range, PLASMA_enum uplo,
                             PLASMA_desc *A,
                             float vl, float vu, int il, int iu,
                             float abstol, int *nbcomputedeig,
                             float *W,
                             PLASMA_desc *T,
                             float *Q, int LDQ,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descQ;
    float *AB;
    float *E, *D;
    int fraction, N, NB, LDAB, status;
    plasma_time_t timelpk, timeaplQ2, timeT;
    plasma_time_t timeB, timeblg, timeaplQ1, timeconv1, timeconv2, timeall;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssyevr_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_ssyevr_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_ssyevr_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "matrix need to be square");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N    = descA.m;
    NB   = min(descA.mb,descA.m);
    LDAB = 2*NB+1;

    /* Allocate workspace for band storage of the band matrix A and for the off diagonal after tridiagonalisation */
    AB = (float *)plasma_shared_alloc(plasma, LDAB*N, PlasmaRealFloat);
    memset( AB, 0, LDAB * N * sizeof(float) );
    if (AB == NULL) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "plasma_shared_alloc(AB) failed");
        plasma_shared_free(plasma, AB);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    D = (float *)plasma_shared_alloc(plasma, N, PlasmaRealDouble);
    E = (float *)plasma_shared_alloc(plasma, N, PlasmaRealDouble);
    if ((E == NULL) || (D == NULL) ) {
        plasma_error("PLASMA_ssyevr_Tile_Async", "plasma_shared_alloc(E) failed");
        plasma_shared_free(plasma, E);
        plasma_shared_free(plasma, D);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    timeB = timeall = plasma_gettime(plasma);

    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */
    /*=======================================
     *  calling Reduction to BAND
     *  then convert matrix to band form
     *=======================================*/
    plasma_dynamic_call_5(plasma_pssytrd_he2hb,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_6( plasma_pssbcpy_t2bl,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        float*, AB,
        int, LDAB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_sync();

    status = sequence->status;
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssyevr","pssytrd_he2hb+pscopy");
        return status;
    }
    /*=======================================
     *  END of calling Reduction to BAND
     *=======================================*/

    timeB = plasma_gettime(plasma) - timeB;
    plasma_printtime("\n  Finish Band red    timing= %lf \n", timeB);

    /*=======================================
     *  calling bulge chasing
     *=======================================*/
    float *TAU2 = NULL;
    float *V2 = NULL;
    float *T2 = NULL;
    int *isuppz = NULL;
    int Vblksiz, blkcnt, LDT, LDV;
    int WANTZ   = 0;
    int blguplo = PlasmaLower;
    /* int NE      = N; // for later use when a portion of the eigenvectors are requested*/
    if( jobz == PlasmaNoVec )
        WANTZ=0;
    else
        WANTZ=2;
    /* Vblksiz correspond to the blocking used when applying V2 to the matrix Q
     * it is similar to IB in LAPACK SORMQR.
     * blkcnt is the number of losange or tile of Vs */
    /* Note that in case PlamaVec requested, the V2 and T2 are stored by the
     * bulgechasing function in a special format:
     * for V2s: it store the V2(LDV,Vblksiz) of each losange in a tile storage meaning
     * that V2_1 is stored then V2_2,..., V2_blkcnt.
     * blkcnt is the number of losange.
     * */
    Vblksiz = min(NB,48);
    LDT     = Vblksiz;
    if( jobz == PlasmaVec ) {
        findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
        TAU2   = (float *) plasma_shared_alloc(plasma,     blkcnt*Vblksiz, PlasmaRealFloat);
        V2     = (float *) plasma_shared_alloc(plasma, LDV*blkcnt*Vblksiz, PlasmaRealFloat);
        T2     = (float *) plasma_shared_alloc(plasma, LDT*blkcnt*Vblksiz, PlasmaRealFloat);
       isuppz  = (int *) plasma_shared_alloc(plasma, 2*N, PlasmaInteger);
        if ( (TAU2 == NULL) || (V2 == NULL) || (T2 == NULL) || (isuppz == NULL)) {
            plasma_error("PLASMA_ssyevr", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            plasma_shared_free(plasma, T2);
            plasma_shared_free(plasma, isuppz);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0,     blkcnt*Vblksiz*sizeof(float));
        memset(V2,   0, LDV*blkcnt*Vblksiz*sizeof(float));
        memset(T2,   0, LDT*blkcnt*Vblksiz*sizeof(float));
    }
    else {
        TAU2   = (float *) plasma_shared_alloc(plasma, 2*N, PlasmaRealFloat);
        V2     = (float *) plasma_shared_alloc(plasma, 2*N, PlasmaRealFloat);
        if ( (TAU2 == NULL) || (V2 == NULL) ) {
            plasma_error("PLASMA_ssyevr", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0, 2*N*sizeof(float));
        memset(V2,   0, 2*N*sizeof(float));
    }

    timeblg = plasma_gettime(plasma);

    plasma_parallel_call_13(plasma_pssytrd_hb2st_v1,
        PLASMA_enum,         blguplo,
        int,                 N,
        int,                 NB,
        int,                 Vblksiz,
        float*, AB,
        int,                 LDAB,
        float*, V2,
        float*, TAU2,
        float*,             D,
        float*,             E,
        int,                 WANTZ,
        PLASMA_sequence*,    sequence,
        PLASMA_request*,     request);
    /* WARNING: If plasma_pssytrd_hb2st is implemented through a dynamic call, don't
     * forget to synchronize */
    plasma_dynamic_sync();
    /*=======================================
     *  END of calling bulge chasing
     *=======================================*/

    timeblg = plasma_gettime(plasma) - timeblg;
    plasma_printtime("  Finish Bulge       timing= %lf \n", timeblg);
    timelpk = plasma_gettime(plasma);

    /*=======================================
     *  calling eigensolver
     *=======================================*/
    plasma_setlapack_multithreads(plasma->world_size);
    /* call eigensolver using lapack routine for our resulting tridiag [D E] */
    /* for internal use */
    /*
    float vl =0.0, vu=0.0;
    int il=1, iu=N;
    int nbcomputedeig=0, NE=N;
    PLASMA_enum range = PlasmaAllVec;
     int fraction = 20;
    if((fraction<100)&&(fraction>0)) {
        range=PlasmaIvec;
        il  = 1;
        iu  = plasma_ceildiv(fraction*N,100);
    }
    */
    int tryrac =1;
    int NE = N;
    if(range==PlasmaIvec) NE = iu-il+1;
    /*Note that for range=V, the size of Q should be Q(LDQ,NZC).
     * right now it is LDQ,N but it could be optimized when NZC becomes a user input.
     * but in this case our interface becomes different from LAPACK.*/
    status = LAPACKE_sstemr( LAPACK_COL_MAJOR, lapack_const(jobz), lapack_const(range),
                             N, D, E, vl, vu, il, iu, nbcomputedeig,
                             W, Q, LDQ, NE, isuppz, &tryrac);
    if(status != 0){
        plasma_error("PLASMA_sstedc","ZSTEDC");
        return status;
    }
    if((nbcomputedeig[0] != NE) && (range != PlasmaAllVec)){
        printf("ZSTEDC failed computed eig %d to equal requested eig %d\n",nbcomputedeig[0],NE);
        return PLASMA_ERR_UNEXPECTED;
    }
    sequence->status = status;
    plasma_setlapack_sequential(plasma);
    /*=======================================
     *  END of calling eigensolver
     *=======================================*/

    timelpk = plasma_gettime(plasma) - timelpk;
    fraction = NE *100 / N;
    plasma_printtime("  Finish Eigensolver timing= %lf  WANTZ %d threads %d fraction %d  NE %d \n" ,timelpk, WANTZ,  plasma->world_size, fraction, NE);

    if (jobz == PlasmaVec){
        /*=======================================
         *  apply Q2 from the bulge
         *=======================================*/
        /* compute T2 */

        timeT = plasma_gettime(plasma);

        plasma_static_call_8(plasma_pslarft_blgtrd,
            int,                 N,
            int,                 NB,
            int,                 Vblksiz,
            float*, V2,
            float*, T2,
            float*, TAU2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeT     = plasma_gettime(plasma) - timeT;
        plasma_printtime("  Finish compute T2  timing= %lf \n", timeT);
        timeaplQ2 = plasma_gettime(plasma);

        /* apply Q2 from Left */
        plasma_static_call_14(plasma_psormqr_blgtrd,
            PLASMA_enum,         PlasmaLeft,
            PLASMA_enum,         PlasmaNoTrans,
            int,                 N,
            int,                 NB,
            int,                 NE,
            int,                 Vblksiz,
            int,                 WANTZ,
            float*, V2,
            float*, T2,
            float*, TAU2,
            float*, Q,
            int,                 LDQ,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeaplQ2 = plasma_gettime(plasma) - timeaplQ2;
        plasma_printtime("  Finish compute Q2  timing= %lf \n", timeaplQ2 );

        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
        /* CASE NB>N, Q1 doesn't need to be applied, only bulge chasing has been done */
        if( NB < N ){
            timeconv1 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_sooplap2tile( descQ, Q, NB, NB, LDQ, NE, 0, 0, N, NE, sequence, request, plasma_desc_mat_free(&(descQ)) );
            }else {
                plasma_siplap2tile( descQ, Q, NB, NB, LDQ, NE, 0, 0, N, NE, sequence, request );
            }

            timeconv1 = plasma_gettime(plasma) - timeconv1;
            timeaplQ1 = plasma_gettime(plasma);

            /* Accumulate the transformations from the first stage */
            if(uplo==PlasmaLower){
                plasma_dynamic_call_7(plasma_psormqr,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n),
                    PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);

            }
            else {
                plasma_dynamic_call_7(plasma_psormlq,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n         ),
                    PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }

            timeaplQ1 = plasma_gettime(plasma) - timeaplQ1;
            plasma_printtime("  Finish compute Q1  timing= %lf \n", timeaplQ1);
            timeconv2 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_sooptile2lap( descQ, Q, NB, NB, LDQ, NE, sequence, request );
                plasma_dynamic_sync();
                plasma_desc_mat_free(&descQ);
            } else {
                plasma_siptile2lap( descQ, Q, NB, NB, LDQ, NE, sequence, request );
                plasma_dynamic_sync();
            }

            timeconv2 = plasma_gettime(plasma) - timeconv2;
            plasma_printtime("  Finish convert     timing= %lf \n", timeconv1+timeconv2);
        } /* END of ( NB < N ) */
    }

    timeall = plasma_gettime(plasma) - timeall;
    plasma_printtime("  Finish full eigenproblem threads %d  N %d  fraction %d  NE %d timeall= %lf \n", plasma->world_size, N, fraction, NE, timeall);

    if( jobz == PlasmaVec ){
        plasma_shared_free(plasma, T2);
        plasma_shared_free(plasma, isuppz);
    }

    plasma_shared_free(plasma, V2);
    plasma_shared_free(plasma, TAU2);
    plasma_shared_free(plasma, E);
    plasma_shared_free(plasma, AB);

    (void)fraction;
    return PLASMA_SUCCESS;
}
