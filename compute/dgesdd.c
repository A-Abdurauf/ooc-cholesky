/**
 *
 * @file dgesdd.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:56 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

#define REAL
/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dgesdd - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors, by using divide-and-conquer method. The SVD is written
 *
 *       A = U * SIGMA * ugate-transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
 *  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns VT = V**T, not V.
 *
 *  The divide and conquer algorithm makes very mild assumptions about
 *  floating point arithmetic. It will work on machines with a guard
 *  digit in add/subtract, or on those binary machines without guard
 *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
 *  without guard digits, but we know of none.
 *
 *  NOTE that the input parameter of this function differ from the
 *  LAPACK DGESDD. We have jobu and jobvt instead of only jobz first
 *  to be consistent with DGESVD and second to allow a computation
 *  of only U or VT.
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec   = 'A'(lapack):  all M columns of U are returned
 *                        in array U;
 *          = PlasmaNoVec = 'N':  no columns of U (no left singular vectors)
 *                        are computed.
 *          = PlasmaSVec  = 'S': the first min(m,n) columns of U (the left
 *                        singular vectors) are returned in the array U;
 *                        NOT SUPPORTTED YET
 *          = PlasmaOVec  = 'O': the first min(m,n) columns of U (the left
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**T.
 *          Intended usage:
 *          = PlasmaVec   = 'A'(lapack): all N rows of V**T are returned
 *                        in the array VT;
 *          = PlasmaNoVec = 'N': no rows of V**T (no right singular vectors)
 *                        are computed.
 *          = PlasmaSVec  = 'S': the first min(m,n) rows of V**T (the right
 *                        singular vectors) are returned in the array VT;
 *                        NOT SUPPORTTED YET
 *          = PlasmaOVec  = 'O': the first min(m,n) rows of V**T (the right
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 *          Note: jobu and jobvt cannot both be PlasmaOVec.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**T (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] S
 *          The double precision singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_dgesdd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**T;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**T (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgesdd_Tile
 * @sa PLASMA_dgesdd_Tile_Async
 * @sa PLASMA_cgesdd
 * @sa PLASMA_dgesdd
 * @sa PLASMA_sgesdd
 *
 ******************************************************************************/
int PLASMA_dgesdd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N,
                  double *A, int LDA,
                  double *S,
                  PLASMA_desc *descT,
                  double *U, int LDU,
                  double *VT, int LDVT)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgesdd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobu != PlasmaNoVec  && jobu !=PlasmaVec) {
        plasma_error("PLASMA_dgesdd", "illegal value of jobu");
        return -1;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_dgesdd", "illegal value of jobvt");
        return -2;
    }
    if (jobvt != jobu) {
        plasma_error("PLASMA_dgesdd", "in this version: jobu should be equal jobvt");
        return -22;
    }
    if (M < 0) {
        plasma_error("PLASMA_dgesdd", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_dgesdd", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_dgesdd", "illegal value of LDA");
        return -6;
    }
    if (LDU < 1) {
        plasma_error("PLASMA_dgesdd", "illegal value of LDU");
        return -9;
    }
    if (LDVT < 1) {
        plasma_error("PLASMA_dgesdd", "illegal value of LDVT");
        return -11;
    }
    /* Quick return */
    if (min(M, N) == 0) {
        return PLASMA_SUCCESS;
    }


    /* Tune NB & IB depending on M & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DGESVD, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgesdd", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N,
                            sequence, &request);
    }



    /* Call the tile interface */
    PLASMA_dgesdd_Tile_Async(jobu, jobvt, &descA, S, descT, U, LDU, VT, LDVT, sequence, &request);


    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }


    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dgesdd_Tile - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Tile equivalent of PLASMA_dgesdd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec   = 'A'(lapack):  all M columns of U are returned
 *                        in array U;
 *          = PlasmaNoVec = 'N':  no columns of U (no left singular vectors)
 *                        are computed.
 *          = PlasmaSVec  = 'S': the first min(m,n) columns of U (the left
 *                        singular vectors) are returned in the array U;
 *                        NOT SUPPORTTED YET
 *          = PlasmaOVec  = 'O': the first min(m,n) columns of U (the left
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**T.
 *          Intended usage:
 *          = PlasmaVec   = 'A'(lapack): all N rows of V**T are returned
 *                        in the array VT;
 *          = PlasmaNoVec = 'N': no rows of V**T (no right singular vectors)
 *                        are computed.
 *          = PlasmaSVec  = 'S': the first min(m,n) rows of V**T (the right
 *                        singular vectors) are returned in the array VT;
 *                        NOT SUPPORTTED YET
 *          = PlasmaOVec  = 'O': the first min(m,n) rows of V**T (the right
 *                        singular vectors) are overwritten on the array A;
 *                        NOT SUPPORTTED YET
 *
 *          Note: jobu and jobvt cannot both be PlasmaOVec.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**T (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[out] S
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[in, out] T
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_dgesdd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**T;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**T (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 *******************************************************************************
 *
 * @return
 *          \return PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dgesdd
 * @sa PLASMA_dgesdd_Tile_Async
 * @sa PLASMA_cgesdd_Tile
 * @sa PLASMA_dgesdd_Tile
 * @sa PLASMA_sgesdd_Tile
 *
 ******************************************************************************/
int PLASMA_dgesdd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt,
                       PLASMA_desc *A,
                       double *S,
                       PLASMA_desc *T,
                       double *U, int LDU,
                       double *VT, int LDVT)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgesdd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dgesdd_Tile_Async(jobu, jobvt, A, S, T, U, LDU, VT, LDVT, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dgesdd_Tile_Async - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Non-blocking equivalent of PLASMA_dgesdd_Tile().
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
 * @sa PLASMA_dgesdd
 * @sa PLASMA_dgesdd_Tile
 * @sa PLASMA_cgesdd_Tile_Async
 * @sa PLASMA_dgesdd_Tile_Async
 * @sa PLASMA_sgesdd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dgesdd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt,
                             PLASMA_desc *A,
                             double *S,
                             PLASMA_desc *T,
                             double *U, int LDU,
                             double *VT, int LDVT,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descU, descVT;
    double *AB;
    double *E;
    int M, N, MINMN, NB, LDAB, i, status, fraction;
    plasma_time_t timelpk, timeaplQ2, timeT;
    plasma_time_t timeB, timeblg, timeaplQ1, timeconv1, timeconv2, timeall;

    plasma_context_t *plasma;
    plasma = plasma_context_self();

    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dgesdd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dgesdd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dgesdd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "invalid fourth descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (jobu != PlasmaNoVec  && jobu != PlasmaVec) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "illegal value of jobu");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "illegal value of jobvt");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (jobvt != jobu) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "in this version: jobu should be equal jobvt");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    timeall = plasma_gettime(plasma);

    PLASMA_enum uplo = descA.m >= descA.n ? PlasmaUpper : PlasmaLower;
    M     = descA.m;
    N     = descA.n;
    MINMN = min(M,N);
    NB    = min(descA.mb,MINMN);
    LDAB  = 3*NB+1;
    /*=======================================
     *  case M<NB or N<NB call lapack
     *=======================================*/
    if( ( M<= NB) || ( N <=NB ) ){
        /* convert the tile descA to lapack A */
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            AB = (double *) plasma_shared_alloc(plasma, M*N, PlasmaRealDouble);
            if (AB == NULL) {
                plasma_error("PLASMA_dgesvd_Tile_Async", "plasma_shared_alloc(AB-0-) failed");
                plasma_shared_free(plasma, AB);
                return PLASMA_ERR_OUT_OF_RESOURCES;
            }
            plasma_dooptile2lap( descA, AB, NB, NB, M, N,  sequence, request);
        } else {
            AB = descA.mat;
            plasma_diptile2lap( descA, AB, NB, NB, M, N,  sequence, request);
        }
        plasma_dynamic_sync();
        /*=======================================
         *  calling LAPACK DGESVD
         *=======================================*/
        timelpk = plasma_gettime(plasma);

        plasma_setlapack_multithreads(plasma->world_size);
        /* call SVD solver using lapack routine  */
        PLASMA_enum jobz  = jobu  == PlasmaVec ? PlasmaAllVec : PlasmaNoVec;
        status = LAPACKE_dgesdd(LAPACK_COL_MAJOR, lapack_const(jobz),
                                M, N, AB, M, S, U, LDU, VT, LDVT);
        if(status != 0){
            plasma_error("PLASMA_dgesvd","DGESVD");
        }
        sequence->status = status;
        plasma_setlapack_sequential(plasma);

        timelpk = plasma_gettime(plasma)-timelpk;
        plasma_printtime("  Finish Eigensolver-lpkonly timing= %lf  threads %d\n",timelpk, plasma->world_size);

        /*=======================================
         *  END of calling SVD solver
         *=======================================*/
        /* convert the lapack to the tile descA */
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            //plasma_dooplap2tile_noalloc( descA, AB, NB, NB,  M, N, 0, 0, M, N, sequence, request);
            plasma_parallel_call_5( plasma_pdlapack_to_tile,
            double*, AB,
            int,                 M,
            PLASMA_desc,         descA,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
            plasma_dynamic_sync();
            free(AB);
        } else {
            plasma_diplap2tile( descA, AB, NB, NB,  M, N, 0, 0, M, N,
                                sequence, request);
        }
        plasma_dynamic_sync();
        return PLASMA_SUCCESS;
    }
    /*=======================================
     *  END OF case M<NB or N<NB
     *=======================================*/
    /*
     * Allocate workspace for band storage of the band matrix A
     * AB looks like:
     *       __________________________________
     * NB   |               zero               |
     *       ----------------------------------
     * NB+1 |               band A             |
     *       ----------------------------------
     * NB   |_______________zero_______________|
     *
     * */
    AB = (double *)plasma_shared_alloc(plasma, LDAB*MINMN, PlasmaRealDouble);
    memset( AB, 0, LDAB * MINMN * sizeof(double) );
    if (AB == NULL) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "plasma_shared_alloc(AB) failed");
        plasma_shared_free(plasma, AB);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    E = (double *)plasma_shared_alloc(plasma, MINMN, PlasmaRealDouble);
    if (E == NULL) {
        plasma_error("PLASMA_dgesdd_Tile_Async", "plasma_shared_alloc(E) failed");
        plasma_shared_free(plasma, E);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    /*=======================================
     *  calling Reduction from DENSE to BAND
     *  then convert matrix to band form
     *=======================================*/
    timeB = plasma_gettime(plasma);

    /*
     * Reduction to BAND bidiagonal form
     * May be further optimized using the algo described in Trefethen
     */
    /* if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
        plasma_dynamic_call_4(plasma_pdgebrd_ge2gb,
            PLASMA_desc, descA,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    /* } */
    /* else { */
    /*     plasma_dynamic_call_4(plasma_pdgebrd_ge2gb_rh, */
    /*         PLASMA_desc, descA, */
    /*         PLASMA_desc, descT, */
    /*         PLASMA_sequence*, sequence, */
    /*         PLASMA_request*, request); */
    /* } */
    //plasma_dynamic_sync();
    plasma_dynamic_call_6( plasma_pdgbcpy_t2bl,
        PLASMA_enum, descA.m >= descA.n ? PlasmaUpper : PlasmaLower,
        PLASMA_desc, descA,
        double*, &(AB[NB]),
        int, LDAB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);
    plasma_dynamic_sync();
    status = sequence->status;
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dgesdd","pdgebrd_ge2gb+pdcopy");
        return status;
    }

    timeB = plasma_gettime(plasma)-timeB;
    plasma_printtime("\n  Finish Band red    timing= %lf \n",timeB);

    /*=======================================
     *  END of calling Reduction to BAND
     *=======================================*/
    /*=======================================
     *  calling Reduction from BAND to bidiag
     *=======================================*/
    double *VQ2   = NULL;
    double *VP2   = NULL;
    double *TAUQ2 = NULL;
    double *TAUP2 = NULL;
    double *TQ2   = NULL;
    double *TP2   = NULL;
    int Vblksiz, blkcnt, LDT, LDV;
    int WANTZ   = 0;
    /* compute fraction of singular vectors */
    /* for future use when a portion of the eigenvectors are requested */
    int ME      = M;
    int NE      = N;
    int MINMNE  = MINMN;
    /*
    int fraction = 100;
    if((fraction<100)&&(fraction>0)) {
        MINMNE  = plasma_ceildiv(fraction*MINMN,100);
        ME      = MINMNE;
        NE      = MINMNE;
    }
    */
    if( jobu == PlasmaNoVec )
        WANTZ=0;
    else
        WANTZ=2;
    /* Vblksiz correspond to the blocking used when applying V2 to the matrix U
     * it is similar to IB in LAPACK DORMQR.
     * blkcnt is the number of diamond or tile of Vs */
    /* Note that in case PlamaVec requested, the V2 and T2 are stored by the
     * bulgechasing function in a special format:
     * for V2s: it store the V2(LDV,Vblksiz) of each diamond in a tile storage meaning
     * that V2_1 is stored then V2_2,..., V2_blkcnt.
     * blkcnt is the number of diamond.
     * */
    Vblksiz = min(NB,48);
    LDT     = Vblksiz;
    /* data for U */
    if( jobu == PlasmaVec ) {
        findVTsiz(MINMN, NB, Vblksiz, &blkcnt, &LDV);
        TAUQ2 = (double *) plasma_shared_alloc(plasma,     blkcnt*Vblksiz, PlasmaRealDouble);
        VQ2   = (double *) plasma_shared_alloc(plasma, LDV*blkcnt*Vblksiz, PlasmaRealDouble);
        TQ2   = (double *) plasma_shared_alloc(plasma, LDT*blkcnt*Vblksiz, PlasmaRealDouble);
        if ( (TAUQ2 == NULL) || (VQ2 == NULL) || (TQ2 == NULL) ) {
            plasma_error("PLASMA_dgesdd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAUQ2);
            plasma_shared_free(plasma, VQ2);
            plasma_shared_free(plasma, TQ2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAUQ2, 0,     blkcnt*Vblksiz*sizeof(double));
        memset(VQ2,   0, LDV*blkcnt*Vblksiz*sizeof(double));
        memset(TQ2,   0, LDT*blkcnt*Vblksiz*sizeof(double));
    }
    else {
        TAUQ2 = (double *) plasma_shared_alloc(plasma, 2*MINMN, PlasmaRealDouble);
        VQ2   = (double *) plasma_shared_alloc(plasma, 2*MINMN, PlasmaRealDouble);
        if ( (TAUQ2 == NULL) || (VQ2 == NULL) ) {
            plasma_error("PLASMA_dgesdd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAUQ2);
            plasma_shared_free(plasma, VQ2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAUQ2, 0, 2*MINMN*sizeof(double));
        memset(VQ2,   0, 2*MINMN*sizeof(double));
    }
    /* data for VT */
    if( jobvt == PlasmaVec ) {
        findVTsiz(MINMN, NB, Vblksiz, &blkcnt, &LDV);
        TAUP2 = (double *) plasma_shared_alloc(plasma,     blkcnt*Vblksiz, PlasmaRealDouble);
        VP2   = (double *) plasma_shared_alloc(plasma, LDV*blkcnt*Vblksiz, PlasmaRealDouble);
        TP2   = (double *) plasma_shared_alloc(plasma, LDT*blkcnt*Vblksiz, PlasmaRealDouble);
        if ( (TAUP2 == NULL) || (VP2 == NULL) || (TP2 == NULL) ) {
            plasma_error("PLASMA_dgesdd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAUP2);
            plasma_shared_free(plasma, VP2);
            plasma_shared_free(plasma, TP2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAUP2, 0,     blkcnt*Vblksiz*sizeof(double));
        memset(VP2,   0, LDV*blkcnt*Vblksiz*sizeof(double));
        memset(TP2,   0, LDT*blkcnt*Vblksiz*sizeof(double));
    }
    else {
        TAUP2 = (double *) plasma_shared_alloc(plasma, 2*MINMN, PlasmaRealDouble);
        VP2   = (double *) plasma_shared_alloc(plasma, 2*MINMN, PlasmaRealDouble);
        if ( (TAUP2 == NULL) || (VP2 == NULL) ) {
            plasma_error("PLASMA_dgesdd", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAUQ2);
            plasma_shared_free(plasma, VQ2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAUP2, 0, 2*MINMN*sizeof(double));
        memset(VP2,   0, 2*MINMN*sizeof(double));
    }
    /*=======================================
     *  calling bulge chasing
     *=======================================*/

    timeblg = plasma_gettime(plasma);

    plasma_parallel_call_16(plasma_pdgebrd_gb2bd_v1,
        PLASMA_enum,         uplo,
        int,                 MINMN,
        int,                 NB,
        int,                 Vblksiz,
        double*, AB,
        int,                 LDAB,
        double*, VQ2,
        double*, TAUQ2,
        double*, VP2,
        double*, TAUP2,
        double*,             S,
        double*,             E,
        int,                 WANTZ,
        int,                 WANTZ,
        PLASMA_sequence*,    sequence,
        PLASMA_request*,     request);
    /* WARNING: If plasma_pdgebrd_gb2bd is implemented through a dynamic call, don't
     * forget to synchronize */
    plasma_dynamic_sync();

    timeblg = plasma_gettime(plasma)-timeblg;
    plasma_printtime("  Finish Bulge       timing= %lf \n",timeblg);

    /*=======================================
     *  END of calling bulge chasing
     *=======================================*/
    /*=======================================
     *  calling SVD solver
     *=======================================*/
    plasma_setlapack_multithreads(plasma->world_size);
    /* call SVD solver using lapack routine for our resulting bidiag [S E] */
    /* call eigensolver using lapack routine for our resulting tridiag [D E] */
    if((jobu == PlasmaNoVec)&&(jobvt == PlasmaNoVec)){
        timelpk = plasma_gettime(plasma);
        status = LAPACKE_dbdsdc(LAPACK_COL_MAJOR, lapack_const(uplo), lapack_const(PlasmaNoVec),
                                MINMN, S, E, NULL, LDU, NULL, LDVT, NULL, NULL);
        timelpk = plasma_gettime(plasma)-timelpk;
    } else {
        // set U and VT to zero
        memset(VT,   0, N*LDVT*sizeof(double));
        memset(U,   0, M*LDU*sizeof(double));
        /* Initialize U(N+1:M,N+1:M) and VT(M+1:N,M+1:N) to Identity */

        for(i=N; i<M; i++){
            U[i+i*LDU] = 1.0;
        }
         for(i=M; i<N; i++){
            VT[i+i*LDVT] = 1.0;
        }
#if defined COMPLEX
        int j;
        double *RU    = NULL;
        double *RVT   = NULL;
        RU   = (double *) plasma_shared_alloc(plasma,     MINMN*MINMN, PlasmaRealDouble);
        RVT  = (double *) plasma_shared_alloc(plasma,     MINMN*MINMN, PlasmaRealDouble);

        timelpk   = plasma_gettime(plasma);
        status = LAPACKE_dbdsdc(LAPACK_COL_MAJOR, lapack_const(uplo), lapack_const(PlasmaIvec),
                                MINMN, S, E, RU, MINMN, RVT, MINMN, NULL, NULL);
        timelpk = plasma_gettime(plasma)-timelpk;
        // copy the real U and VT matrix to the complex
        //LAPACKE_dlacp2(LAPACK_COL_MAJOR, lapack_const(uplo), MINMN, MINMN, RU, MINMN, U, LDU );
        //LAPACKE_dlacp2(LAPACK_COL_MAJOR, lapack_const(uplo), MINMN, MINMN, RVT, MINMN, VT, LDVT );
        for(j=0; j<MINMN; j++){
            for(i=0; i<MINMN; i++){
                U[i+j*LDU] = (double) RU[i+j*MINMN];
            }
        }
        for(j=0; j<MINMN; j++){
            for(i=0; i<MINMN; i++){
                VT[i+j*LDVT] = (double) RVT[i+j*MINMN];
            }
        }
        plasma_shared_free(plasma, RU);
        plasma_shared_free(plasma, RVT);
#else
        timelpk   = plasma_gettime(plasma);
        status = LAPACKE_dbdsdc(LAPACK_COL_MAJOR, lapack_const(uplo), lapack_const(PlasmaIvec),
                                MINMN, S, E, U, LDU, VT, LDVT, NULL, NULL);
        timelpk = plasma_gettime(plasma)-timelpk;
#endif
    }
    if(status != 0){
        plasma_error("PLASMA_dbdedc","ZBDEDC");
        /*return status;*/
    }
    sequence->status = status;
    plasma_setlapack_sequential(plasma);
    fraction = (MINMNE*100)/N;
    plasma_printtime("  Finish Eigensolver timing= %lf  WANTZ %d threads %d fraction %d  MINMNE %d \n" ,timelpk, WANTZ,  plasma->world_size, fraction, MINMNE);

    /*=======================================
     *  END of calling SVD solver
     *=======================================*/
    /*=======================================
     *  generate U from the bulge
     *=======================================*/
    if (jobu == PlasmaVec){
        timeT = plasma_gettime(plasma);

        /* compute T2 */
        plasma_static_call_8(plasma_pdlarft_blgtrd,
            int,                 MINMN,
            int,                 NB,
            int,                 Vblksiz,
            double*, VQ2,
            double*, TQ2,
            double*, TAUQ2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeT = plasma_gettime(plasma)-timeT;
        plasma_printtime("  Finish compute TU2  timing= %lf \n",timeT);
        timeaplQ2 = plasma_gettime(plasma);

        /* apply Q2 from Left */
        plasma_static_call_14(plasma_pdormqr_blgtrd,
            PLASMA_enum,         PlasmaLeft,
            PLASMA_enum,         PlasmaNoTrans,
            int,                 MINMN,
            int,                 NB,
            int,                 MINMNE,
            int,                 Vblksiz,
            int,                 WANTZ,
            double*, VQ2,
            double*, TQ2,
            double*, TAUQ2,
            double*, U,
            int,                 LDU,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeaplQ2 = plasma_gettime(plasma)-timeaplQ2;
        plasma_printtime("  Finish compute U2  timing= %lf \n",timeaplQ2);

        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
        /* CASE NB>N, Q1 doesn't need to be applied, only bulge chasing has been done */
        if( NB < N ){
            timeconv1 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooplap2tile( descU, U, NB, NB, LDU, ME, 0, 0, M, ME, sequence, request, plasma_desc_mat_free(&(descU)) );
            } else {
                plasma_diplap2tile( descU, U, NB, NB, LDU, ME, 0, 0, M, ME, sequence, request);
            }

            timeconv1 = plasma_gettime(plasma)-timeconv1;
            timeaplQ1 = plasma_gettime(plasma);

            /* Accumulate the transformations from the first stage */
            if(M<N){
                plasma_dynamic_call_7(plasma_pdormqr,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descU, descU.mb, 0, descU.m-descU.mb, descU.n),
                    PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }
            else {
                plasma_dynamic_call_7(plasma_pdormqr,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, descA,
                    PLASMA_desc, descU,
                    PLASMA_desc, descT,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }

            timeaplQ1 = plasma_gettime(plasma)-timeaplQ1;
            plasma_printtime("  Finish compute U1  timing= %lf \n",timeaplQ1);
            timeconv2 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooptile2lap( descU, U, NB, NB, LDU, ME, sequence, request );
                plasma_dynamic_sync();
                plasma_desc_mat_free(&descU);
            } else {
                plasma_diptile2lap( descU, U, NB, NB, LDU, ME, sequence, request );
                plasma_dynamic_sync();
            }

            timeconv2 = plasma_gettime(plasma)-timeconv2;
            plasma_printtime("  Finish convert  U  timing= %lf \n",timeconv1+timeconv2);
        } /* END of ( NB < N ) */
    }
    /*=======================================
     *  END of calling computing U
     *=======================================*/
    /*=======================================
     *  generate VT from the bulge
     *=======================================*/
    if (jobvt == PlasmaVec){
        timeT = plasma_gettime(plasma);

        /* compute T2 */
        plasma_static_call_8(plasma_pdlarft_blgtrd,
            int,                 MINMN,
            int,                 NB,
            int,                 Vblksiz,
            double*, VP2,
            double*, TP2,
            double*, TAUP2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeT = plasma_gettime(plasma)-timeT;
        plasma_printtime("  Finish compute TV2  timing= %lf \n",timeT);
        timeaplQ2 = plasma_gettime(plasma);

        /* apply Q2 from Left */
        plasma_static_call_14(plasma_pdormqr_blgtrd,
            PLASMA_enum,         PlasmaRight,
            PLASMA_enum,         PlasmaTrans,
            int,                 MINMN,
            int,                 NB,
            int,                 MINMNE,
            int,                 Vblksiz,
            int,                 WANTZ,
            double*, VP2,
            double*, TP2,
            double*, TAUP2,
            double*, VT,
            int,                 LDVT,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);

        timeaplQ2 = plasma_gettime(plasma)-timeaplQ2;
        plasma_printtime("  Finish compute V2 timing= %lf \n",timeaplQ2);

        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
        /* CASE NB>N, Q1 doesn't need to be applied, only bulge chasing has been done */
        if( NB < N ){
            timeconv1 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, NE, N, sequence, request, plasma_desc_mat_free(&(descVT)) );
            } else {
                plasma_diplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, NE, N, sequence, request);
            }

            timeconv1 = plasma_gettime(plasma)-timeconv1;
            timeaplQ1 = plasma_gettime(plasma);

            /* Accumulate the transformations from the first stage */
            if(M<N){
                plasma_dynamic_call_7(plasma_pdormlq,
                    PLASMA_enum, PlasmaRight,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, descA,
                    PLASMA_desc, descVT,
                    PLASMA_desc, descT,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }
            else {
                plasma_dynamic_call_7(plasma_pdormlq,
                    PLASMA_enum, PlasmaRight,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descVT,0, descVT.nb, descVT.m, descVT.n-descVT.nb),
                    PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }

            timeaplQ1 = plasma_gettime(plasma)-timeaplQ1;
            plasma_printtime("  Finish compute V1  timing= %lf \n",timeaplQ1);
            timeconv2 = plasma_gettime(plasma);

            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooptile2lap( descVT, VT, NB, NB, LDVT, N, sequence, request );
                plasma_dynamic_sync();
                plasma_desc_mat_free(&descVT);
            } else {
                plasma_diptile2lap( descVT, VT, NB, NB, LDVT, N, sequence, request );
                plasma_dynamic_sync();
            }

            timeconv2 = plasma_gettime(plasma)-timeconv2;
            plasma_printtime("  Finish convert VT  timing= %lf \n",timeconv1+timeconv2);
        } /* END of ( NB < N ) */
    }
    /*=======================================
     *  END of calling computing VT
     *=======================================*/
    timeall = plasma_gettime(plasma)-timeall;
    plasma_printtime("  Finish full eigenproblem threads %d  N %d  fraction %d  MINMNE %d  timeall= %lf \n", plasma->world_size, MINMN, fraction, MINMNE, timeall);

    if( jobu  == PlasmaVec )
        plasma_shared_free(plasma, TQ2);
    if( jobvt == PlasmaVec )
        plasma_shared_free(plasma, TP2);
    plasma_shared_free(plasma, VQ2);
    plasma_shared_free(plasma, TAUQ2);
    plasma_shared_free(plasma, VP2);
    plasma_shared_free(plasma, TAUP2);
    plasma_shared_free(plasma, E);
    plasma_shared_free(plasma, AB);

    (void)fraction;
    return PLASMA_SUCCESS;
}
