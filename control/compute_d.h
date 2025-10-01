/**
 *
 * @file compute_d.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:03:01 2016
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_ddesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_dooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = plasma_desc_init(                                           \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pdlapack_to_tile,                                        \
        double*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

#define plasma_dooplap2tile_affinity( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = plasma_desc_init(                                           \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( plasma_desc_mat_alloc_affinity( &(descA) , (plasma) ) ) {                 \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pdlapack_to_tile_affinity,                                   \
        double*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

#define plasma_diplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = plasma_desc_init(                                         \
        PlasmaRealDouble, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_dgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), (seq), (req));


#define plasma_dooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    plasma_parallel_call_5(plasma_pdtile_to_lapack,                 \
                           PLASMA_desc,         (descA),            \
                           double*, (A),                \
                           int,                 (lm),               \
                           PLASMA_sequence*,    (seq),              \
                           PLASMA_request*,     (req));

#define plasma_dooptile2lap_affinity( descA, A, mb, nb, lm, ln, seq, req)    \
    plasma_parallel_call_5(plasma_pdtile_to_lapack_affinity,                 \
                           PLASMA_desc,         (descA),            \
                           double*, (A),                \
                           int,                 (lm),               \
                           PLASMA_sequence*,    (seq),              \
                           PLASMA_request*,     (req));

#define plasma_diptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    PLASMA_dgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), (seq), (req));


#define plasma_dooplap2tile_noalloc( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    plasma_parallel_call_5(                                             \
        plasma_pdlapack_to_tile,                                        \
        double*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pdtradd  (plasma_context_t *plasma);
void plasma_pdgelqf (plasma_context_t *plasma);
void plasma_pdgemm  (plasma_context_t *plasma);
void plasma_pdgeqrf (plasma_context_t *plasma);
void plasma_pdgerbb (plasma_context_t *plasma);
void plasma_pdgetmi2(plasma_context_t *plasma);
void plasma_pdgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pdsymm  (plasma_context_t *plasma);
void plasma_pdsyrk  (plasma_context_t *plasma);
void plasma_pdsyr2k (plasma_context_t *plasma);
#endif
void plasma_pdlacpy (plasma_context_t *plasma);
void plasma_pdlag2s (plasma_context_t *plasma);
void plasma_pdlange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pdlansy (plasma_context_t *plasma);
#endif
void plasma_pdlansy (plasma_context_t *plasma);
void plasma_pdpack  (plasma_context_t *plasma);
void plasma_pdplgsy (plasma_context_t *plasma);
void plasma_pdplgsy (plasma_context_t *plasma);
void plasma_pdpltmg(plasma_context_t *plasma);
void plasma_pdpotrf (plasma_context_t *plasma);
#ifdef PLASMA_WITH_CUDA
void plasma_pdpotrf_gpu (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_dist (plasma_context_t *plasma);
#endif
#ifdef PLASMA_WITH_HPCASIA24
#ifdef PLASMA_WITH_CUDA
void plasma_pdpotrf_gpu_async_copy (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table_all_managed (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table_no_free (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_static_table (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_dist (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table_dist (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table_all_managed_dist (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_table_no_free_dist (plasma_context_t *plasma);
void plasma_pdpotrf_gpu_reuse_data_static_table_dist (plasma_context_t *plasma);
#endif
#endif
void plasma_pdshift (plasma_context_t *plasma);
void plasma_pdsymm  (plasma_context_t *plasma);
void plasma_pdsyrk  (plasma_context_t *plasma);
void plasma_pdsyr2k (plasma_context_t *plasma);
void plasma_pdtrmm  (plasma_context_t *plasma);
void plasma_pdtrsm  (plasma_context_t *plasma);
void plasma_pdtrsmpl(plasma_context_t *plasma);
void plasma_pdtrsmrv(plasma_context_t *plasma);
void plasma_pdorglq (plasma_context_t *plasma);
void plasma_pdorgqr (plasma_context_t *plasma);
void plasma_pdorgqrrh(plasma_context_t *plasma);
void plasma_pdormlq (plasma_context_t *plasma);
void plasma_pdormqr (plasma_context_t *plasma);
void plasma_pdunpack(plasma_context_t *plasma);
void plasma_pdgebrd_gb2bd_v1(plasma_context_t *plasma);
void plasma_pdsytrd_hb2st_v1(plasma_context_t *plasma);
void plasma_pdormqr_blgtrd(plasma_context_t *plasma);
void plasma_pdlarft_blgtrd(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int plasma_dshift(plasma_context_t *plasma, int m, int n, double *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_dpltmg_condex(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_dpltmg_house( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
#ifdef REAL
void plasma_pdlaed0_quark(int icompq, char range, int id1, int n, double *D, double *E, int il, int iu, double vl, double vu, double *Q, int LDQ, double *qstore, int lwork, double *WORK, double *WORK2, int LDWORK, int *IWORK, int *localdata, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaed1_quark(char range, int n, int n1, double *D, int il, int iu, double vl, double vu, double *Q, int LDQ, int *INDXQ, double *beta, double *work, double *work2, int *iwork, int *K1, int *K2, int last_merge, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
#ifdef COMPLEX
void plasma_pdlag2z_quark(int m, int n, const double *R, int ldr, double *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pdtradd_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_tl2row_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdbarrier_row2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, double *D, double *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgebrd_gb2bd_v1_quark(PLASMA_enum uplo, int MINMN, int NB, int Vblksiz,
                                   double *A, int LDA,
                                   double *VQ, double *TAUQ,
                                   double *VP, double *TAUP,
                                   double *D, double *E, int WANTZ, int WANTP,
                                   PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgebrd_ge2gb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgemm_quark(PLASMA_enum transA, PLASMA_enum transB, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgeqp3_quark( PLASMA_desc A, int *jpvt, double *tau,
                           double *work, double *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, double *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_nopiv_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV, PLASMA_desc W, int *Wi, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, double *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdgbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, double *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, double *D, double *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsygst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pdsymm_quark(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pdsytrd_hb2st_v1_quark(PLASMA_enum uplo, int N, int NB, int Vblksiz, double *A, int LDA, double *V, double *TAU, double *D, double *E, int WANTZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsytrd_he2hb_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlag2s_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlange_quark(PLASMA_enum norm, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pdlansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pdlansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlantr_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlascal_quark( PLASMA_enum uplo, double alpha, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaset_quark( PLASMA_enum uplo, double alpha, double beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaset2_quark(PLASMA_enum uplo, double alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaset_identity_quark(int n, double *A, int lda, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaswp_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlaswpc_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdplgsy_quark(double bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdplgsy_quark(double bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_quark(PLASMA_enum mtxtype, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_fiedler_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_toeppd_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_circul_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_chebvand_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpltmg_hankel_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pdpotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef PLASMA_WITH_CUDA
void plasma_pdpotrf_gpu_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
#ifdef PLASMA_WITH_HPCASIA24
#ifdef PLASMA_WITH_CUDA
void plasma_pdpotrf_gpu_async_copy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_all_managed_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_no_free_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_static_table_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_all_managed_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_table_no_free_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdpotrf_gpu_reuse_data_static_table_dist_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
#endif
void plasma_pdshift_quark(int, int, int, double *, int *, int, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pdstedc_quark(PLASMA_enum compz, int matsiz, double *D, double *E, double *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdswaps_quark(int n, int *perm, double *Z, int LDZ, double *work, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsymm_quark(PLASMA_enum side, PLASMA_enum uplo, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, const int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrsmrv_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, double alpha, PLASMA_desc A, PLASMA_desc W, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdorgtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdormlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);


/* Subroutines to be added in lapacke */
//double PLASMA_FCALL(dlanst,DLANST)(char *type, int *n, double *D, double *E);
//void   PLASMA_FCALL(dlascl,DLASCL)(char *type, int *i1, int *i2,
//                                   double *scale, double *scale2,
//                                   int *n, int *i4, double *D, int *m, int *info);
//void   PLASMA_FCALL(dlamrg,DLAMRG)(int *K, int *n2, double *D, int *id1, int *id2, int *INDXQ);
//
//#ifdef REAL
//double PLASMA_FCALL(dlapy2,ZLAPY2)(double *f1, double *f2);
//double PLASMA_FCALL(dlamc3,ZLAMC3)(double *f1, double *f2);
//void               PLASMA_FCALL(dlaed4,ZLAED4)(const int *n, const int *i, const double *D,
//                                               const double *Z, double *delta, const double *beta,
//                                               double *Di, int *info);
//#endif

