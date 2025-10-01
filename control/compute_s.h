/**
 *
 * @file compute_s.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:03:01 2016
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_sdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_sooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = plasma_desc_init(                                           \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pslapack_to_tile,                                        \
        float*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

#define plasma_siplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = plasma_desc_init(                                         \
        PlasmaRealFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_sgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), (seq), (req));


#define plasma_sooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    plasma_parallel_call_5(plasma_pstile_to_lapack,                 \
                           PLASMA_desc,         (descA),            \
                           float*, (A),                \
                           int,                 (lm),               \
                           PLASMA_sequence*,    (seq),              \
                           PLASMA_request*,     (req));

#define plasma_siptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    PLASMA_sgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), (seq), (req));


#define plasma_sooplap2tile_noalloc( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    plasma_parallel_call_5(                                             \
        plasma_pslapack_to_tile,                                        \
        float*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pstradd  (plasma_context_t *plasma);
void plasma_psgelqf (plasma_context_t *plasma);
void plasma_psgemm  (plasma_context_t *plasma);
void plasma_psgeqrf (plasma_context_t *plasma);
void plasma_psgerbb (plasma_context_t *plasma);
void plasma_psgetmi2(plasma_context_t *plasma);
void plasma_psgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pssymm  (plasma_context_t *plasma);
void plasma_pssyrk  (plasma_context_t *plasma);
void plasma_pssyr2k (plasma_context_t *plasma);
#endif
void plasma_pslacpy (plasma_context_t *plasma);
void plasma_pslag2d (plasma_context_t *plasma);
void plasma_pslange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pslansy (plasma_context_t *plasma);
#endif
void plasma_pslansy (plasma_context_t *plasma);
void plasma_pspack  (plasma_context_t *plasma);
void plasma_psplgsy (plasma_context_t *plasma);
void plasma_psplgsy (plasma_context_t *plasma);
void plasma_pspltmg(plasma_context_t *plasma);
void plasma_pspotrf (plasma_context_t *plasma);
void plasma_psshift (plasma_context_t *plasma);
void plasma_pssymm  (plasma_context_t *plasma);
void plasma_pssyrk  (plasma_context_t *plasma);
void plasma_pssyr2k (plasma_context_t *plasma);
void plasma_pstrmm  (plasma_context_t *plasma);
void plasma_pstrsm  (plasma_context_t *plasma);
void plasma_pstrsmpl(plasma_context_t *plasma);
void plasma_pstrsmrv(plasma_context_t *plasma);
void plasma_psorglq (plasma_context_t *plasma);
void plasma_psorgqr (plasma_context_t *plasma);
void plasma_psorgqrrh(plasma_context_t *plasma);
void plasma_psormlq (plasma_context_t *plasma);
void plasma_psormqr (plasma_context_t *plasma);
void plasma_psunpack(plasma_context_t *plasma);
void plasma_psgebrd_gb2bd_v1(plasma_context_t *plasma);
void plasma_pssytrd_hb2st_v1(plasma_context_t *plasma);
void plasma_psormqr_blgtrd(plasma_context_t *plasma);
void plasma_pslarft_blgtrd(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int plasma_sshift(plasma_context_t *plasma, int m, int n, float *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_spltmg_condex(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_spltmg_house( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
#ifdef REAL
void plasma_pslaed0_quark(int icompq, char range, int id1, int n, float *D, float *E, int il, int iu, float vl, float vu, float *Q, int LDQ, float *qstore, int lwork, float *WORK, float *WORK2, int LDWORK, int *IWORK, int *localdata, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaed1_quark(char range, int n, int n1, float *D, int il, int iu, float vl, float vu, float *Q, int LDQ, int *INDXQ, float *beta, float *work, float *work2, int *iwork, int *K1, int *K2, int last_merge, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
#ifdef COMPLEX
void plasma_pslag2c_quark(int m, int n, const float *R, int ldr, float *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pstradd_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_tl2row_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psbarrier_row2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgebrd_gb2bd_v1_quark(PLASMA_enum uplo, int MINMN, int NB, int Vblksiz,
                                   float *A, int LDA,
                                   float *VQ, float *TAUQ,
                                   float *VP, float *TAUP,
                                   float *D, float *E, int WANTZ, int WANTP,
                                   PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgebrd_ge2gb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgemm_quark(PLASMA_enum transA, PLASMA_enum transB, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgeqp3_quark( PLASMA_desc A, int *jpvt, float *tau,
                           float *work, float *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_psgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, float *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_nopiv_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV, PLASMA_desc W, int *Wi, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, float *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psgbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, float *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssygst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pssymm_quark(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pssytrd_hb2st_v1_quark(PLASMA_enum uplo, int N, int NB, int Vblksiz, float *A, int LDA, float *V, float *TAU, float *D, float *E, int WANTZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssytrd_he2hb_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslag2d_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslange_quark(PLASMA_enum norm, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pslansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pslansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslantr_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslascal_quark( PLASMA_enum uplo, float alpha, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaset_quark( PLASMA_enum uplo, float alpha, float beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaset2_quark(PLASMA_enum uplo, float alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaset_identity_quark(int n, float *A, int lda, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaswp_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaswpc_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psplgsy_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psplgsy_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_quark(PLASMA_enum mtxtype, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_fiedler_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_toeppd_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_circul_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_chebvand_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspltmg_hankel_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pspotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psshift_quark(int, int, int, float *, int *, int, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_psstedc_quark(PLASMA_enum compz, int matsiz, float *D, float *E, float *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psswaps_quark(int n, int *perm, float *Z, int LDZ, float *work, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssymm_quark(PLASMA_enum side, PLASMA_enum uplo, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pssyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, const int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrsmrv_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, float alpha, PLASMA_desc A, PLASMA_desc W, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psorgtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_psormlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);


/* Subroutines to be added in lapacke */
//float PLASMA_FCALL(slanst,SLANST)(char *type, int *n, float *D, float *E);
//void   PLASMA_FCALL(slascl,SLASCL)(char *type, int *i1, int *i2,
//                                   float *scale, float *scale2,
//                                   int *n, int *i4, float *D, int *m, int *info);
//void   PLASMA_FCALL(slamrg,SLAMRG)(int *K, int *n2, float *D, int *id1, int *id2, int *INDXQ);
//
//#ifdef REAL
//float PLASMA_FCALL(slapy2,ZLAPY2)(float *f1, float *f2);
//float PLASMA_FCALL(slamc3,ZLAMC3)(float *f1, float *f2);
//void               PLASMA_FCALL(slaed4,ZLAED4)(const int *n, const int *i, const float *D,
//                                               const float *Z, float *delta, const float *beta,
//                                               float *Di, int *info);
//#endif

