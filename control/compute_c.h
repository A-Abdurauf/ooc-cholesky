/**
 *
 * @file compute_c.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:03:01 2016
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_cdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_cooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = plasma_desc_init(                                           \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pclapack_to_tile,                                        \
        PLASMA_Complex32_t*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

#define plasma_ciplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = plasma_desc_init(                                         \
        PlasmaComplexFloat, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_cgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), (seq), (req));


#define plasma_cooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    plasma_parallel_call_5(plasma_pctile_to_lapack,                 \
                           PLASMA_desc,         (descA),            \
                           PLASMA_Complex32_t*, (A),                \
                           int,                 (lm),               \
                           PLASMA_sequence*,    (seq),              \
                           PLASMA_request*,     (req));

#define plasma_ciptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    PLASMA_cgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), (seq), (req));


#define plasma_cooplap2tile_noalloc( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    plasma_parallel_call_5(                                             \
        plasma_pclapack_to_tile,                                        \
        PLASMA_Complex32_t*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pctradd  (plasma_context_t *plasma);
void plasma_pcgelqf (plasma_context_t *plasma);
void plasma_pcgemm  (plasma_context_t *plasma);
void plasma_pcgeqrf (plasma_context_t *plasma);
void plasma_pcgerbb (plasma_context_t *plasma);
void plasma_pcgetmi2(plasma_context_t *plasma);
void plasma_pcgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pchemm  (plasma_context_t *plasma);
void plasma_pcherk  (plasma_context_t *plasma);
void plasma_pcher2k (plasma_context_t *plasma);
#endif
void plasma_pclacpy (plasma_context_t *plasma);
void plasma_pclag2z (plasma_context_t *plasma);
void plasma_pclange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pclanhe (plasma_context_t *plasma);
#endif
void plasma_pclansy (plasma_context_t *plasma);
void plasma_pcpack  (plasma_context_t *plasma);
void plasma_pcplghe (plasma_context_t *plasma);
void plasma_pcplgsy (plasma_context_t *plasma);
void plasma_pcpltmg(plasma_context_t *plasma);
void plasma_pcpotrf (plasma_context_t *plasma);
void plasma_pcshift (plasma_context_t *plasma);
void plasma_pcsymm  (plasma_context_t *plasma);
void plasma_pcsyrk  (plasma_context_t *plasma);
void plasma_pcsyr2k (plasma_context_t *plasma);
void plasma_pctrmm  (plasma_context_t *plasma);
void plasma_pctrsm  (plasma_context_t *plasma);
void plasma_pctrsmpl(plasma_context_t *plasma);
void plasma_pctrsmrv(plasma_context_t *plasma);
void plasma_pcunglq (plasma_context_t *plasma);
void plasma_pcungqr (plasma_context_t *plasma);
void plasma_pcungqrrh(plasma_context_t *plasma);
void plasma_pcunmlq (plasma_context_t *plasma);
void plasma_pcunmqr (plasma_context_t *plasma);
void plasma_pcunpack(plasma_context_t *plasma);
void plasma_pcgebrd_gb2bd_v1(plasma_context_t *plasma);
void plasma_pchetrd_hb2st_v1(plasma_context_t *plasma);
void plasma_pcunmqr_blgtrd(plasma_context_t *plasma);
void plasma_pclarft_blgtrd(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions
 **/
int plasma_cshift(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_cpltmg_condex(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_cpltmg_house( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
#ifdef REAL
void plasma_pslaed0_quark(int icompq, char range, int id1, int n, float *D, float *E, int il, int iu, float vl, float vu, float *Q, int LDQ, float *qstore, int lwork, float *WORK, float *WORK2, int LDWORK, int *IWORK, int *localdata, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslaed1_quark(char range, int n, int n1, float *D, int il, int iu, float vl, float vu, float *Q, int LDQ, int *INDXQ, float *beta, float *work, float *work2, int *iwork, int *K1, int *K2, int last_merge, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
#ifdef COMPLEX
void plasma_pslag2c_quark(int m, int n, const float *R, int ldr, PLASMA_Complex32_t *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pctradd_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_Complex32_t beta, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_tl2row_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcbarrier_row2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgbrdb_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgebrd_gb2bd_v1_quark(PLASMA_enum uplo, int MINMN, int NB, int Vblksiz,
                                   PLASMA_Complex32_t *A, int LDA,
                                   PLASMA_Complex32_t *VQ, PLASMA_Complex32_t *TAUQ,
                                   PLASMA_Complex32_t *VP, PLASMA_Complex32_t *TAUP,
                                   float *D, float *E, int WANTZ, int WANTP,
                                   PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgebrd_ge2gb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgemm_quark(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgeqp3_quark( PLASMA_desc A, int *jpvt, PLASMA_Complex32_t *tau,
                           PLASMA_Complex32_t *work, float *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, PLASMA_Complex32_t *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_nopiv_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV, PLASMA_desc W, int *Wi, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_Complex32_t *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcgbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_Complex32_t *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchbrdt_quark(PLASMA_enum uplo, PLASMA_desc A, float *D, float *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchegst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pchemm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcherk_quark(PLASMA_enum uplo, PLASMA_enum trans, float alpha, PLASMA_desc A, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcher2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, float beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pchetrd_hb2st_v1_quark(PLASMA_enum uplo, int N, int NB, int Vblksiz, PLASMA_Complex32_t *A, int LDA, PLASMA_Complex32_t *V, PLASMA_Complex32_t *TAU, float *D, float *E, int WANTZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pchetrd_he2hb_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclag2z_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclange_quark(PLASMA_enum norm, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pclanhe_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pclansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclantr_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, float *work, float *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclascal_quark( PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaset_quark( PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_Complex32_t beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaset2_quark(PLASMA_enum uplo, PLASMA_Complex32_t alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaset_identity_quark(int n, PLASMA_Complex32_t *A, int lda, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaswp_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclaswpc_quark(PLASMA_desc B, const int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcplghe_quark(float bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcplgsy_quark(PLASMA_Complex32_t bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_quark(PLASMA_enum mtxtype, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_fiedler_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_toeppd_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_circul_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_chebvand_quark( PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpltmg_hankel_quark( PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pcpotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcshift_quark(int, int, int, PLASMA_Complex32_t *, int *, int, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pcstedc_quark(PLASMA_enum compz, int matsiz, float *D, float *E, PLASMA_Complex32_t *Z, int LDZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcswaps_quark(int n, int *perm, PLASMA_Complex32_t *Z, int LDZ, PLASMA_Complex32_t *work, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcsymm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_Complex32_t beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex32_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, const int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrsmrv_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc W, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcungtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pcunmlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);


/* Subroutines to be added in lapacke */
//float PLASMA_FCALL(clanst,CLANST)(char *type, int *n, float *D, float *E);
//void   PLASMA_FCALL(clascl,CLASCL)(char *type, int *i1, int *i2,
//                                   float *scale, float *scale2,
//                                   int *n, int *i4, PLASMA_Complex32_t *D, int *m, int *info);
//void   PLASMA_FCALL(clamrg,CLAMRG)(int *K, int *n2, float *D, int *id1, int *id2, int *INDXQ);
//
//#ifdef REAL
//PLASMA_Complex32_t PLASMA_FCALL(clapy2,ZLAPY2)(PLASMA_Complex32_t *f1, PLASMA_Complex32_t *f2);
//PLASMA_Complex32_t PLASMA_FCALL(clamc3,ZLAMC3)(PLASMA_Complex32_t *f1, PLASMA_Complex32_t *f2);
//void               PLASMA_FCALL(claed4,ZLAED4)(const int *n, const int *i, const float *D,
//                                               const float *Z, float *delta, const float *beta,
//                                               float *Di, int *info);
//#endif

