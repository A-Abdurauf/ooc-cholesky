/**
 *
 * @file core_sgbtype2cb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2012-12-15
 * @generated s Fri Apr  1 11:02:31 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

//#define AU(m,n) &(A[(m) + LDA*(n)])
//#define AL(m,n) &(A[(m) + LDA*(n)])
#define AL(m_, n_) (A + NB + LDA * (n_) + ((m_)-(n_)))
#define AU(m_, n_) (A + NB + LDA * (n_) + ((m_)-(n_)+NB))
#define VQ(m)     (VQ + (m))
#define VP(m)     (VP + (m))
#define TAUQ(m)   (TAUQ + (m))
#define TAUP(m)   (TAUP + (m))

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_sgbtype2cb is a kernel that will operate on a region (triangle) of data
 *  bounded by st and ed. This kernel apply the right update remaining from the
 *  type1 and this later will create a bulge so it eliminate the first column of
 *  the created bulge and do the corresponding Left update.
 *
 *  All detail are available on technical report or SC11 paper.
 *  Azzam Haidar, Hatem Ltaief, and Jack Dongarra. 2011.
 *  Parallel reduction to condensed forms for symmetric eigenvalue problems
 *  using aggregated fine-grained and memory-aware kernels. In Proceedings
 *  of 2011 International Conference for High Performance Computing,
 *  Networking, Storage and Analysis (SC '11). ACM, New York, NY, USA, ,
 *  Article 8 , 11 pages.
 *  http://doi.acm.org/10.1145/2063384.2063394
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A.
 *
 * @param[in] NB
 *          The size of the band.
 *
 * @param[in, out] A
 *          A pointer to the matrix A of size (3*NB+1)-by-N.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix A. LDA >= max(1,3*NB+1)
 *
 * @param[in, out] VP
 *          TODO: Check and fix doc
 *          float array, dimension N if eigenvalue only
 *          requested or (LDV*blkcnt*Vblksiz) if Eigenvectors requested
 *          The Householder reflectors of the previous type 1 are used here
 *          to continue update then new one are generated to eliminate the
 *          bulge and stored in this array.
 *
 * @param[in, out] TAUP
 *          TODO: Check and fix doc
 *          float array, dimension (N).
 *          The scalar factors of the Householder reflectors of the previous
 *          type 1 are used here to continue update then new one are generated
 *          to eliminate the bulge and stored in this array.
 *
 * @param[in, out] VQ
 *          TODO: Check and fix doc
 *          float array, dimension N if eigenvalue only
 *          requested or (LDV*blkcnt*Vblksiz) if Eigenvectors requested
 *          The Householder reflectors of the previous type 1 are used here
 *          to continue update then new one are generated to eliminate the
 *          bulge and stored in this array.
 *
 * @param[in, out] TAUQ
 *          TODO: Check and fix doc
 *          float array, dimension (N).
 *          The scalar factors of the Householder reflectors of the previous
 *          type 1 are used here to continue update then new one are generated
 *          to eliminate the bulge and stored in this array.
 *
 * @param[in] st
 *          A pointer to the start index where this kernel will operate.
 *
 * @param[in] ed
 *          A pointer to the end index where this kernel will operate.
 *
 * @param[in] sweep
 *          The sweep number that is eliminated. it serve to calculate the
 *          pointer to the position where to store the Vs and Ts.
 *
 * @param[in] Vblksiz
 *          constant which correspond to the blocking used when applying the Vs.
 *          it serve to calculate the pointer to the position where to store the
 *          Vs and Ts.
 *
 * @param[in] WANTZ
 *          constant which indicate if Eigenvalue are requested or both
 *          Eigenvalue/Eigenvectors.
 *
 * @param[in] WORK
 *          Workspace of size nb.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
/***************************************************************************
 *          TYPE 2-BAND-bidiag Lower/Upper columnwise-Householder
 ***************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgbtype2cb = PCORE_sgbtype2cb
#define CORE_sgbtype2cb PCORE_sgbtype2cb
#endif
void CORE_sgbtype2cb(PLASMA_enum uplo, int N, int NB,
                     float *A, int LDA,
                     float *VQ, float *TAUQ,
                     float *VP, float *TAUP,
                     int st, int ed, int sweep, int Vblksiz, int WANTZ,
                     float *WORK)
{
    float ctmp;
    int i, J1, J2, len, lem, LDX;
    int blkid, vpos, taupos, tpos;

    LDX = LDA-1;
    J1  = ed+1;
    J2  = min(ed+NB,N-1);
    lem = ed-st+1;
    len = J2-J1+1;

    if( uplo == PlasmaUpper ) {
        /* ========================
        *       UPPER CASE
        * ========================*/
        if( len > 0 ) {
            if( WANTZ == 0 ) {
                vpos   = ((sweep+1)%2)*N + st;
                taupos = ((sweep+1)%2)*N + st;
            } else {
                findVTpos(N, NB, Vblksiz, sweep, st,
                          &vpos, &taupos, &tpos, &blkid);
            }
            /* Apply remaining Left commming from type1/3_upper */
            ctmp = (*TAUQ(taupos));
            LAPACKE_slarfx_work(LAPACK_COL_MAJOR, lapack_const(PlasmaLeft),
                                lem, len, VQ(vpos), ctmp, AU(st, J1), LDX, WORK);
        }

        if( len > 1 ) {
            if( WANTZ == 0 ) {
                vpos   = ((sweep+1)%2)*N + J1;
                taupos = ((sweep+1)%2)*N + J1;
            } else {
                findVTpos(N,NB,Vblksiz,sweep,J1, &vpos, &taupos, &tpos, &blkid);
            }

            /* Remove the top row of the created bulge */
            *VP(vpos) = 1.;
            for(i=1; i<len; i++){
                *VP(vpos+i)     = (*AU(st, J1+i));
                *AU(st, J1+i) = 0.;
            }
            /* Eliminate the row  at st */
            ctmp = (*AU(st, J1));
            LAPACKE_slarfg_work(len, &ctmp, VP(vpos+1), 1, TAUP(taupos) );
            *AU(st, J1) = ctmp;
            /*
             * Apply Right on A(J1:J2,st+1:ed)
             * We decrease len because we start at row st+1 instead of st.
             * row st is the row that has been revomved;
             */
            lem = lem-1;
            ctmp = *TAUP(taupos);
            LAPACKE_slarfx_work(LAPACK_COL_MAJOR, lapack_const(PlasmaRight),
                                lem, len, VP(vpos), ctmp, AU(st+1, J1), LDX, WORK);
        }
    }else{
        /* ========================
         *       LOWER CASE
         * ========================*/
        if( len > 0 ) {
            if( WANTZ == 0 ) {
                vpos   = ((sweep+1)%2)*N + st;
                taupos = ((sweep+1)%2)*N + st;
            } else {
                findVTpos(N, NB, Vblksiz, sweep, st,
                          &vpos, &taupos, &tpos, &blkid);
            }
            /* Apply remaining Right commming from type1/3_lower */
            ctmp = (*TAUP(taupos));
            LAPACKE_slarfx_work(LAPACK_COL_MAJOR, lapack_const(PlasmaRight),
                                len, lem, VP(vpos), ctmp, AL(J1, st), LDX, WORK);
        }
        if( len > 1 ) {
            if( WANTZ == 0 ) {
                vpos   = ((sweep+1)%2)*N + J1;
                taupos = ((sweep+1)%2)*N + J1;
            } else {
                findVTpos(N,NB,Vblksiz,sweep,J1, &vpos, &taupos, &tpos, &blkid);
            }

            /* Remove the first column of the created bulge */
            *VQ(vpos) = 1.;
            memcpy(VQ(vpos+1), AL(J1+1, st), (len-1)*sizeof(float));
            memset(AL(J1+1, st), 0, (len-1)*sizeof(float));
            /* Eliminate the col  at st */
            LAPACKE_slarfg_work(len, AL(J1, st), VQ(vpos+1), 1, TAUQ(taupos) );
            /*
             * Apply left on A(J1:J2,st+1:ed)
             * We decrease len because we start at col st+1 instead of st.
             * col st is the col that has been revomved;
             */
            lem = lem-1;
            ctmp = (*TAUQ(taupos));
            LAPACKE_slarfx_work(LAPACK_COL_MAJOR, lapack_const(PlasmaLeft),
                                len, lem, VQ(vpos), ctmp, AL(J1, st+1), LDX, WORK);
        }
    }
    /* end of uplo case */
    return;
}
/***************************************************************************/
#undef AU
#undef AL
#undef VQ
#undef VP
#undef TAUQ
#undef TAUP
