/**
 *
 * @file pcgebrd_gb2bd_v1.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2012-12-15
 * @generated c Fri Apr  1 11:02:59 2016
 *
 **/
#include "common.h"
#undef REAL
#define COMPLEX

/*
 * The storage of the band matrix A is described below:
 * note that we keep same description as Lapack meaning that
 * a lower band start at diagonal (row=0) and below diag are
 * at row i such as
 * For Lower: A(i-j   , j) = A(i,j)
 * For upper: A(NB+i-j, j) = A(i,j)
 *       ----------------------------------
 * NB+1 |               band A             |
 *       ----------------------------------
 * For our reduction from band bidiag to bidiag
 * a bulge is created on both side of the Lower/upper
 * band and thus we allocate NB more on each side.
 *
 * LDAB  = 3*NB+1;
 * AB looks like:
 *       __________________________________
 * NB   |               zero               |
 *       ----------------------------------
 * NB+1 |               band A             |
 *       ----------------------------------
 * NB   |_______________zero_______________|
 *
 * */
/* For Lapack storage */
//#define AU(m,n) &(A[(m) + LDA*(n)])
//#define AL(m,n) &(A[(m) + LDA*(n)])
/* For band storage */
#define AL(m_, n_) (A + NB + LDA * (n_) + ((m_)-(n_)))
#define AU(m_, n_) (A + NB + LDA * (n_) + ((m_)-(n_)+NB))
/***************************************************************************/
/**
 *  Parallel bulge chasing column-wise - static scheduling
 **/
/***************************************************************************/
void plasma_pcgebrd_gb2bd_v1(plasma_context_t *plasma)
{
    int my_core_id = PLASMA_RANK;
    int cores_num  = plasma->world_size;
    /*===========================*/
    int grsiz = 1;
    int MINMN, NB, LDA, Vblksiz, WANTZ;
    float *D;
    float *E;
    PLASMA_enum uplo;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *VQ, *VP;
    PLASMA_Complex32_t *TAUQ, *TAUP;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    /*===========================
     *   local variables
     *===========================*/
    int nbtiles;
    PLASMA_Complex32_t *WORK;
    int coreid, sweepid, myid, shift, stt, st, ed, stind, edind;
    int blklastind, colpt, stepercol;
    int i,j,m,k;
    int thgrsiz, thgrnb, thgrid, thed;
    int colblktile,maxrequiredcores,colpercore,allcoresnb;

    /* variables used when this function is standalone */
#ifdef COMPLEX
    int standalonework = 0 ;
#endif

    plasma_unpack_args_16(uplo, MINMN, NB, Vblksiz, A, LDA, VQ, TAUQ, VP, TAUP, D, E, WANTZ, i, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;


    /* Quick return */
    /*********************
     *    CASE MINMN = 0     *
     *********************
     */
    if (MINMN == 0){
        return;
    }
    /*********************
     *    CASE NB = 0    *
     *********************
     */
    /*
     * The matrix is diagonal. We just copy it (abs for complex) and store it in D
     * sequential code here is better so only core 0 will work
     */
    if ( NB == 0 ) {
        if( my_core_id == 0 ){
            memset(E, 0, (MINMN-1)*sizeof(float));
#ifdef COMPLEX
            if( uplo == PlasmaLower ){
                /* PlasmaLower */
                for (i=0; i<MINMN; i++){
                    D[i]  = cabsf(*AL(i,i));
                }
            }else{
                /* PlasmaUpper */
                for (i=0; i<MINMN; i++){
                    D[i]  = cabsf(*AU(i,i));
                }
            }
#else
            if( uplo == PlasmaLower ){
                /* PlasmaLower */
                for (i=0; i<MINMN; i++){
                    D[i]  = *AL(i,i);
                }
            }else{
                /* PlasmaUpper */
                for (i=0; i<MINMN; i++){
                    D[i]  = *AU(i,i);
                }
            }
#endif
        }
        return;
    }
    /*********************
     *    CASE NB = 1    *
     *********************
     */
    /*
     * The matrix is already Bidiagonal. We do not need to apply the bulge chasing.
     * We make diagonal and superdiagonal elements real, and store them in
     * D and E.
     */
    /* =======================================================
     * NOTE in CASE USED STANDALONE FUNCTION
     * =======================================================*/
    /* in case where the first stage has been done we are sure
     * that all element of E are real so no need to go through
     * making the off diagonal element real just copy them to D and E.
     * However in case where this function is used for only
     * band matrices meaning only stage2 want to be called then it requires
     * to make all off-diagonal elements real and so remove this portion below
     * and uncomment the NB=1 case below.*/

    /* sequential code here so only core 0 will work */
    if ( NB == 1 ) {
        if( my_core_id == 0 ) {

            /* Case used when this function is standalone and so the off diagonal
             * element require to be real. (Only for complex) */
#ifdef COMPLEX
            if( standalonework != 0 ) {
                plasma_error("pcgebrd_gb2bd_v1", "standalone is not supported\n");
                plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
                return;
            }
            else
#endif
            {
                /* Case used when first stage has been applied or all real cases */
                if( uplo == PlasmaLower ){
                    for (i=0; i < MINMN-1; i++) {
                        D[i] = crealf(*AL(i,i));
                        E[i] = crealf(*AL(i+1,i));
                    }
                    D[MINMN-1] = crealf(*AL(MINMN-1,MINMN-1));
                }
                else {
                    for (i=0; i < MINMN-1; i++) {
                        D[i] = crealf(*AU(i,i));
                        E[i] = crealf(*AU(i,i+1));
                    }
                    D[MINMN-1] = crealf(*AU(MINMN-1,MINMN-1));
                }
            }
        }
        return;
    }

    /*********************
     * General case:     *
     *********************
     */
    /*
     * As I store V in the V vector there are overlap between
     * tasks so shift is now 4 where group need to be always
     * multiple of 2 (or shift=5 if not multiple of 2),
     * because as example if grs=1 task 2 from
     * sweep 2 can run with task 6 sweep 1., but task 2 sweep 2
     * will overwrite the V of tasks 5 sweep 1 which are used by
     * task 6, so keep in mind that group need to be multiple of 2,
     * and thus tasks 2 sweep 2 will never run with task 6 sweep 1.
     * OR if we allocate V as V(N,2) and we switch between the storing of
     * sweep's like odd in V(N,1) and even in V(N,2) then no overlap and so
     * shift is 3.
     * when storing V in matrix style, shift could be back to 3.
     * */

    /* Some tunning for the bulge chasing code
     * see technical report for details */
    nbtiles    = plasma_ceildiv(MINMN,NB);
    allcoresnb = cores_num;
    /* grsiz   = 2;*/
    if( WANTZ == 0 ) {
        shift = 3;
    } else {
        shift = 3; /* it was 5 before see above for explanation*/
    }

    if( grsiz == 1 )
        colblktile = 1;
    else
        colblktile = grsiz/2;

    maxrequiredcores = max( nbtiles/colblktile, 1 );
    colpercore = colblktile*NB;
    allcoresnb = min( allcoresnb, maxrequiredcores );
    thgrsiz = MINMN;
    #if defined (ENABLE_DEBUG)
    if(my_core_id==0){
    if(cores_num > maxrequiredcores)
    {
        printf("==================================================================================\n");
        printf("  WARNING only %3d threads are required to run this test optimizing cache reuse\n",maxrequiredcores);
        printf("==================================================================================\n");
    }
    printf("  Static bulgechasing version v9_9col threads  %4d   threads_used  %4d   N %5d      NB %5d    grs %4d thgrsiz %4d  WANTZ %4d\n",cores_num, allcoresnb, MINMN, NB, grsiz,thgrsiz,WANTZ);
    }
    #endif

    /* allocate small workspace for the kernel to avoid an alloc/dealloc at each call */
    WORK    = (PLASMA_Complex32_t*) plasma_private_alloc(plasma, NB, PlasmaComplexFloat);
    /* Initialize static scheduler progress table */
    ss_init(2*nbtiles+shift+cores_num+10, 1, 0);

    /* main bulge chasing code */
    i = shift/grsiz;
    stepercol =  i*grsiz == shift ? i:i+1;
    i       = (MINMN-1)/thgrsiz;
    thgrnb  = i*thgrsiz == (MINMN-1) ? i:i+1;
    for (thgrid = 1; thgrid<=thgrnb; thgrid++){
        stt  = (thgrid-1)*thgrsiz+1;
        thed = min( (stt + thgrsiz -1), (MINMN-1));
        for (i = stt; i <= MINMN-1; i++){
            ed = min(i,thed);
            if(stt>ed) break;
            for (m = 1; m <=stepercol; m++){
                st=stt;
                for (sweepid = st; sweepid <=ed; sweepid++){

                    for (k = 1; k <=grsiz; k++){
                        myid = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + k;
                        if(myid%2 ==0){
                            colpt      = (myid/2)*NB+1+sweepid-1;
                            stind      = colpt-NB+1;
                            edind      = min(colpt,MINMN);
                            blklastind = colpt;
                        } else {
                            colpt      = ((myid+1)/2)*NB + 1 +sweepid -1 ;
                            stind      = colpt-NB+1;
                            edind      = min(colpt,MINMN);
                            if( (stind>=edind-1) && (edind==MINMN) )
                                blklastind=MINMN;
                            else
                                blklastind=0;
                        }
                        coreid = (stind/colpercore)%allcoresnb;

                        if(my_core_id==coreid) {
                            if(myid==1) {

                                ss_cond_wait(myid+shift-1, 0, sweepid-1);
                                CORE_cgbtype1cb(uplo, MINMN, NB, A, LDA, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);
                                ss_cond_set(myid, 0, sweepid);

                                if(blklastind >= (MINMN-1)) {
                                    for (j = 1; j <= shift; j++)
                                        ss_cond_set(myid+j, 0, sweepid);
                                }
                            } else {
                                ss_cond_wait(myid-1,       0, sweepid);
                                ss_cond_wait(myid+shift-1, 0, sweepid-1);
                                if(myid%2 == 0){
                                    CORE_cgbtype2cb(uplo, MINMN, NB, A, LDA, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);
                                }else{
                                    CORE_cgbtype3cb(uplo, MINMN, NB, A, LDA, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);
                                }
                                ss_cond_set(myid, 0, sweepid);
                                if(blklastind >= (MINMN-1)) {
                                    for (j = 1; j <= shift+allcoresnb; j++)
                                        ss_cond_set(myid+j, 0, sweepid);
                                }
                            } /* END if myid==1 */
                        } /* END if my_core_id==coreid */

                        if(blklastind >= (MINMN-1)) {
                            stt++;
                            break;
                        }
                    } /* END for k=1:grsiz */
                } /* END for sweepid=st:ed */
            } /* END for m=1:stepercol */
        } /* END for i=1:MINMN-1 */
    } /* END for thgrid=1:thgrnb */

    /* finalize static sched */
    ss_finalize();

    /*================================================
     *  store resulting diag and lower diag D and E
     *  note that D and E are always real after the bulgechasing
     *================================================*/
    /* sequential code here so only core 0 will work */
    if( my_core_id == 0 ) {
        if( uplo == PlasmaLower ){
            for (i=0; i < MINMN-1; i++) {
                D[i] = crealf(*AL(i,i));
                E[i] = crealf(*AL(i+1,i));
            }
            D[MINMN-1] = crealf(*AL(MINMN-1,MINMN-1));
        }
        else {
            for (i=0; i < MINMN-1; i++) {
                D[i] = crealf(*AU(i,i));
                E[i] = crealf(*AU(i,i+1));
            }
            D[MINMN-1] = crealf(*AU(MINMN-1,MINMN-1));
        }
    }

    plasma_private_free(plasma, WORK);
    return;
} /* END FUNCTION */
/***************************************************************************/

/***************************************************************************/
/**
 *  Parallel bulge chasing Reduction from BAND Bidiagonal
 *  column-wise - dynamic scheduling
 **/
/***************************************************************************/
void plasma_pcgebrd_gb2bd_v1_quark(PLASMA_enum uplo, int MINMN, int NB, int Vblksiz,
                                   PLASMA_Complex32_t *A, int LDA,
                                   PLASMA_Complex32_t *VQ, PLASMA_Complex32_t *TAUQ,
                                   PLASMA_Complex32_t *VP, PLASMA_Complex32_t *TAUP,
                                   float *D, float *E, int WANTZ, int WANTP,
                                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    /*===========================
     *   local variables
     *===========================*/
    int INgrsiz, INthgrsiz;
    int myid, grsiz, shift=3, stt, st, ed, stind, edind;
    int blklastind, colpt, PCOL, ACOL, MCOL;
    int stepercol, mylastid, grnb, grid;
    int *DEP,*MAXID;
    int i, sweepid, m;
    int thgrsiz, thgrnb, thgrid, thed;
    /* variables used when this function is standalone */
#ifdef COMPLEX
    int standalonework = 0 ;
#endif

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    /* Quick return */
    if (MINMN == 0) {
        return;
    }
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    /*
     * The matrix is diagonal. We just copy it (abs for complex) and store it in D
     * sequential code here is better so only core 0 will work
     */
    if ( NB == 0 ) {
        memset(E, 0, (MINMN-1)*sizeof(float));
#ifdef COMPLEX
            if( uplo == PlasmaLower ){
                /* PlasmaLower */
                for (i=0; i<MINMN; i++){
                    D[i]  = cabsf(*AL(i,i));
                }
            }else{
                /* PlasmaUpper */
                for (i=0; i<MINMN; i++){
                    D[i]  = cabsf(*AU(i,i));
                }
            }
#else
            if( uplo == PlasmaLower ){
                /* PlasmaLower */
                for (i=0; i<MINMN; i++){
                    D[i]  = *AL(i,i);
                }
            }else{
                /* PlasmaUpper */
                for (i=0; i<MINMN; i++){
                    D[i]  = *AU(i,i);
                }
            }
#endif
        return;
    }
    /*
     * Barrier is used because the bulge have to wait until
     * the reduction to band has been finish.
     * otherwise, I can remove this BARRIER when I integrate
     * the function dependencies link inside the reduction to
     * band. Keep in mind the case when NB=1, where no bulge-chasing.
     */
    /***************************************************************/
    QUARK_Barrier(plasma->quark);
    /***************************************************************/
    /*
     * Case NB=1 ==> matrix is already Bidiagonal. no need to bulge.
     * We make diagonal and superdiagonal elements real, and store them in
     * D and E.
     * For Q, PT: ZSCAL should be done if WANTQ == 1.
     */
    /* =======================================================
     * NOTE in CASE USED STANDALONE FUNCTION
     * =======================================================*/
    /* in case where the first stage has been done we are sure
     * that all element of E are real so no need to go through
     * making the off diagonal element real just copy them to D and E.
     * However in case where this function is used for only
     * band matrices meaning only stage2 want to be called then it requires
     * to make all off-diagonal elements real and so remove this portion below
     * and uncomment the NB=1 case below.*/
    /* sequential code here so only core 0 will work */
    if ( NB == 1 ) {
        /* Case used when this function is standalone and so the off diagonal
         * element require to be real. (Only for complex) */
#ifdef COMPLEX
        if( standalonework != 0 ) {
                plasma_error("pcgebrd_gb2bd_v1", "standalone is not supported\n");
                plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
                return;
        }
        else
#endif
            {
                /* Case used when first stage has been applied or all real cases */
                if( uplo == PlasmaLower ){
                    for (i=0; i < MINMN-1; i++) {
                        D[i] = crealf(*AL(i,i));
                        E[i] = crealf(*AL(i+1,i));
                    }
                    D[MINMN-1] = crealf(*AL(MINMN-1,MINMN-1));
                }
                else {
                    for (i=0; i < MINMN-1; i++) {
                        D[i] = crealf(*AU(i,i));
                        E[i] = crealf(*AU(i,i+1));
                    }
                    D[MINMN-1] = crealf(*AU(MINMN-1,MINMN-1));
                }
            }
        return;
    }
    /*********************
     * General case:     *
     *********************
     */
    /***************************************************************************
     *                       START BULGE CHASING CODE
     **************************************************************************/
    /* General case NB > 1 && N > NB */
    DEP   = (int *)                plasma_shared_alloc(plasma, MINMN+1, PlasmaInteger      );
    MAXID = (int *)                plasma_shared_alloc(plasma, MINMN+1, PlasmaInteger      );
    memset(MAXID,0,(MINMN+1)*sizeof(int));
    QUARK_Barrier(plasma->quark);

    /*
     * Initialisation of local parameter. those parameter should be
     * input or tuned parameter.
     */
    INgrsiz = 1;
    if( NB > 160 ) {
        INgrsiz = 2;
    }
    else if( NB > 100 ) {
        if( MINMN < 5000 )
            INgrsiz = 2;
        else
            INgrsiz = 4;
    } else {
        INgrsiz = 6;
    }
    INthgrsiz = MINMN;

    grsiz   = INgrsiz;
    thgrsiz = INthgrsiz;
    if( grsiz   == 0 ) grsiz   = 6;
    if( thgrsiz == 0 ) thgrsiz = MINMN;

    i = shift/grsiz;
    stepercol =  i*grsiz == shift ? i:i+1;

    i       = (MINMN-1)/thgrsiz;
    thgrnb  = i*thgrsiz == (MINMN-1) ? i:i+1;

    #if defined (ENABLE_DEBUG)
    printf("  Dynamic bulgechasing version v9_9col threads  %4d   threads_used  %4d   N %5d      NB %5d    grs %4d thgrsiz %4d  WANTZ %4d\n",plasma->world_size, plasma->world_size, MINMN, NB, grsiz,thgrsiz,WANTZ);
    #endif

    for (thgrid = 1; thgrid<=thgrnb; thgrid++){
        stt  = (thgrid-1)*thgrsiz+1;
        thed = min( (stt + thgrsiz -1), (MINMN-1));
        for (i = stt; i <= MINMN-1; i++){
            ed=min(i,thed);
            if(stt>ed)break;
            for (m = 1; m <=stepercol; m++){
                st=stt;
                for (sweepid = st; sweepid <=ed; sweepid++){
                    /* PCOL:  dependency on the ID of the master of the group of the previous column.  (Previous Column:PCOL). */
                    /* ACOL:  dependency on the ID of the master of the previous group of my column.   (Acctual  Column:ACOL). (it is 0(NULL) for myid=1) */
                    /* MCOL:  OUTPUT dependency on the my ID, to be used by the next ID. (My Column: MCOL). I am the master of this group. */
                    myid     = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + 1;
                    mylastid = myid+grsiz-1;
                    PCOL     = mylastid+shift-1;  /* to know the dependent ID of the previous column. need to know the master of its group*/
                    MAXID[sweepid] = myid;
                    PCOL     = min(PCOL,MAXID[sweepid-1]); /* for the last columns, we might do only 1 or 2 kernel, so the PCOL will be wrong. this is to force it to the last ID of the previous col.*/
                    grnb     = PCOL/grsiz;
                    grid     = grnb*grsiz == PCOL ? grnb:grnb+1;
                    PCOL     = (grid-1)*grsiz +1; /* give me the ID of the master of the group of the previous column.*/
                    ACOL     = myid-grsiz;
                    if(myid==1)ACOL=0;
                    MCOL     = myid;

                    QUARK_CORE_cbrdalg1(
                        plasma->quark, &task_flags,
                        uplo, MINMN, NB, A, LDA,
                        VQ, TAUQ, VP, TAUP,
                        Vblksiz, WANTZ,
                        i, sweepid, m, grsiz,
                        &DEP[PCOL], &DEP[ACOL], &DEP[MCOL] );

                    if(mylastid%2 ==0){
                        blklastind      = (mylastid/2)*NB+1+sweepid-1;
                    }else{
                        colpt      = ((mylastid+1)/2)*NB + 1 +sweepid -1 ;
                        stind      = colpt-NB+1;
                        edind      = min(colpt,MINMN);
                        if( (stind>=edind-1) && (edind==MINMN) )
                            blklastind=MINMN;
                        else
                            blklastind=0;
                    }
                    if(blklastind >= (MINMN-1))  stt=stt+1;
                } /* END for sweepid=st:ed    */
            } /* END for m=1:stepercol */
        } /* END for i=1:MINMN-2      */
    } /* END for thgrid=1:thgrnb     */

     /*
     * Barrier used only for now, to be sure that everything
     * is done before copying the D and E and free workspace.
     * this will be removed later when D and E are directly filled
     * during the bulge process.
     */
    QUARK_Barrier(plasma->quark);
    plasma_shared_free(plasma, (void*) DEP);
    plasma_shared_free(plasma, (void*) MAXID);
    /*================================================
     *  store resulting diag and lower diag D and E
     *  note that D and E are always real after the bulgechasing
     *================================================*/
    /* sequential code here so only core 0 will work */
    memset(D, 0, MINMN*sizeof(float));
    memset(E, 0, (MINMN-1)*sizeof(float));
    if( uplo == PlasmaLower ){
        for (i=0; i < MINMN-1; i++) {
            D[i] = crealf(*AL(i,i));
            E[i] = crealf(*AL(i+1,i));
        }
        D[MINMN-1] = crealf(*AL(MINMN-1,MINMN-1));
    }
    else {
        for (i=0; i < MINMN-1; i++) {
            D[i] = crealf(*AU(i,i));
            E[i] = crealf(*AU(i,i+1));
        }
        D[MINMN-1] = crealf(*AU(MINMN-1,MINMN-1));
    }
    return;
} /* END FUNCTION */
/***************************************************************************/
#undef AL
#undef AU
