/**
 *
 * @file sgecfi2.h
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated s Fri Apr  1 11:02:55 2016
 *
 **/

#ifndef SGECFI2_H
#define SGECFI2_H

#define ipt_call( name, m1, n1, mb, nb ) \
  ipt_s##name(plasma, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_s##name(plasma, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_s##name(plasma, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_s##name(plasma, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

#define ipt_cal2( name, m1, n1, mb, nb ) \
  ipt_s##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_s##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_s##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_s##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

/* one transformation */
#define ipt_srm2rrrb(  plasma, m, n, A, mb, nb, seq, req) ipt_scm2ccrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_srrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_sccrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_scm2ccrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_sccrb2cm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_sccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_scrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

#define ipt_scrrb2rrrb(plasma, m, n, A, mb, nb, seq, req) ipt_sccrb2rcrb((plasma), (m), (n), (A), (mb), (nb), (seq), (req));
#define ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, seq, req) ipt_sccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_srrrb2crrb(plasma, m, n, A, mb, nb, seq, req) ipt_sccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_sccrb2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 2 transformations */
#define ipt_srm2crrb(  plasma, m, n, A, mb, nb, seq, req) ipt_scm2rcrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_scrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_srcrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_scm2rcrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srcrb2cm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_sccrb2rrrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srrrb2ccrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_scrrb2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srcrb2crrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_scm2crrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_scrrb2cm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srcrb2rm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srm2rcrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 3 transformations */
int ipt_scm2rrrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srrrb2cm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_sccrb2rm  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srm2ccrb  (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 4 transformations */
int ipt_scm2rm    (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_srm2cm    (plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);


int ipt_spanel2all(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_sall2panel(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_spanel2tile(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_stile2panel(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
#endif /* SGECFI2_H */
