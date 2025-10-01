/**
 *
 * @file dgecfi2.h
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
 * @generated d Fri Apr  1 11:02:55 2016
 *
 **/

#ifndef DGECFI2_H
#define DGECFI2_H

#define ipt_call( name, m1, n1, mb, nb ) \
  ipt_d##name(plasma, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_d##name(plasma, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_d##name(plasma, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_d##name(plasma, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

#define ipt_cal2( name, m1, n1, mb, nb ) \
  ipt_d##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_d##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_d##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_d##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

/* one transformation */
#define ipt_drm2rrrb(  plasma, m, n, A, mb, nb, seq, req) ipt_dcm2ccrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_drrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_dccrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_dcm2ccrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dccrb2cm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_dccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dcrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

#define ipt_dcrrb2rrrb(plasma, m, n, A, mb, nb, seq, req) ipt_dccrb2rcrb((plasma), (m), (n), (A), (mb), (nb), (seq), (req));
#define ipt_drcrb2ccrb(plasma, m, n, A, mb, nb, seq, req) ipt_dccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_drrrb2crrb(plasma, m, n, A, mb, nb, seq, req) ipt_dccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_dccrb2rcrb(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 2 transformations */
#define ipt_drm2crrb(  plasma, m, n, A, mb, nb, seq, req) ipt_dcm2rcrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_dcrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_drcrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_dcm2rcrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drcrb2cm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_dccrb2rrrb(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drrrb2ccrb(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dcrrb2rcrb(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drcrb2crrb(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_dcm2crrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dcrrb2cm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drcrb2rm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drm2rcrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 3 transformations */
int ipt_dcm2rrrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drrrb2cm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dccrb2rm  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drm2ccrb  (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 4 transformations */
int ipt_dcm2rm    (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_drm2cm    (plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);


int ipt_dpanel2all(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dall2panel(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dpanel2tile(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_dtile2panel(plasma_context_t *plasma, int m, int n, double *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
#endif /* DGECFI2_H */
