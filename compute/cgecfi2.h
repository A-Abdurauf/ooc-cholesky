/**
 *
 * @file cgecfi2.h
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
 * @generated c Fri Apr  1 11:02:55 2016
 *
 **/

#ifndef CGECFI2_H
#define CGECFI2_H

#define ipt_call( name, m1, n1, mb, nb ) \
  ipt_c##name(plasma, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_c##name(plasma, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_c##name(plasma, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_c##name(plasma, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

#define ipt_cal2( name, m1, n1, mb, nb ) \
  ipt_c##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_c##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_c##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_c##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

/* one transformation */
#define ipt_crm2rrrb(  plasma, m, n, A, mb, nb, seq, req) ipt_ccm2ccrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_crrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_cccrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_ccm2ccrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_cccrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_cccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_ccrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

#define ipt_ccrrb2rrrb(plasma, m, n, A, mb, nb, seq, req) ipt_cccrb2rcrb((plasma), (m), (n), (A), (mb), (nb), (seq), (req));
#define ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, seq, req) ipt_cccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_crrrb2crrb(plasma, m, n, A, mb, nb, seq, req) ipt_cccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_cccrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 2 transformations */
#define ipt_crm2crrb(  plasma, m, n, A, mb, nb, seq, req) ipt_ccm2rcrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_ccrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_crcrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_ccm2rcrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crcrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_cccrb2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crrrb2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_ccrrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crcrb2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_ccm2crrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_ccrrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crcrb2rm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crm2rcrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 3 transformations */
int ipt_ccm2rrrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crrrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_cccrb2rm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crm2ccrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 4 transformations */
int ipt_ccm2rm    (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_crm2cm    (plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);


int ipt_cpanel2all(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_call2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_cpanel2tile(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_ctile2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
#endif /* CGECFI2_H */
