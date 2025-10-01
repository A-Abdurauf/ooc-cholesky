/**
 *
 * @file async.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_ASYNC_H_
#define _PLASMA_ASYNC_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
int  plasma_request_fail(PLASMA_sequence *sequence, PLASMA_request *request, int error);
int  plasma_sequence_create(plasma_context_t *plasma, PLASMA_sequence **sequence);
int  plasma_sequence_destroy(plasma_context_t *plasma, PLASMA_sequence *sequence);
int  plasma_sequence_wait(plasma_context_t *plasma, PLASMA_sequence *sequence);
void plasma_sequence_flush(Quark *quark, PLASMA_sequence *sequence, PLASMA_request *request, int status);

#ifdef __cplusplus
}
#endif

#endif
