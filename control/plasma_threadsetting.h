/**
 *
 * @file plasma_threadsetting.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_THREADSETTING_H_
#define _PLASMA_THREADSETTING_H_

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************//**
 *  Internal routines
 **/
void plasma_setlapack_multithreads(int numthreads);
void plasma_setlapack_sequential(plasma_context_t *plasma);
void plasma_psetlapack_numthreads(plasma_context_t *plasma);
/***************************************************************************/
#ifdef __cplusplus
}
#endif

#endif
