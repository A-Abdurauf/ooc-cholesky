/**
 *
 * @file core_blas_dag.h
 *
 *  PLASMA macros to create dag object when using Quark scheduler
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CORE_BLAS_DAG_H_
#define _PLASMA_CORE_BLAS_DAG_H_

/**
 * Macros to set the DAG properties
 */
#if defined(QUARK_DOT_DAG_ENABLE) /* || 1 */
#define DAG_SET_PROPERTIES( _name, _color )                            \
  QUARK_Task_Flag_Set(task_flags, TASK_LABEL, (intptr_t)(_name));      \
  QUARK_Task_Flag_Set(task_flags, TASK_COLOR, (intptr_t)(_color));
#else
#define DAG_SET_PROPERTIES( _name, _color )
#endif

#define DAG_CORE_ASUM   DAG_SET_PROPERTIES( "ASUM"  , "white"   )
#define DAG_CORE_AXPY   DAG_SET_PROPERTIES( "AXPY"  , "white"   )
#define DAG_CORE_GEADD  DAG_SET_PROPERTIES( "GEADD" , "white"   )
#define DAG_CORE_GELQT  DAG_SET_PROPERTIES( "GELQT" , "green"   )
#define DAG_CORE_GEMM   DAG_SET_PROPERTIES( "GEMM"  , "yellow"  )
#define DAG_CORE_GEMV   DAG_SET_PROPERTIES( "GEMV"  , "blue"    )
#define DAG_CORE_GEQRT  DAG_SET_PROPERTIES( "GEQRT" , "green"   )
#define DAG_CORE_GESSM  DAG_SET_PROPERTIES( "GESSM" , "cyan"    )
#define DAG_CORE_GETRF  DAG_SET_PROPERTIES( "GETRF" , "green"   )
#define DAG_CORE_GETRIP DAG_SET_PROPERTIES( "GETRIP", "white"   )
#define DAG_CORE_HEMM   DAG_SET_PROPERTIES( "HEMM"  , "white"   )
#define DAG_CORE_HER2K  DAG_SET_PROPERTIES( "HER2K" , "white"   )
#define DAG_CORE_HERFB  DAG_SET_PROPERTIES( "HERFB" , "pink"    )
#define DAG_CORE_HERK   DAG_SET_PROPERTIES( "HERK"  , "yellow"  )
#define DAG_CORE_LACPY  DAG_SET_PROPERTIES( "LACPY" , "white"   )
#define DAG_CORE_LAG2C  DAG_SET_PROPERTIES( "LAG2C" , "white"   )
#define DAG_CORE_LAG2Z  DAG_SET_PROPERTIES( "LAG2Z" , "white"   )
#define DAG_CORE_LANGE  DAG_SET_PROPERTIES( "LANGE" , "white"   )
#define DAG_CORE_LANHE  DAG_SET_PROPERTIES( "LANHE" , "white"   )
#define DAG_CORE_LANSY  DAG_SET_PROPERTIES( "LANSY" , "white"   )
#define DAG_CORE_LASCL  DAG_SET_PROPERTIES( "LASCL" , "white"   )
#define DAG_CORE_LASET  DAG_SET_PROPERTIES( "LASET" , "orange"  )
#define DAG_CORE_LASSQ  DAG_SET_PROPERTIES( "LASSQ" , "white"   )
#define DAG_CORE_LASWP  DAG_SET_PROPERTIES( "LASWP" , "orange"  )
#define DAG_CORE_LATRO  DAG_SET_PROPERTIES( "LATRO" , "white"   )
#define DAG_CORE_LAUUM  DAG_SET_PROPERTIES( "LAUUM" , "white"   )
#define DAG_CORE_PLGHE  DAG_SET_PROPERTIES( "PLGHE" , "white"   )
#define DAG_CORE_PLGSY  DAG_SET_PROPERTIES( "PLGSY" , "white"   )
#define DAG_CORE_PLRNT  DAG_SET_PROPERTIES( "PLRNT" , "white"   )
#define DAG_CORE_POTRF  DAG_SET_PROPERTIES( "POTRF" , "green"   )
#define DAG_CORE_SHIFT  DAG_SET_PROPERTIES( "SHIFT" , "white"   )
#define DAG_CORE_SHIFTW DAG_SET_PROPERTIES( "SHIFTW", "white"   )
#define DAG_CORE_SSSSM  DAG_SET_PROPERTIES( "SSSSM" , "yellow"  )
#define DAG_CORE_STEDC  DAG_SET_PROPERTIES( "STEDC" , "yellow"  )
#define DAG_CORE_STEQR  DAG_SET_PROPERTIES( "STEQR" , "yellow"  )
#define DAG_CORE_SWPAB  DAG_SET_PROPERTIES( "SWPAB" , "white"   )
#define DAG_CORE_SYMM   DAG_SET_PROPERTIES( "SYMM"  , "white"   )
#define DAG_CORE_SYR2K  DAG_SET_PROPERTIES( "SYR2K" , "white"   )
#define DAG_CORE_SYRK   DAG_SET_PROPERTIES( "SYRK"  , "red"     )
#define DAG_CORE_TRMM   DAG_SET_PROPERTIES( "TRMM"  , "cyan"    )
#define DAG_CORE_TRSM   DAG_SET_PROPERTIES( "TRSM"  , "cyan"    )
#define DAG_CORE_TRTRI  DAG_SET_PROPERTIES( "TRTRI" , "white"   )
#define DAG_CORE_TSLQT  DAG_SET_PROPERTIES( "TSLQT" , "red"     )
#define DAG_CORE_TSMLQ  DAG_SET_PROPERTIES( "TSMLQ" , "yellow"  )
#define DAG_CORE_TSMQR  DAG_SET_PROPERTIES( "TSMQR" , "yellow"  )
#define DAG_CORE_TSQRT  DAG_SET_PROPERTIES( "TSQRT" , "red"     )
#define DAG_CORE_TSTRF  DAG_SET_PROPERTIES( "TSTRF" , "red"     )
#define DAG_CORE_TTLQT  DAG_SET_PROPERTIES( "TTLQT" , "pink"    )
#define DAG_CORE_TTMLQ  DAG_SET_PROPERTIES( "TTMLQ" , "magenta" )
#define DAG_CORE_TTMQR  DAG_SET_PROPERTIES( "TTMQR" , "magenta" )
#define DAG_CORE_TTQRT  DAG_SET_PROPERTIES( "TTQRT" , "pink"    )
#define DAG_CORE_UNMLQ  DAG_SET_PROPERTIES( "UNMLQ" , "cyan"    )
#define DAG_CORE_UNMQR  DAG_SET_PROPERTIES( "UNMQR" , "cyan"    )


#define DAG_CORE_LAED0_BETAAPPROX      DAG_SET_PROPERTIES( "Rank-1 Approx."    , "white"  )
#define DAG_CORE_LAED1_PIPELINED       DAG_SET_PROPERTIES( "COPY_Q_TO_Q2+\nLAED4+COMP_W", "blue"   )
#define DAG_CORE_LAED2_COMPUTEK        DAG_SET_PROPERTIES( "COMPUTE_K"         , "red"    )
#define DAG_CORE_LAED2_COMPRESSQ       DAG_SET_PROPERTIES( "COPY_Q_TO_Q2"      , "deepskyblue" )
#define DAG_CORE_LAED2_COPYDEF         DAG_SET_PROPERTIES( "COPY_Q2_TO_Q"      , "white"   )
#define DAG_CORE_LAED3_COMPVEC         DAG_SET_PROPERTIES( "COMPUTE_VECTORS"   , "purple" )
#define DAG_CORE_LAED3_COMPW           DAG_SET_PROPERTIES( "COMPUTE_W"         , "cyan"   )
#define DAG_CORE_LAED3_FREE            DAG_SET_PROPERTIES( "LAED3_FREE"        , "white"  )
#define DAG_CORE_LAED3_PIPELINED       DAG_SET_PROPERTIES( "COMPUTE+UPDATE VEC", "green"  )
#define DAG_CORE_LAED3_REDUCEW         DAG_SET_PROPERTIES( "REDUCE_W"          , "orange" )
#define DAG_CORE_LAED3_UPDATEVECTORS   DAG_SET_PROPERTIES( "UPDATE_VECTORS"    , "green"  )
#define DAG_CORE_LAED4                 DAG_SET_PROPERTIES( "LAED4"             , "blue"   )
#define DAG_CORE_FAKE_DC               DAG_SET_PROPERTIES( "COPY_SYNC"         , "white"  )

/**
 * Macros to set the trace properties
 */
#if defined(PLASMA_EZTRACE)
#include <eztrace.h>
#include <ev_codes.h>
#include "../core_blas-eztrace/coreblas_ev_codes.h"

#define plasma_trace_setflags( _flags, _NAME )                          \
    QUARK_Task_Flag_Set( _flags, TASK_START_CODE, (intptr_t)(FUT_COREBLAS(_NAME))); \
    QUARK_Task_Flag_Set( _flags, TASK_STOP_CODE,  (intptr_t)(FUT_COREBLAS(STOP )));

#else /* defined(PLASMA_EZTRACE) */

#define plasma_trace_setflags( _flags, _NAME ) do {} while(0)

#endif  /* defined(PLASMA_EZTRACE) */


/**
 * Macros to decide if we trace per kernel, or per function
 */
#if !defined(TRACE_BY_KERNEL) && !defined(TRACE_BY_FUNCTION)
#define TRACE_BY_KERNEL
#endif

#if defined(TRACE_BY_KERNEL)

#define plasma_trace_by_kernel( _flags, _NAME )  plasma_trace_setflags( _flags, _NAME )
#define plasma_gendag_by_kernel( _flags, _NAME ) DAG_CORE_ ## _NAME

#else

#define plasma_trace_by_kernel( _flags, _NAME )  do {} while(0)
#define plasma_gendag_by_kernel( _flags, _NAME ) do {} while(0)

#endif /* defined(TRACE_BY_KERNEL) */

#if defined(TRACE_BY_FUNCTION)

#define plasma_trace_by_function( _flags, _NAME )  plasma_trace_setflags( _flags, _NAME )
#define plasma_gendag_by_function( _flags, _NAME ) DAG_CORE_ ## _NAME

#else

#define plasma_trace_by_function( _flags, _NAME )  do {} while(0)
#define plasma_gendag_by_function( _flags, _NAME ) do {} while(0)

#endif /* defined(TRACE_BY_FUNCTION) */


#define plasma_profile_by_kernel( _flags, _NAME )       \
    {                                                   \
        plasma_gendag_by_kernel( _flags, _NAME );       \
        plasma_trace_by_kernel( _flags, _NAME );        \
    }

#define plasma_profile_by_function( _flags, _NAME )     \
    {                                                   \
        plasma_gendag_by_function( _flags, _NAME );     \
        plasma_trace_by_function( _flags, _NAME );      \
    }

#endif /* _PLASMA_CORE_BLAS_DAG_H_ */
