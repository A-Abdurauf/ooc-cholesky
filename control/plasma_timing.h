/**
 *
 * @file plasma_timing.h
 *
 *  PLASMA internal timing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @date 2011-05-15
 *
 **/
#ifndef _PLASMA_TIMING_H_
#define _PLASMA_TIMING_H_

typedef double plasma_time_t;

#if defined(PLASMA_ENABLE_TIMER)

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#include <time.h>
#include <sys/timeb.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME         ft;
    unsigned __int64 tmpres = 0;
    static int       tzflag;

    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);
        tmpres |=  ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |=  ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres /= 10;  /*convert into microseconds*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;

        tv->tv_sec  = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }
    if (NULL != tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime     = _daylight;
    }
    return 0;
}

#else  /* Non-Windows */
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

/**
 *****************************************************************************
 *
 *  plasma_gettime - Routines to compute the time of internal subroutines. This
 *  function is defined only when the flag PLASMA_ENABLE_TIMER is enabled.
 *
 *******************************************************************************
 *
 * @return
 *          \retval The time in seconds
 *
 *****************************************************************************/
static inline double plasma_gettime(plasma_context_t *plasma)
{
    struct timeval tp;

    plasma_dynamic_sync();
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

#define plasma_printtime(fmt, ...) do { fprintf(stderr, fmt, ##__VA_ARGS__); } while(0)

#else

/**
 *****************************************************************************
 *
 *  plasma_gettime - Routines to compute the time of internal subroutines. This
 *  function is defined only when the flag PLASMA_ENABLE_TIMER is enabled.
 *
 *******************************************************************************
 *
 * @return
 *          \retval The time in seconds
 *
 *****************************************************************************/
static inline double plasma_gettime(plasma_context_t *plasma)
{
    return 0.;
}

#define plasma_printtime(...) do {  } while(0)

#endif /* defined(PLASMA_ENABLE_TIMER) */

#endif /* _PLASMA_TIMING_H_ */
