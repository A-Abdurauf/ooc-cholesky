/**
 *
 * @file eztrace_convert_coreblas.c
 *
 *  PLASMA core_blas tracing kernels
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * This file provides the functions to generate the trace
 * in function of the events.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 * Easy way to add new kernel:
 *   1 - Add the enum coreblas_ev_codes.h: COREBLAS_*kernelname* (Don't forget to check if it already exist or not)
 *   2 - Add the line to initialize your new kernel in the eztrace_convert_coreblas_init function
 *
 **/
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <GTG.h>
#include <ev_codes.h>
#include <eztrace_list.h>
#include <eztrace_convert.h>
#include "coreblas_ev_codes.h"

#if defined( _WIN32 ) || defined( _WIN64 )
#define my_getenv(var, str) {                        \
        int len = 512;                               \
        int ret;                                     \
        str = (char*)malloc(len * sizeof(char));     \
        ret = GetEnvironmentVariable(var, str, len); \
        if (ret == 0) {                              \
            free(str);                               \
            str = NULL;                              \
        }                                            \
    }

#define my_cleanenv(str) if (str != NULL) free(str);

#else /* Other OS systems */

#define my_getenv(var, str)  envstr = getenv(var);
#define my_cleanenv(str)

#endif

int get_convert_format() {
    char *envstr;
    int rc = 0;

    my_getenv("PLASMA_TRACE_SEQUENCE", envstr);
    if ( envstr != NULL ) {
        rc = atoi(envstr);
    }
    my_cleanenv(envstr);
    return rc;
}

#ifndef min
#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef max
#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

#define COREBLAS_STATE      "ST_Thread"

#define COREBLAS_THREADS_MAX 4096


#define COREBLAS_NAME(ev) coreblas_array[ (int)(COREBLAS_GET_CODE(ev)) ].name
#define GTG_RANDOM gtg_get_random_color()

#define COREBLAS_INIT_EVENT( _idx_, _name_, _color_ )     \
    coreblas_array[_idx_].name  = _name_;                 \
    coreblas_array[_idx_].color = _color_;                \
    coreblas_array[_idx_].nb    = 0;                      \


/*
 * Check EZTrace version for compatibility
 */
#if !defined(EZTRACE_API_VERSION) || !(EZTRACE_API_VERSION > 0x00000400)
#error "EZTrace 0.7 or greater is required"
#endif

/*
 * Data for each event
 *   @name: name of the kernel (several kernels can use the same name, for
 *          example: gemm & gemm_f1, or laswp & laswpc
 *   @color: color assigne to the kernel
 *   @nb: number of calls to the kernel (with eztrace_stats)
 *   @sum: Total time spent executing this kernel (with eztrace_stats)
 *   @min: Minimal execution time of this kernel (with eztrace_stats)
 *   @max: Maximal execution time of this kernel (with eztrace_stats)
 */
typedef struct coreblas_s {
    char *name;
    gtg_color_t color;
    int  nb;
    double sum;
    double min;
    double max;
} coreblas_t;

static int         coreblas_array_initialized = 0;
static coreblas_t  coreblas_array[COREBLAS_NBMAX_EVENTS];
static gtg_color_t colors_array[20];

/*
 * Keep information on thread status.
 *   @tid: Thread Id
 *   @active: Thread is active/inactive
 *   @nbtask: Number of tasks pushed to threads when using pop/push
 *   @lasttime: Start time of this task
 */
typedef struct coreblas_thrdstate_s {
    unsigned int tid;
    int  active;
    int  nbtask;
    double lasttime;
} coreblas_thrdstate_t;

static coreblas_thrdstate_t *thrdstate = NULL;
static int nbtrhd = 0;

static inline gtg_color_t gtg_get_random_color(){
    static int i = -1;
    i = (i+1)%20;
    return colors_array[i];
}

static int get_short_tid( unsigned int thread_id )
{
    int tid;

    for (tid=0; tid<nbtrhd; tid++) {
        if ( thrdstate[tid].tid == thread_id )
            break;
    }

    /* Thread not found, we add it */
    if ( tid == nbtrhd ) {
        if ( tid < COREBLAS_THREADS_MAX ) {
            thrdstate[nbtrhd].tid = thread_id;
            nbtrhd++;
        }
        else {
            fprintf(stderr, "Too many threads, increase COREBLAS_THREADS_MAX and recompile\n");
            return -1;
        }
    }
    return tid;
}

/*
 * Case with priority and request need to be handle correctly in GTG
 * before to be included here
 */
#define MAX_SEQUENCE 100
typedef struct sequence_s {
    uint64_t id;
    char *name;
} sequence_t;

sequence_t *seqtab;
static int use_sequence = 0;

void sequenceInit(){
    seqtab = (sequence_t*)malloc((MAX_SEQUENCE+1) * sizeof(sequence_t));
    memset(seqtab, 0, (MAX_SEQUENCE+1) * sizeof(sequence_t));

    seqtab[0].id   = -1;
    seqtab[0].name = "SequenceOutOfRange";
    seqtab++;
}

void sequenceDestroy(){
    int i=0;
    while( i < MAX_SEQUENCE && seqtab[i].id != 0)
    {
        free(seqtab[i].name);
        i++;
    }
    seqtab--;
    free(seqtab);
}

int getSequence(uint64_t seq)
{
    int i=0;

    while ( (i < MAX_SEQUENCE)
            && (seqtab[i].id != 0)
            && (seqtab[i].id != seq) )
        i++;

    if (i < MAX_SEQUENCE)
    {
        if ( seqtab[i].id == seq )
        {
            return i;
        }
        else
        {
            seqtab[i].id = seq;
            if ( asprintf(&(seqtab[i].name), "Sequence%03d", i) < 0 ) {
                fprintf(stderr, "Failed to create new sequence name\n");
                exit(-1);
            }

            addEntityValue(seqtab[i].name, COREBLAS_STATE, seqtab[i].name, colors_array[i%20] );
            return i;
        }
    } else {
        fprintf(stderr, "WARNING: Too many sequences, you need to increase the limit and recompile\n");
        return -1;
    }
}

void handle_coreblas_start(eztrace_event_t *ev)
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    if ( use_sequence && GET_NBPARAMS(ev) > 0 ) {
        char *seqname = seqtab[ getSequence(GET_PARAM(ev, 1)) ].name;
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, seqname );
    } else {
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, COREBLAS_NAME(ev) );
    }
}

void
handle_coreblas_stop ()
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, "wait");
}

int
eztrace_convert_coreblas_init()
{
    int i;
    use_sequence = get_convert_format();

    if ( coreblas_array_initialized == 0 ) {

        /* Initialize the colors_array */
        colors_array[ 0] = GTG_RED;
        colors_array[ 1] = GTG_GREEN;
        colors_array[ 2] = GTG_BLUE;
        colors_array[ 3] = GTG_WHITE;
        colors_array[ 4] = GTG_TEAL;
        colors_array[ 5] = GTG_DARKGREY;
        colors_array[ 6] = GTG_YELLOW;
        colors_array[ 7] = GTG_PURPLE;
        colors_array[ 8] = GTG_LIGHTBROWN;
        colors_array[ 9] = GTG_DARKBLUE;
        colors_array[10] = GTG_PINK;
        colors_array[11] = GTG_DARKPINK;
        colors_array[12] = GTG_SEABLUE;
        colors_array[13] = GTG_KAKI;
        colors_array[14] = GTG_REDBLOOD;
        colors_array[15] = GTG_BROWN;
        colors_array[16] = GTG_GRENAT;
        colors_array[17] = GTG_ORANGE;
        colors_array[18] = GTG_MAUVE;
        colors_array[19] = GTG_LIGHTPINK;

        /* First initialization to fill in the gap */
        for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
            coreblas_array[i].name  = "";
            coreblas_array[i].color = GTG_RANDOM;
            coreblas_array[i].nb    = -1;
            coreblas_array[i].sum   = 0.;
            coreblas_array[i].min   = 999999999999.;
            coreblas_array[i].max   = 0.;
        }

        /*
         * Some kernels have a predefined color to keep them from one figure to another.
         * Those which has never been assigned a color, use a random one
         */
        COREBLAS_INIT_EVENT(COREBLAS_FREE,  "free",   GTG_DARKGREY  );
        COREBLAS_INIT_EVENT(COREBLAS_FOO,   "foo",    GTG_DARKGREY  );
        COREBLAS_INIT_EVENT(COREBLAS_GEMM,  "gemm",   GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_HERK,  "herk",   GTG_WHITE     );
        COREBLAS_INIT_EVENT(COREBLAS_SYRK,  "syrk",   GTG_WHITE     );
        COREBLAS_INIT_EVENT(COREBLAS_HEMM,  "hemm",   GTG_DARKPINK  );
        COREBLAS_INIT_EVENT(COREBLAS_SYMM,  "symm",   GTG_DARKPINK  );
        COREBLAS_INIT_EVENT(COREBLAS_TRMM,  "trmm",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_TRSM,  "trsm",   GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_HER2K, "her2k",  GTG_PINK      );
        COREBLAS_INIT_EVENT(COREBLAS_SYR2K, "syr2k",  GTG_PINK      );
        COREBLAS_INIT_EVENT(COREBLAS_GEMV,  "gemv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_GBMV,  "gbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HEMV,  "hemv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HBMV,  "hbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HPMV,  "hpmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SYMV,  "symv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SBMV,  "sbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SPMV,  "spmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TRMV,  "trmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TBMV,  "tbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TPMV,  "tpmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TRSV,  "trsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TBSV,  "tbsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TPSV,  "tpsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_GER,   "ger",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_GERU,  "geru",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_GERC,  "gerc",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HER,   "her",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HPR,   "hpr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HER2,  "her2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HPR2,  "hpr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SYR,   "syr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SPR,   "spr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SYR2,  "syr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SPR2,  "spr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_ROTG,  "rotg",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROTMG, "rotmg",  GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROT,   "rot",    GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROTM,  "rotm",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_SWAP,  "swap",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_SCAL,  "scal",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_COPY,  "copy",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_AXPY,  "axpy",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_GEADD, "geadd",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_DOT,   "dot",    GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_DOTU,  "dotu",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_DOTC,  "dotc",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_XDOT,  "xdot",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_NRM2,  "nrm2",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_ASUM,  "asum",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_AMAX,  "amax",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LACPY, "lacpy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANGE, "lange",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANHE, "lanhe",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANSY, "lansy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFB, "larfb",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFT, "larft",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_LASWP, "laswp",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_LAUUM, "lauum",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_POTRF, "potrf",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_TRTRI, "trtri",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LASET, "laset",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LASSQ, "lassq",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_GELQT, "gelqt",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_GEQRT, "geqrt",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_GESSM, "gessm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_GETRF, "getrf",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_LATRO, "latro",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_SSSSM, "ssssm",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TITRO, "titro",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_TRBMM, "trbmm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TRGMM, "trgmm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TSLQT, "tslqt",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_TSMLQ, "tsmlq",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TSMQR, "tsmqr",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TSQRT, "tsqrt",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_TSRFB, "tsrfb",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TSTRF, "tstrf",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TTLQT, "ttlqt",  GTG_REDBLOOD  );
        COREBLAS_INIT_EVENT(COREBLAS_TTMLQ, "ttmlq",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TTMQR, "ttmqr",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TTQRT, "ttqrt",  GTG_REDBLOOD  );
        COREBLAS_INIT_EVENT(COREBLAS_TTRFB, "ttrfb",  GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_UNMLQ, "unmlq",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_UNMQR, "unmqr",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_GETRIP,"getrip", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLGHE, "plghe",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLGSY, "plgsy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SHIFT, "shift",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SHIFTW,"shiftw", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SWPAB, "swpab",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLRNT, "plrnt",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PEMV,  "pemv",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_BRDALG,"brdalg", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_TRDALG,"trdalg", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_HEGST, "hegst",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SYGST, "sygst",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_HERFB, "herfb",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SYRFB, "syrfb",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFG,        "larfg",      GTG_BLUE     );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_INIT,   "qp3_init",   GTG_BLUE     );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_NORMS,  "qp3_norms",  GTG_YELLOW   );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_PIVOT,  "qp3_pivot",  GTG_REDBLOOD );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_UPDATE, "qp3_update", GTG_GREEN    );
        COREBLAS_INIT_EVENT(COREBLAS_PIVOT_UPDATE, "pivot_update", GTG_GREEN    );
        COREBLAS_INIT_EVENT(COREBLAS_SETVAR, "setvar", GTG_ORANGE   );
        COREBLAS_INIT_EVENT(COREBLAS_STEDC,  "stedc",  GTG_YELLOW   );
        COREBLAS_INIT_EVENT(COREBLAS_STEQR,  "steqr",  GTG_YELLOW   );
        COREBLAS_INIT_EVENT(COREBLAS_LASCL,  "lascl",  GTG_YELLOW   );

        COREBLAS_INIT_EVENT(COREBLAS_LAED2_COMPRESSQ,     "COPY_Q_TO_Q2",    GTG_SEABLUE );
        COREBLAS_INIT_EVENT(COREBLAS_LAED2_COMPUTEK,      "COMPUTE_K",       GTG_RED     );
        COREBLAS_INIT_EVENT(COREBLAS_LAED2_COPYDEF,       "COPY_Q2_TO_Q",    GTG_WHITE   );
        COREBLAS_INIT_EVENT(COREBLAS_LAED3_COMPVEC,       "COMPUTE_VECTORS", GTG_PURPLE  );
        COREBLAS_INIT_EVENT(COREBLAS_LAED3_WSCOPY,        "LAED3_WSCOPY",    GTG_PURPLE  );
        COREBLAS_INIT_EVENT(COREBLAS_LAED3_COMPW,         "COMPUTE_W",       GTG_BLUE    );
        COREBLAS_INIT_EVENT(COREBLAS_LAED3_REDUCEW,       "REDUCE W",        GTG_ORANGE  );
        COREBLAS_INIT_EVENT(COREBLAS_LAED3_UPDATEVECTORS, "UPDATE_VECTORS",  GTG_GREEN   );
        COREBLAS_INIT_EVENT(COREBLAS_LAED4,               "LAED4",           GTG_BLUE    );
        COREBLAS_INIT_EVENT(COREBLAS_LAG2Z,               "LAG2Z",           GTG_TEAL    );

        COREBLAS_INIT_EVENT(COREBLAS_TYPE1, "TYPE1", GTG_TEAL );
        COREBLAS_INIT_EVENT(COREBLAS_TYPE2, "TYPE2", GTG_TEAL );
        COREBLAS_INIT_EVENT(COREBLAS_TYPE3, "TYPE3", GTG_TEAL );

        coreblas_array_initialized = 1;
    }

    /*
     * Init only for trace conversion, not for statistics
     */
    if(get_mode() == EZTRACE_CONVERT) {
        for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
            if ( coreblas_array[i].nb == -1 )
                continue;

            addEntityValue( coreblas_array[i].name,
                            COREBLAS_STATE,
                            coreblas_array[i].name,
                            coreblas_array[i].color );
        }

        /* plasma coreblas */
        addEntityValue ("wait", COREBLAS_STATE, "wait", GTG_BLACK );

        if (use_sequence)
            sequenceInit();
    }
    return 0;
}

void
eztrace_convert_coreblas_finalize()
{

}


int
handle_coreblas_events(eztrace_event_t *ev)
{
    double time;

    if ( thrdstate == NULL ) {
        thrdstate = (coreblas_thrdstate_t*)malloc(COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
        memset( thrdstate, 0, COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
    }

    if ( IS_A_COREBLAS_EV(ev) ) {
        switch (LITL_READ_GET_CODE(ev)) {
        case FUT_COREBLAS(STOP)  :
        {
            int tid = get_short_tid( CUR_THREAD_ID );
            assert(tid != -1);

            thrdstate[tid].nbtask--;

            if ( thrdstate[tid].nbtask < 0 ) {
                fprintf(stderr, "WARNING: The end of a state appears before the beginning\n");
                thrdstate[tid].nbtask = 0;
                return 0;
            }

            if ( thrdstate[tid].nbtask ==  0 ) {
                time = ( CURRENT - thrdstate[tid].lasttime );

                /* Check that we have an existing state */
                assert(  coreblas_array[ thrdstate[tid].active ].nb >= 0 );

                if( coreblas_array[ thrdstate[tid].active ].nb == 0 ) {
                    coreblas_array[ thrdstate[tid].active ].sum = 0.;
                    coreblas_array[ thrdstate[tid].active ].max = 0.;
                    coreblas_array[ thrdstate[tid].active ].min = 999999999999.;
                }
                coreblas_array[ thrdstate[tid].active ].nb++;
                coreblas_array[ thrdstate[tid].active ].sum += time;
                coreblas_array[ thrdstate[tid].active ].max = max( coreblas_array[ thrdstate[tid].active ].max, time );
                coreblas_array[ thrdstate[tid].active ].min = min( coreblas_array[ thrdstate[tid].active ].min, time );

                thrdstate[tid].active = 0;
                thrdstate[tid].lasttime = 0;

                if( get_mode() == EZTRACE_CONVERT ) {
                    handle_coreblas_stop(ev);
                }
            }
            return 1;
        }
        break;

        default: /* All the different states */
        {
            int ev_code = COREBLAS_GET_CODE(ev);
            int tid = get_short_tid( CUR_THREAD_ID );

            assert(tid != -1);
            thrdstate[tid].nbtask++;

            if ( thrdstate[tid].active != 0 ) {
                fprintf(stderr,
                        "WARNING: Thread %d change to state %s (%d) before to stop previous state %s (%d)\n"
                        "         This might happen if tracing functions\n",
                        (int)CUR_THREAD_ID,
                        coreblas_array[ ev_code ].name, ev_code,
                        coreblas_array[ thrdstate[tid].active ].name, thrdstate[tid].active );
            }
            else {
                thrdstate[tid].active = ev_code;
                thrdstate[tid].lasttime = CURRENT;

                if( get_mode() == EZTRACE_CONVERT ) {
                    handle_coreblas_start(ev);
                }
            }
            return 1;
        }
        }
    }
    else
        return 0;
}

/*
 * Print the results of statistics.
 */
void print_coreblas_stats() {
    int i;

    printf ( "\nCoreblas Module:\n");
    printf (   "-----------\n");

    for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
        if ( coreblas_array[ i ].nb > 0 ) {
            printf ( "%s : %d calls\n"
                     "\tAverage time: %.3f ms\n"
                     "\tMaximun time: %.3f ms\n"
                     "\tMinimun time: %.3f ms\n",
                     coreblas_array[ i ].name,
                     coreblas_array[ i ].nb,
                     coreblas_array[ i ].sum / (double)(coreblas_array[ i ].nb),
                     coreblas_array[ i ].max,
                     coreblas_array[ i ].min);
        }
    }
}

struct eztrace_convert_module coreblas_module;

void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
  coreblas_module.api_version = EZTRACE_API_VERSION;

  /* Specify the initialization function.
   * This function will be called once all the plugins are loaded
   * and the trace is started.
   * This function usually declared StateTypes, LinkTypes, etc.
   */
  coreblas_module.init = eztrace_convert_coreblas_init;

  /* Specify the function to call for handling an event
   */
  coreblas_module.handle = handle_coreblas_events;

  /* Specify the function to call for handling an event when eztrace_stats is called
   */
  coreblas_module.handle_stats = handle_coreblas_events;

  /* Print the results of statistics
   */
  coreblas_module.print_stats = print_coreblas_stats;

  /* Specify the module prefix */
  coreblas_module.module_prefix = COREBLAS_EVENTS_ID;

  if ( asprintf(&coreblas_module.name, "coreblas") < 0 ) {
      fprintf(stderr, "Failed to create module name\n");
      exit(-1);
  }
  if ( asprintf(&coreblas_module.description, "Module for kernels used in PLASMA (BLAS, LAPACK and coreblas)") < 0 ) {
      fprintf(stderr, "Failed to create module description\n");
      exit(-1);
  }

  coreblas_module.token.data = &coreblas_module;

  /* Register the module to eztrace_convert */
  eztrace_convert_register_module(&coreblas_module);
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
    if(use_sequence && (get_mode() == EZTRACE_CONVERT)) {
        sequenceDestroy();
    }
}
