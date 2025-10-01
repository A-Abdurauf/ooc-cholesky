###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.8.0
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

# Add the headers from src directory
INC := -I../control/ $(INC)

ZSRC=	qwrapper_zgeqrt.c qwrapper_ztsqrt.c qwrapper_ztsmqr.c			        \
	qwrapper_zgelqt.c qwrapper_ztslqt.c qwrapper_ztsmlq.c			        \
		          qwrapper_zttqrt.c qwrapper_zttmqr.c			        \
	                  qwrapper_zttlqt.c qwrapper_zttmlq.c			        \
	qwrapper_zgetrf_incpiv.c qwrapper_ztstrf.c qwrapper_zgessm.c qwrapper_zssssm.c  \
	qwrapper_zpotrf.c qwrapper_ztrtri.c qwrapper_zlauum.c                           \
	qwrapper_zgemm.c  qwrapper_zsymm.c  qwrapper_zsyrk.c qwrapper_ztrsm.c           \
	qwrapper_ztrmm.c  qwrapper_zsyr2k.c qwrapper_zgeadd.c qwrapper_ztradd.c	        \
	qwrapper_zlacpy.c qwrapper_zlange.c qwrapper_zlansy.c qwrapper_zlantr.c	        \
	qwrapper_dzasum.c qwrapper_ztrasm.c qwrapper_zunmqr.c qwrapper_zunmlq.c	        \
	qwrapper_zhemm.c  qwrapper_zherk.c  qwrapper_zlanhe.c qwrapper_zher2k.c         \
	qwrapper_zlag2c.c                                                               \
	qwrapper_zplssq.c qwrapper_zgessq.c qwrapper_zhessq.c qwrapper_zsyssq.c qwrapper_ztrssq.c        \
	qwrapper_zswpab.c qwrapper_zgetrip.c qwrapper_zshift.c                          \
	qwrapper_zplrnt.c qwrapper_zplghe.c qwrapper_zplgsy.c			        \
	qwrapper_zpltmg.c qwrapper_zpltmg_fiedler.c qwrapper_zpltmg_toeppd.c		\
	qwrapper_zpltmg_chebvand.c qwrapper_zpltmg_circul.c qwrapper_zpltmg_hankel.c    \
	qwrapper_zlaswp.c qwrapper_zgetrf.c                                             \
	qwrapper_zgetrf_rectil.c qwrapper_zgetrf_reclap.c                               \
	qwrapper_zlaset.c qwrapper_zlaset2.c                                            \
	qwrapper_zherfb.c qwrapper_zhegst.c                                             \
	qwrapper_ztsmqr_hetra1.c qwrapper_ztsmqr_corner.c                               \
	qwrapper_ztsmlq_hetra1.c qwrapper_ztsmlq_corner.c                               \
	qwrapper_zpamm.c                                                                \
        qwrapper_zlatro.c                                                               \
	qwrapper_zgemv.c qwrapper_zgemv_tile.c qwrapper_zgemm_tile.c                    \
	qwrapper_zsetvar.c                                                              \
	\
	qwrapper_zswap.c  		                                                \
	qwrapper_dlag2z.c  		                                                \
	qwrapper_zlascl.c qwrapper_zlascal.c                                            \
	qwrapper_zsteqr.c                                                               \
	qwrapper_zstedc.c                                                               \
	\
	qwrapper_dlaed0_betaapprox.c  			                                \
	qwrapper_dlaed1_pipelined.c  		                                        \
	qwrapper_dlaed2_compressq.c                                                     \
	qwrapper_dlaed2_computeK.c                                                      \
	qwrapper_dlaed2_copyvectorsback.c   qwrapper_dDC_fakedep.c 	                \
	qwrapper_dlaed3_computeW.c                                                      \
	qwrapper_dlaed3_computevectors.c                                                \
	qwrapper_dlaed3_freebigwork.c  	                                                \
	qwrapper_dlaed3_pipelined.c                                                     \
	qwrapper_dlaed3_reduceW.c                                                       \
	qwrapper_dlaed3_updatevectors.c                                                 \
	qwrapper_dlaed4.c                                                               \
	qwrapper_zlaset_identity.c  			                                \
        \
	qwrapper_zgeqp3_init.c qwrapper_zgeqp3_larfg.c qwrapper_zgeqp3_norms.c	        \
	qwrapper_zgeqp3_pivot.c qwrapper_zgeqp3_update.c	                        \
        \
        qwrapper_zbrdalg1.c qwrapper_ztrdalg1.c                                         \
	\
	qwrapper_zgetrf_nopiv.c                                                         \
	qwrapper_zlacpy_pivot.c qwrapper_zgeqp3_tntpiv.c

-include Makefile.src

ALLSRC = qwrapper_alloc.c       \
        qwrapper_pivot_update.c \
	$(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBCOREBLASQW)

$(LIBCOREBLASQW): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
