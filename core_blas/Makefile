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

ZSRC=	core_zgeqrt.c core_ztsqrt.c core_ztsmqr.c			\
	core_zgelqt.c core_ztslqt.c core_ztsmlq.c			\
		      core_zttqrt.c core_zttmqr.c			\
		      core_zttlqt.c core_zttmlq.c			\
	core_zgetrf_incpiv.c core_ztstrf.c core_zgessm.c core_zssssm.c  \
	core_zpotrf.c core_ztrtri.c core_zlauum.c                       \
	core_zgemm.c  core_zsymm.c  core_zsyrk.c core_ztrsm.c           \
	core_ztrmm.c  core_zsyr2k.c core_zgeadd.c core_ztradd.c		\
	core_zlacpy.c core_zlange.c core_zlansy.c core_zlantr.c		\
	core_dzasum.c core_ztrasm.c core_zunmqr.c core_zunmlq.c		\
	core_zhemm.c  core_zherk.c  core_zlanhe.c core_zher2k.c         \
	core_zlag2c.c                                                   \
	core_zgessq.c core_zhessq.c core_zsyssq.c core_ztrssq.c         \
	core_zswpab.c core_zgetrip.c core_zshift.c                      \
	core_zplrnt.c core_zplghe.c core_zplgsy.c		        \
	core_zpltmg.c core_zpltmg_fiedler.c core_zpltmg_toeppd.c	\
	core_zpltmg_condex.c core_zpltmg_chebvand.c core_zpltmg_circul.c\
	core_zpltmg_hankel.c			                        \
	core_zlaswp.c core_zgetrf.c                                     \
	core_zgetrf_rectil.c core_zgetrf_reclap.c                       \
	core_zlaset.c core_zlaset2.c                                    \
	core_zherfb.c core_zhegst.c                                     \
	core_ztsmqr_hetra1.c core_ztsmqr_corner.c                       \
	core_ztsmlq_hetra1.c core_ztsmlq_corner.c                       \
	core_zparfb.c core_zpamm.c core_zpemv.c                         \
	core_zlarfy.c core_zlatro.c core_zlarfb_gemm.c                  \
	core_zgemv.c                                                    \
	core_zsetvar.c                                                  \
        \
	core_zgeqp3_init.c core_zgeqp3_larfg.c core_zgeqp3_norms.c	\
	core_zgeqp3_pivot.c core_zgeqp3_update.c	                \
        \
	core_zhbtype1cb.c core_zhbtype2cb.c core_zhbtype3cb.c   \
	core_zgbtype1cb.c core_zgbtype2cb.c core_zgbtype3cb.c   \
        core_zbrdalg1.c core_ztrdalg1.c                         \
	core_zhbelr.c core_zhbrce.c core_zhblrx.c               \
	core_zgbelr.c core_zgbrce.c core_zgblrx.c               \
	core_zlarfx_tbrd.c                                      \
	\
	core_zgetf2_nopiv.c core_zgetrf_nopiv.c                 \
	core_zlacpy_pivot.c core_zgeqp3_tntpiv.c                \
        \
	core_zswap.c core_zlascl.c core_zlascal.c core_dlag2z.c \
	core_zsteqr.c core_zstedc.c      			\
	core_dlapst.c      		                	\
        core_dlaed2_compressq.c                                 \
	core_dlaed2_computeK.c                                  \
	core_dlaed3_computeW.c                                  \
	core_dlaed3_computevectors.c       			\
	core_dlaed3_merge.c                    			\
	core_dlaed3_reduceW.c                                   \
	core_dlaed3_updatevectors.c       			\
	core_dlaed4.c                   			\
	core_dlaed0_betaapprox.c

-include Makefile.src

ALLHDR = primes.h gkkleader.h

ALLSRC = global.c	                \
	core_pivot_update.c             \
        primes.c gkkleader.c            \
        descriptor.c                    \
	$(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBCOREBLAS)

$(LIBCOREBLAS): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
