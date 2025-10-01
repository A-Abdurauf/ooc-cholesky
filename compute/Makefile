###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.8.0
# @author Jakub Kurzak
# @date 2010-11-15
#
###

# Overwritten in make.inc
PLASMA_DIR = ..
include ../Makefile.internal

INC := -I$(PLASMA_DIR)/control/ $(INC)

ZHDR =  zgecfi2.h

ZSRC =  \
	zgels.c  zgelqf.c zgelqs.c zgeqrf.c zgeqrs.c		\
	zgesv_incpiv.c  zgetrf_incpiv.c zgetrs_incpiv.c		\
	zgesv.c  zgetrf.c zgetrs.c zlaswp.c zlaswpc.c		\
	zposv.c  zpotrf.c zpotrs.c				\
	ztrsm.c  zgemm.c  zsymm.c zsyrk.c ztrsmpl.c ztrsmrv.c	\
        zgeadd.c ztradd.c                                       \
	zpotri.c ztrtri.c zlauum.c ztrmm.c zsyr2k.c		\
	zgetri.c						\
	zlange.c zlansy.c zlantr.c zlaset.c zlacpy.c zlascal.c	\
	zunmqr.c zunmlq.c zungqr.c zunglq.c			\
	zhemm.c  zherk.c  zlanhe.c zher2k.c			\
	zcposv.c zcgesv.c zcungesv.c				\
	zgecfi2.c zgecfi.c zgetmi.c				\
	zplrnt.c zpltmg.c zplghe.c zplgsy.c				\
	zhegv.c zhegvd.c  zheev.c zheevd.c zheevr.c zhegst.c zhetrd.c	\
	zgesvd.c zgesdd.c zgebrd.c 				\
	zgecon.c zpocon.c					\
	zgeqp3.c						\
	zgetrf_nopiv.c zgetrf_tntpiv.c                          \
	zstedc.c					        \
	\
	pzgelqf.c pzgeqrf.c pzpotrf.c pztrtri.c pzlauum.c		\
	pzgelqfrh.c pzgeqrfrh.c						\
	pztrsm.c  pzgemm.c  pzsymm.c  pzsyrk.c pzsyr2k.c                \
	pztrsmpl.c pztrsmrv.c	                                        \
	pzlange.c pzlansy.c pzlantr.c				        \
	pzlacpy.c pztradd.c pzgetrf_incpiv.c pztrmm.c			\
	pzunmqr.c pzunmlq.c pzungqr.c pzunglq.c				\
	pzunmqrrh.c pzungqrrh.c						\
	pzunmlqrh.c pzunglqrh.c						\
	pzlanhe.c pzhemm.c pzherk.c pzher2k.c				\
	pzlag2c.c pzlaset.c pzlaset2.c pzlascal.c			\
	pzpack.c  pzshift.c pzgetmi2.c					\
	pzpltmg.c pzplghe.c pzplgsy.c					\
	pzpltmg_condex.c pzpltmg_fiedler.c pzpltmg_toeppd.c		\
	pzpltmg_circul.c pzpltmg_chebvand.c pzpltmg_hankel.c		\
	pzpltmg_house.c                                                 \
	pzgetrf_reclap.c pzgetrf_rectil.c pzlaswp.c pzlaswpc.c		\
	pzbarrier.c							\
	pzgebrd_ge2gb.c pzgebrd_gb2bd_v1.c				\
	pzhetrd_he2hb.c pzhetrd_hb2st_v1.c				\
	pzhegst.c							\
	pzlarft_blgtrd.c pzunmqr_blgtrd.c				\
	pzhbcpy_t2bl.c pzgbcpy_t2bl.c					\
	pzgeqp3.c               					\
	pzgetrf_nopiv.c pzgetrf_tntpiv.c                                \
	pzswaps.c pzlaset_identity.c pzstedc.c                          \
	\
	pdlaed0.c pdlaed1.c pdlag2z.c

-include Makefile.src

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBPLASMA)

$(LIBPLASMA): $(HEADERS) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%.o : %.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

zgecfi2.c zgecfi.c : $(ZHDR)
cgecfi2.c cgecfi.c : $(CHDR)
dgecfi2.c dgecfi.c : $(DHDR)
sgecfi2.c sgecfi.c : $(SHDR)

.PHONY: clean cleanall
