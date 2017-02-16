# @Author: Jenkins Alec <alec>
# @Date:   2012-07-18T05:54:06-07:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-02-13T18:59:38-08:00

CC=g++

WARNINGS=-W -Wall -Wstrict-prototypes -Wmissing-prototypes -Waggregate-return \
    -Wcast-align -Wcast-qual -Wnested-externs -Wshadow -Wbad-function-cast \
    -Wwrite-strings

SRCFOLDER = ~/UCSB/cpp_files/kiss_fft130

CFLAGS=-O3 -I $(SRCFOLDER) -I $(SRCFOLDER)/tools $(WARNINGS)
CFLAGS+=-ffast-math -fomit-frame-pointer
#CFLAGS+=-funroll-loops
#CFLAGS+=-march=prescott
#CFLAGS+= -mtune=native
# TIP: try adding -openmp or -fopenmp  to enable OPENMP directives and use of multiple cores
#CFLAGS+=-fopenmp
CFLAGS+= $(CFLAGADD)

ifeq "$(NFFT)" ""
 NFFT=1800
endif
ifeq "$(NUMFFTS)" ""
 NUMFFTS=10000
endif

ifeq "$(DATATYPE)" ""
 DATATYPE=double
endif

TESTCPP=testcpp_kiss_$(DATATYPE)

TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)

ifeq  "$(DATATYPE)" "int16_t"
 TYPEFLAGS=-DFIXED_POINT=16
endif

ifeq  "$(DATATYPE)" "int32_t"
 TYPEFLAGS=-DFIXED_POINT=32
endif

ifeq  "$(DATATYPE)" "simd"
 TYPEFLAGS=-DUSE_SIMD=1 -msse
endif

SRCFILES= $(SRCFOLDER)/kiss_fft.c $(SRCFOLDER)/tools/kiss_fftnd.c $(SRCFOLDER)/tools/kiss_fftr.c $(SRCFOLDER)/test/pstats.c $(SRCFOLDER)/tools/kfc.c $(SRCFOLDER)/tools/kiss_fftndr.c

$(TESTCPP): testcpp.cc $(SRCFILES)
	$(CC) -o $@ $(CFLAGS) $(TYPEFLAGS)  $+ -lm
