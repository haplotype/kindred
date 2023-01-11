#***************************************************************************
#       user configuration
# set mpi = yes for the parallel version
# set debug = yes to debug
DEBUG ?= no
# end user configuration;
#***************************************************************************

#SHELL = /bin/sh

CC = g++

ifeq ($(strip $(DEBUG)), yes)
    CFLAGS += -ggdb -DDEBUG -D THPOOL_DEBUG
else
    CFLAGS += -O3
endif

CFLAGS += -Wall -m64 -fno-math-errno  -ffast-math

HTSINC?=/usr/local/include/htslib
HTSLIB?=/usr/local/lib
CFLAGS +=-I$(HTSINC) -g
#LDFLAGS =-L$(HTSLIB)  -lm -lhts # -lz -llzma -lbz2
LDFLAGS = -lm -lhts # -lz -llzma -lbz2
#LDFLAGS += /usr/local/lib/libhts.a
#LDFLAGS += /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a
LDFLAGS += -lpthread 

OBJS = kindred.o lsfit.o thpool.o

all:	kindred 

.cpp.o : ;  $(CC) $(CFLAGS) -c $<

$(OBJS): %.o: %.cpp

kindred: kindred.o lsfit.o thpool.o; $(CC) $(OBJS) $(LDFLAGS) -o kindred
#test: test_vcf.o;  $(CC) test_vcf.o $(LDFLAGS) -o test
#kde: test_kde.o kde.o;  $(CC) test_kde.o kde.o $(LDFLAGS) -o kde  
#trd: test_thread.o;  $(CC) test_thread.o $(LDFLAGS) -o tmt  

kindred.o: kindred.cpp 
thpool.o: thpool.cpp thpool.h
lsfit.o: lsfit.cpp lsfit.h 
#test_kde.o: test_kde.cpp kde.cpp kde.h
#test_thread.o: test_thread.cpp


clean:
	rm -f kindred kindred.o lsfit.o thpool.o
#	rem -f kde.o test test_kde.o test_vcf.o kde test_kde.o 

