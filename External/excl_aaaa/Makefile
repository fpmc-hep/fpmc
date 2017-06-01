####################################################
# M. Saimpert-interface for sm_aaaa coupling 
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE 
SOFLAGS=-fPIC
GSL=`gsl-config --cflags --libs`
INS=-Icommons/

all: excl_aaaa.a libexcl_aaaa.so

commons=commons/functions.o commons/helicity_amplitudes.o
commons_objects=$(commons:%=%)
objects=sm_sqme_aaaa/sqme.o \
	bsmf_sqme_aaaa/sqme.o \
	bsmv_sqme_aaaa/sqme.o \
	eft_sqme_aaaa/sqme.o \
	resonances0even_sqme_aaaa/sqme.o \
	resonances0even_sqme_az/sqme.o \
	resonances0even_sqme_ww/sqme.o \
	resonances0even_sqme_zz/sqme.o \
	resonances2_sqme_aaaa/sqme.o
f_0_objects=$(objects:%=%)

excl_aaaa.a:$(f_0_objects) $(commons_objects)
	echo $(f_0_objects)
	ar -r $@ $(f_0_objects) $(commons_objects)

libexcl_aaaa.so:$(f_0_objects) $(commons_objects)
	$(CC) $(SOFLAGS) -shared $(f_0_objects) $(commons_objects) -o $@ $(GSL) $(COMMON)

.cpp.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $(GSL) $(COMMON)

clean:
	rm -rf src/*.o test/*.o *.a *.so

