####################################################
# C. Baldenegro-interface for resonances0even_aazz coupling 
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE 
SOFLAGS=-fPIC
GSL=`gsl-config --cflags --libs`

all: resonances0evenoh_sqme_zz.a libresonances0evenoh_sqme_zz.so

objects=functions.o helicity_amplitudes.o sqme.o
f_0_objects=$(objects:%=src/%)

resonances0evenoh_sqme_zz.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

libresonances0evenoh_sqme_zz.so:$(f_0_objects)
	$(CC) $(SOFLAGS) -shared $(f_0_objects) -o $@ $(GSL)

.cpp.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $(GSL)

clean:
	rm -rf src/*.o test/*.o *.a *.so testcall testcall2

