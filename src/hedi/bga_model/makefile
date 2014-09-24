CC=gcc
macCFLAGS= -fPIC
files = src/bga.o src/fit.o src/spectrum.o src/optim.o src/utils.o
mac_include=-I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/Library/Python/2.7/site-packages/numpy/core/include
mac_libs= -lbga -lgsl -lcblas -framework Python -lsndfile
mac_libs_dir=-L./


sftor_include=-I/usr/include/python2.7/ -I/auto/fhome/hsoula/include
sftor_libs= -lbga -lgsl -lgslcblas -lpython2.7
sftor_libs_dir= -I/usr/include/python2.7/ -I/auto/fhome/hsoula/include

all: libbga.a find_fund.out pyBGA.so

find_fund.out: find_fund.c libbga.a
	$(CC) -o find_fund.out find_fund.c $(mac_libs)  $(mac_libs_dir)

pyBGA.c: pyBGA.pyx
	cython pyBGA.pyx

pyBGA.o: pyBGA.c
	$(CC) $(macCFLAGS) $(mac_include) -c -o pyBGA.o pyBGA.c 

pyBGA.so: pyBGA.o libbga.a	
	gcc -shared pyBGA.o -o pyBGA.so $(mac_libs)  $(mac_libs_dir)

$(filter %.o,$(files)): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

libbga.a: $(files)
	ar rus libbga.a $(files)

clean:
	rm *.o 