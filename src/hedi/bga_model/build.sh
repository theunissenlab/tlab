cython pyBGA.pyx ; 
gcc -O3 -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/Library/Python/2.7/site-packages/numpy/core/include -c -o pyBGA.o pyBGA.c 
make
gcc -shared pyBGA.o -o pyBGA.so bga.o spectrum.o -framework Python -lgsl -lcblas -lbga -L./	