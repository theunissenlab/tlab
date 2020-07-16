LEVMARLIB=-llevmar -L/Users/hsoula/Devel/levmar-2.6/ -llapack -lcblas -latlas -lf2c -L /opt/local/lib/
GSLLIB=-lgsl -lcblas -lm

all: smSLc.out

smSL.out: sm.h sm.o smSL.c wave.o
	gcc -o smSL.out smSL.c sm.o wave.o $(GSLLIB) -lsndfile

smSLc.out: sm.h sm.o smSLc.c wave.o
	gcc -o smSLc.out smSLc.c sm.o wave.o $(GSLLIB) -lsndfile

sm.o: sm.c sm.h
	gcc -g -c sm.c 
wave.o: wave.c wave.h
	gcc -g -c wave.c

clean:
	rm -f *.out *~ *.o