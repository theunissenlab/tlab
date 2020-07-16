LEVMARLIB=-llevmar -L/Users/hsoula/Devel/levmar-2.6/ -llapack -lcblas -latlas -lf2c -L /opt/local/lib/
GSLLIB=-lgsl -lcblas -lm

all: smMM.out smGA.out

smGA.out: smGA.c sm.h sm.o ../ga/ga.o ../ga/ga.h ../ga/utils.h ../ga/utils.o
	gcc -g -o smGA.out smGA.c sm.o ../ga/ga.o ../ga/utils.o -lsndfile $(GSLLIB)

smMM.out: smMM.c sm.o sm.h
	gcc -g -o smMM.out smMM.c sm.o -lsndfile $(GSLLIB)

sm.o: sm.c sm.h
	gcc -g -c sm.c 

clean:
	rm -f *.out *~ *.o