test.out: test.c utils.o ga.o
	gcc -g -o test.out test.c utils.o ga.o 

ga.o: ga.h ga.c 
	gcc -g -c ga.c 
utils.o: utils.h utils.c 
	gcc -g -c utils.c

clean:
	rm *.o *.out