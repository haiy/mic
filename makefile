mic:mic.o
	gcc -o mic mic.c -g -lm
clean:
	rm mic.o
