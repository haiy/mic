CFLAGS=-g -lm -Wall
mic_set.exe:main.o mic.o mic_matrix.o  para_chk.o read_content.o
	cc -o mic_set main.o  mic.o mic_matrix.o  para_chk.o read_content.o $(CFLAGS)
main.o:main.c main.h
	cc -c main.c $(CFLAGS)
mic.o:mic.c mic.h
	cc -c mic.c $(CFLAGS)
mic_matrix.o:mic_matrix.c mic_matrix.h
	cc -c mic_matrix.c $(CFLAGS)
para_chk.o:para_chk.c para_chk.h
	cc -c para_chk.c $(CFLAGS)
read_content.o:read_content.c read_content.h
	cc -c read_content.c $(CFLAGS)
clean:
	rm main *.o


