CFLAGS = -O3 -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -std=gnu99 -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wdeclaration-after-statement -Wbad-function-cast -Wnested-externs
FUNCTIONS_O = solve.o matrix.o check.o
HEADERS = solve.h returns.h matrix.h check.h
MATH = -lm

%.o: %.c $(HEADERS)
	gcc -c $(CFLAGS) $<
%.out: %.o $(FUNCTIONS_O)
	gcc $(CFLAGS) $^ -o $@ $(MATH) 
all: a.out

archive:
	zip Vatchenko_MY.zip Makefile *.c *.h 

rmarchive:
	rm Vatchenko_MY.zip

clean:
	rm -f *.o *.out all
