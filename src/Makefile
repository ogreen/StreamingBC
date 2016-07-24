
CC=gcc -std=gnu9x
CXX=g++
LDLIBS=-lm -lrt

PREFIX=stinger

STINGER_SRC=$(PREFIX)/stinger.c $(PREFIX)/stinger-utils.c $(PREFIX)/timer.c $(PREFIX)/xmalloc.c $(PREFIX)/x86-full-empty.c

CFLAGS+=-g -I$(PREFIX) -fopenmp -g -w

SRC = streamingbc.c streamingbc_init.c streamingbc_internal.c streamingbc_aux.c

SRCSBC= sbcMain.c $(SRC) $(STINGER_SRC)

main: $(SRCSBC)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

.PHONY:	clean
clean:
	rm -f main
