
CC=gcc -std=gnu9x
CXX=g++
LDLIBS=-lm -lrt

PREFIX=stinger
OBJ_PREFIX=stinger

STINGER_SRC=$(PREFIX)/stinger.c $(PREFIX)/stinger-utils.c $(PREFIX)/timer.c $(PREFIX)/xmalloc.c $(PREFIX)/x86-full-empty.c


CFLAGS+=-g -O2 -I$(PREFIX) -fopenmp -g -w


SRC = streaming_insertions.c streaming_deletions.c operation.c streamingbc_aux.c

SRCSBC= sbcMain.c $(SRC) $(STINGER_SRC) dsUtils.c


main: $(SRCSBC)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)


.PHONY:	clean
clean:
	rm -f main \
		$(MAINPL) $(GENSTREAMSPL)
