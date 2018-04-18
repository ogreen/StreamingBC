
CC=gcc -std=gnu9x
CXX=g++ -std=gnu9x
LDLIBS=-lm -lrt

PREFIX=stinger

STINGER_SRC=$(PREFIX)/stinger.c $(PREFIX)/stinger-utils.c $(PREFIX)/timer.c $(PREFIX)/xmalloc.c $(PREFIX)/x86-full-empty.c

CFLAGS+=-g -I$(PREFIX) -O3 -fopenmp -g -w
CFLAGS+=-I ${BOOST_ROOT}/include -I include

SRC = streamingbc.c streamingbc_init.c streamingbc_internal.c streamingbc_internal_fine.c streamingbc_aux.c

SRCSBC= sbcMain.c sbcMain_utils.c $(SRC) $(STINGER_SRC)

main: $(SRCSBC)
	$(CXX) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

.PHONY:	clean
clean:
	rm -f main
