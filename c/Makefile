CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99 -fPIC
LDFLAGS = -lm

LIB_SRC = universalfft.c
LIB_OBJ = universalfft.o

.PHONY: all lib demo test clean

all: lib demo

lib: libuniversalfft.so

libuniversalfft.so: $(LIB_OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

$(LIB_OBJ): $(LIB_SRC) universalfft.h
	$(CC) $(CFLAGS) -c -o $@ $<

demo: demo.c $(LIB_OBJ)
	$(CC) $(CFLAGS) -o ufft_demo demo.c $(LIB_OBJ) $(LDFLAGS)

test: demo
	./ufft_demo

clean:
	rm -f *.o *.so ufft_demo *.csv
