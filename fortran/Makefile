FC      = gfortran
FFLAGS  = -O2 -Wall -std=f2008 -fimplicit-none
MOD_OBJ = universalfft.o

.PHONY: all demo test clean

all: demo

universalfft.o: universalfft.f90
	$(FC) $(FFLAGS) -c -o $@ $<

demo: demo.f90 $(MOD_OBJ)
	$(FC) $(FFLAGS) -o ufft_demo_f demo.f90 $(MOD_OBJ)

test: demo
	mkdir -p ../tests/data
	./ufft_demo_f ../tests/data

clean:
	rm -f *.o *.mod ufft_demo_f *.csv
